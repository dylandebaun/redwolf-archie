#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx[length(idx)] == length(args)) {
    stop(sprintf("Missing value for %s", flag))
  }
  args[idx[length(idx)] + 1]
}

parse_bool <- function(x, default = TRUE) {
  if (is.null(x)) return(default)
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

# --- read the column contract from sim_meta.json (no more ms.sh parsing) ---
read_sim_meta <- function(path = "sim_meta.json") {
  if (!file.exists(path)) {
    stop(sprintf("sim_meta.json not found at '%s'. Run emit_sim_meta.py first.", path))
  }
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    return(jsonlite::fromJSON(path))
  }
  txt <- paste(readLines(path), collapse = " ")
  keys <- regmatches(txt, gregexpr('"([a-z_]+)"\\s*:\\s*-?[0-9]+', txt))[[1]]
  out <- list()
  for (kv in keys) {
    k <- sub('"([a-z_]+)".*', "\\1", kv)
    v <- as.integer(sub('.*:\\s*(-?[0-9]+).*', "\\1", kv))
    out[[k]] <- v
  }
  out
}

train_path <- get_arg("--train", "../simulations/training_data.txt")
test_path <- get_arg("--test", "../simulations/test_data.txt")
meta_path <- get_arg("--sim-meta", Sys.getenv("SIM_META", "sim_meta.json"))
model_rdata <- get_arg("--model-rdata", "trained_model.Rdata")
out_prefix <- get_arg("--out-prefix", "threshold_calibration")
use_existing_model <- parse_bool(get_arg("--use-existing-model", "true"), default = TRUE)
include_inbetween <- parse_bool(get_arg("--include-inbetween", "false"), default = FALSE)
fpr_targets_str <- get_arg("--target-fpr", "0.01,0.05,0.10")
fpr_targets <- as.numeric(strsplit(fpr_targets_str, ",")[[1]])
fpr_targets <- fpr_targets[is.finite(fpr_targets) & fpr_targets >= 0 & fpr_targets <= 1]
if (length(fpr_targets) == 0) fpr_targets <- c(0.01, 0.05, 0.10)

message("Loading training/test data...")
train <- data.frame(fread(train_path, fill = TRUE))
test <- data.frame(fread(test_path, fill = TRUE))

# --- column layout straight from the contract ---
meta <- read_sim_meta(meta_path)
archaic_col <- meta$archaic_col
not_archaic_col <- meta$not_archaic_col
inbetween_col <- meta$inbetween_col
feature_cols <- 1:meta$feature_count
# optional (simulations/simulate_archie.py sets it): leave out IFS bin N+1,
# "derived in all N focal haps" -- no between-haplotype ancestry signal, and the
# outgroup-polarized sims can't reproduce its real-data count.
if (isTRUE(as.integer(meta$drop_all_derived_bin) == 1L)) {
  feature_cols <- setdiff(feature_cols, meta$n_focal_haps + 1)
}

message(sprintf(
  "From sim_meta.json: n_focal=%d, feature_count=%d, label cols (arch/not/inbtwn)=%d/%d/%d.",
  meta$n_focal_haps, meta$feature_count, archaic_col, not_archaic_col, inbetween_col
))

if (max(c(archaic_col, not_archaic_col, inbetween_col)) > ncol(test)) {
  stop(sprintf(
    "Test data has %d columns but label columns need up to %d",
    ncol(test), max(c(archaic_col, not_archaic_col, inbetween_col))
  ))
}

train_lr <- train[, c(feature_cols, archaic_col)]
names(train_lr)[ncol(train_lr)] <- "label"
train_lr$label <- as.integer(round(train_lr$label))

model <- NULL
model_source <- "trained"
if (use_existing_model && file.exists(model_rdata)) {
  message(sprintf("Loading model from %s", model_rdata))
  env <- new.env(parent = emptyenv())
  load(model_rdata, envir = env)
  if (exists("model", envir = env, inherits = FALSE) && inherits(get("model", envir = env), "glm")) {
    model <- get("model", envir = env)
    model_source <- "loaded"
  } else {
    objs <- ls(env)
    for (nm in objs) {
      obj <- get(nm, envir = env)
      if (inherits(obj, "glm")) {
        model <- obj
        model_source <- sprintf("loaded:%s", nm)
        break
      }
    }
  }
}

if (is.null(model)) {
  message("No usable existing model found; training glm on training_data...")
  model <- glm(label ~ ., family = binomial(link = "logit"), data = train_lr)
  model_source <- "trained"
}

if (!all(names(train_lr)[1:(ncol(train_lr)-1)] %in% names(test))) {
  names(test) <- paste0("V", seq_len(ncol(test)))
}

test_x <- test[, feature_cols, drop = FALSE]
truth <- as.integer(round(test[[archaic_col]]))
not_arch <- as.integer(round(test[[not_archaic_col]]))
inbetween <- as.integer(round(test[[inbetween_col]]))

ok <- complete.cases(test_x) & is.finite(truth) & truth %in% c(0, 1)
if (!include_inbetween) {
  ok <- ok & (inbetween == 0)
}

if (sum(ok) == 0) stop("No valid test rows after filtering.")

if (sum(!ok) > 0) {
  message(sprintf("Dropping %d test rows (malformed/missing/filtered).", sum(!ok)))
}

probs <- as.numeric(plogis(predict.glm(model, newdata = test_x[ok, , drop = FALSE])))
y <- truth[ok]

if (length(unique(y)) < 2) {
  stop("Filtered test set does not contain both positive and negative classes.")
}

thresholds <- seq(0, 1, by = 0.001)

metric_rows <- lapply(thresholds, function(t) {
  pred <- as.integer(probs >= t)
  tp <- sum(pred == 1 & y == 1)
  fp <- sum(pred == 1 & y == 0)
  tn <- sum(pred == 0 & y == 0)
  fn <- sum(pred == 0 & y == 1)

  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  recall <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  fpr <- if ((fp + tn) > 0) fp / (fp + tn) else NA_real_
  f1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA_real_
  }
  bal_acc <- if (!is.na(recall) && !is.na(specificity)) (recall + specificity) / 2 else NA_real_
  youden_j <- if (!is.na(recall) && !is.na(fpr)) recall - fpr else NA_real_

  data.frame(
    threshold = t, tp = tp, fp = fp, tn = tn, fn = fn,
    precision = precision, recall = recall, specificity = specificity,
    fpr = fpr, f1 = f1, balanced_accuracy = bal_acc, youden_j = youden_j,
    stringsAsFactors = FALSE
  )
})

metrics <- do.call(rbind, metric_rows)

best_f1 <- metrics[which.max(replace(metrics$f1, is.na(metrics$f1), -Inf)), , drop = FALSE]
best_j <- metrics[which.max(replace(metrics$youden_j, is.na(metrics$youden_j), -Inf)), , drop = FALSE]

pick_fpr_row <- function(target) {
  m <- metrics[!is.na(metrics$fpr) & metrics$fpr <= target, , drop = FALSE]
  if (nrow(m) == 0) {
    return(data.frame(target_fpr = target, threshold = NA_real_, fpr = NA_real_, recall = NA_real_, precision = NA_real_, f1 = NA_real_, stringsAsFactors = FALSE))
  }
  max_recall <- max(m$recall, na.rm = TRUE)
  recall_candidates <- m[!is.na(m$recall) & m$recall == max_recall, , drop = FALSE]
  if (nrow(recall_candidates) == 0) recall_candidates <- m
  max_f1 <- max(recall_candidates$f1, na.rm = TRUE)
  f1_candidates <- recall_candidates[!is.na(recall_candidates$f1) & recall_candidates$f1 == max_f1, , drop = FALSE]
  if (nrow(f1_candidates) == 0) f1_candidates <- recall_candidates
  top <- f1_candidates[order(f1_candidates$threshold, decreasing = TRUE), , drop = FALSE][1, , drop = FALSE]
  data.frame(
    target_fpr = target, threshold = top$threshold, fpr = top$fpr,
    recall = top$recall, precision = top$precision, f1 = top$f1,
    stringsAsFactors = FALSE
  )
}

fpr_recs <- do.call(rbind, lapply(fpr_targets, pick_fpr_row))

n_pos <- sum(y == 1)
n_neg <- sum(y == 0)
ranks <- rank(probs, ties.method = "average")
auc <- (sum(ranks[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)

summary_tbl <- data.frame(
  model_source = model_source,
  train_rows = nrow(train),
  test_rows_total = nrow(test),
  test_rows_used = sum(ok),
  include_inbetween = include_inbetween,
  n_focal_haps = meta$n_focal_haps,
  feature_count = length(feature_cols),
  archaic_label_col = archaic_col,
  not_archaic_col = not_archaic_col,
  inbetween_col = inbetween_col,
  auc = auc,
  best_f1_threshold = best_f1$threshold,
  best_f1_value = best_f1$f1,
  best_f1_precision = best_f1$precision,
  best_f1_recall = best_f1$recall,
  best_youden_threshold = best_j$threshold,
  best_youden_j = best_j$youden_j,
  stringsAsFactors = FALSE
)

metrics_out <- paste0(out_prefix, ".metrics.tsv")
summary_out <- paste0(out_prefix, ".summary.tsv")
fpr_out <- paste0(out_prefix, ".recommended_by_fpr.tsv")

fwrite(metrics, metrics_out, sep = "\t")
fwrite(summary_tbl, summary_out, sep = "\t")
fwrite(fpr_recs, fpr_out, sep = "\t")

cat("\n=== Threshold Calibration Summary ===\n")
cat(sprintf("Model source: %s\n", model_source))
cat(sprintf("Training rows: %d\n", nrow(train)))
cat(sprintf("Test rows used: %d / %d\n", sum(ok), nrow(test)))
cat(sprintf("Feature columns: %d\n", length(feature_cols)))
cat(sprintf("Label columns (archaic/not/inbetween): %d/%d/%d\n", archaic_col, not_archaic_col, inbetween_col))
cat(sprintf("AUC: %.4f\n", auc))
cat(sprintf("Best F1 threshold: %.3f (F1=%.4f, precision=%.4f, recall=%.4f)\n", best_f1$threshold, best_f1$f1, best_f1$precision, best_f1$recall))
cat(sprintf("Best Youden threshold: %.3f (J=%.4f, recall=%.4f, fpr=%.4f)\n", best_j$threshold, best_j$youden_j, best_j$recall, best_j$fpr))
cat("\nRecommended thresholds by target FPR:\n")
print(fpr_recs, row.names = FALSE)
cat("\nWrote files:\n")
cat(sprintf("  %s\n", normalizePath(metrics_out, mustWork = FALSE)))
cat(sprintf("  %s\n", normalizePath(summary_out, mustWork = FALSE)))
cat(sprintf("  %s\n", normalizePath(fpr_out, mustWork = FALSE)))
