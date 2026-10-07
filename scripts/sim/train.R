library(data.table)

# --- read the column contract from sim_meta.json (no more ms.sh parsing) ---
read_sim_meta <- function(path = "sim_meta.json") {
  if (!file.exists(path)) {
    stop(sprintf("sim_meta.json not found at '%s'. Run emit_sim_meta.py first.", path))
  }
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    return(jsonlite::fromJSON(path))
  }
  # dependency-free fallback for the flat JSON we emit
  txt <- paste(readLines(path), collapse = " ")
  out <- list()
  for (kv in regmatches(txt, gregexpr('"([a-z_]+)"\\s*:\\s*-?[0-9]+', txt))[[1]]) {
    k <- sub('"([a-z_]+)".*', "\\1", kv)
    out[[k]] <- as.integer(sub('.*:\\s*(-?[0-9]+).*', "\\1", kv))
  }
  for (kv in regmatches(txt, gregexpr('"([a-z_]+)"\\s*:\\s*"[^"]+"', txt))[[1]]) {
    k <- sub('"([a-z_]+)".*', "\\1", kv)
    out[[k]] <- sub('.*:\\s*"([^"]+)".*', "\\1", kv)
  }
  out
}

meta <- read_sim_meta(Sys.getenv("SIM_META", "sim_meta.json"))

ARCHAIC_LABEL_COL <- meta$archaic_col      # = 2N + 9
feature_cols      <- 1:meta$feature_count  # = 1 .. 2N + 8
# optional (simulations/simulate_archie.py sets it): leave out IFS bin N+1,
# "derived in all N focal haps" -- no between-haplotype ancestry signal, and the
# outgroup-polarized sims can't reproduce its real-data count.
if (isTRUE(as.integer(meta$drop_all_derived_bin) == 1L)) {
  feature_cols <- setdiff(feature_cols, meta$n_focal_haps + 1)
}

message(sprintf(
  "From sim_meta.json: n_focal=%d, feature_count=%d, archaic label col=%d, total_cols=%d, model_dir=%s.",
  meta$n_focal_haps, meta$feature_count, ARCHAIC_LABEL_COL, meta$total_cols, meta$model_dir
))

# --- training data ---
train_path <- file.path(meta$model_dir, "training_data.txt")
train <- data.frame(fread(train_path))

# sanity: the contract must match the actual file width
if (ncol(train) != meta$total_cols) {
  stop(sprintf(
    "training_data.txt has %d columns but sim_meta.json expects %d. Re-run sims or emit_sim_meta.py.",
    ncol(train), meta$total_cols
  ))
}

lr <- train[, c(feature_cols, ARCHAIC_LABEL_COL)] # features + archaic one-hot label
names(lr)[ncol(lr)] <- "label"
lr$label <- as.integer(round(lr$label))

model <- glm(label ~ ., family = binomial(link = "logit"), data = lr) # train
save.image(file.path(meta$model_dir,"trained_model.Rdata")) # reload later with load("trained_model.Rdata")


# --- predict on test data ---
test_path <- file.path(meta$model_dir, "test_data.txt")
test <- data.frame(fread(test_path, fill = TRUE))
test_x <- test[, feature_cols, drop = FALSE]
ok <- complete.cases(test_x)
if (!all(ok)) {
  warning(sprintf("Dropping %d malformed test rows before prediction.", sum(!ok)))
}
predicted <- rep(NA_real_, nrow(test_x))
predicted[ok] <- plogis(predict.glm(model, newdata = test_x[ok, , drop = FALSE]))

test$predicted <- predicted

# diagnostics keyed off the contract (prop-archaic col), not a hardcoded V72
prop_col <- meta$prop_archaic_col
hist(test$predicted[which(test[[prop_col]] == 1)])
hist(test$predicted[which(test[[prop_col]] < 0.05)])
hist(test[[prop_col]])
