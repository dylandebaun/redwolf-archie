#!/usr/bin/env Rscript
# apply_model.R — score real ArchIE window features with the trained model and
# call windows at the calibrated thresholds.
#
#   Rscript scripts/apply_model.R \
#     --stats-glob 'results/x/data/window_stats/chr*.stats.tsv.gz' \
#     --model-rdata results/x/model/trained_model.Rdata \
#     --thresholds  results/x/model/thresholds.recommended_by_fpr.tsv \
#     --sim-meta    results/x/model/sim_meta.json \
#     --out-prefix  results/x/calls/x
#
# The positive class is the simulation label (calc_stats_ms_rw.py): red-wolf-major
# ancestry, i.e. >70% of the window NOT from the coyote/donor lineage.
#
# Writes:
#   <prefix>.windows.tsv.gz      one row per haplotype x window: sample, hap, chrom,
#                                start, end, prob, call_fpr<target> (0/1) per threshold
#   <prefix>.per_sample.tsv      fraction of windows called, per sample and threshold
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  args[idx[length(idx)] + 1]
}
stats_glob <- get_arg("--stats-glob")
model_rdata <- get_arg("--model-rdata")
thr_path <- get_arg("--thresholds")
meta_path <- get_arg("--sim-meta")
out_prefix <- get_arg("--out-prefix")
if (any(sapply(list(stats_glob, model_rdata, thr_path, meta_path, out_prefix), is.null))) {
  stop("required: --stats-glob --model-rdata --thresholds --sim-meta --out-prefix")
}

meta <- if (requireNamespace("jsonlite", quietly = TRUE)) jsonlite::fromJSON(meta_path) else {
  txt <- paste(readLines(meta_path), collapse = " ")
  list(feature_count = as.integer(sub('.*"feature_count"\\s*:\\s*([0-9]+).*', "\\1", txt)))
}
n_feat <- as.integer(meta$feature_count)

env <- new.env(parent = emptyenv())
load(model_rdata, envir = env)
model <- get("model", envir = env)
thr <- fread(thr_path)[!is.na(threshold)]

files <- Sys.glob(stats_glob)
if (length(files) == 0) stop("no files match ", stats_glob)
# zcat instead of fread("*.gz"): fread needs R.utils for gz, which the archie env lacks
real <- rbindlist(lapply(files, function(f) fread(cmd = paste("zcat", shQuote(f)), header = FALSE)))
if (ncol(real) != n_feat + 4) {
  stop(sprintf("stats files have %d columns; model expects %d features + 4 id columns", ncol(real), n_feat))
}
setnames(real, (n_feat + 1):(n_feat + 4), c("sample", "chrom", "start", "end"))
# two consecutive rows per sample per window = the two haplotypes
real[, hap := seq_len(.N), by = .(sample, chrom, start)]

real[, prob := plogis(predict.glm(model, newdata = real[, 1:n_feat]))]
out <- real[, .(sample, hap, chrom, start, end, prob)]
for (i in seq_len(nrow(thr))) {
  col <- sprintf("call_fpr%g", thr$target_fpr[i])
  out[, (col) := as.integer(prob >= thr$threshold[i])]
}

dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
win_out <- paste0(out_prefix, ".windows.tsv")
fwrite(out, win_out, sep = "\t")
system2("gzip", c("-f", shQuote(win_out)))

call_cols <- grep("^call_fpr", names(out), value = TRUE)
# windows with no usable sites have NaN features -> NaN prob; leave them out of
# the per-sample summary (n_scored = windows that got a probability)
per_sample <- out[, c(list(n_windows = .N, n_scored = sum(!is.na(prob)), mean_prob = mean(prob, na.rm = TRUE)),
                      lapply(.SD, mean, na.rm = TRUE)),
                  by = sample, .SDcols = call_cols]
fwrite(per_sample, paste0(out_prefix, ".per_sample.tsv"), sep = "\t")

cat(sprintf("Scored %d haplotype-windows from %d files\n", nrow(out), length(files)))
cat("Thresholds:\n"); print(thr[, .(target_fpr, threshold)], row.names = FALSE)
cat("Fraction of haplotype-windows called, per sample:\n"); print(per_sample, digits = 3)
cat(sprintf("Wrote %s.gz and %s.per_sample.tsv\n", win_out, out_prefix))
