# =====================================================================
# recover_coverage.R
# -----------------------------------------------------------------
# Post-hoc recomputation of two-sample interval coverage for the HAR
# and bootstrap layers of the Camino-C run, fixing the sign-convention
# bug present in HTDV <= 0.1.0 (HAR/boot reported delta = mean(x1) -
# mean(x2) while the Stan layer reported delta = alpha2 - alpha1; the
# coverage check used +true_delta against a HAR/boot CI centred at
# -true_delta, which produced spuriously low coverage for delta > 0).
#
# Inputs (read from this directory):
#   simstudy_raw.rds        full per-(cell, rep, layer) data frame
# Outputs:
#   simstudy_summary_v2.csv recomputed aggregated summary
#   simstudy_benchmarks_v2.rds  recomputed B1/B2/B3 lists
#
# Run from R with:
#   source("recover_coverage.R")
# or from the shell with:
#   Rscript recover_coverage.R
# =====================================================================

suppressPackageStartupMessages({
  library(HTDV)
})

.resolve_script_dir <- function() {
  frame_ofile <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
  if (!is.null(frame_ofile) && nzchar(frame_ofile))
    return(dirname(normalizePath(frame_ofile)))
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m) == 1L)
    return(dirname(normalizePath(sub("^--file=", "", m))))
  getwd()
}
out_dir   <- .resolve_script_dir()
raw_file  <- file.path(out_dir, "simstudy_raw.rds")
summ_file <- file.path(out_dir, "simstudy_summary_v2.csv")
bench_file <- file.path(out_dir, "simstudy_benchmarks_v2.rds")

if (!file.exists(raw_file)) {
  stop("simstudy_raw.rds not found in ", out_dir)
}

message("Reading raw results from: ", raw_file)
res <- readRDS(raw_file)

message(sprintf("Total rows: %d", nrow(res)))
message(sprintf("Layers: %s", paste(unique(res$layer), collapse = ", ")))

# ---- Sign correction for HAR and bootstrap rows ----------------------
# The Bayesian layer is unaffected (its delta is alpha2 - alpha1 by
# construction). The reject decision for HAR/boot is sign-symmetric
# (|z| > z_crit; or 0 outside [q_lo, q_hi]) so it stays as recorded.
# What we recompute is the CI orientation and the coverage flag.

flip <- res$layer %in% c("har", "boot")
if (any(flip)) {
  message(sprintf("Flipping sign on %d HAR/boot rows.", sum(flip)))
  est_old   <- res$estimate[flip]
  ci_lo_old <- res$ci_lo[flip]
  ci_hi_old <- res$ci_hi[flip]
  res$estimate[flip] <- -est_old
  res$ci_lo[flip]    <- -ci_hi_old
  res$ci_hi[flip]    <- -ci_lo_old
}

true_delta <- res$delta * res$sigma_inf
ci_lo <- res$ci_lo
ci_hi <- res$ci_hi

res$covered <- ifelse(is.finite(ci_lo) & is.finite(ci_hi),
                      (true_delta >= ci_lo) & (true_delta <= ci_hi),
                      NA)
res$ci_length <- ifelse(is.finite(ci_lo) & is.finite(ci_hi),
                        ci_hi - ci_lo, NA_real_)

# ---- Re-aggregate ----------------------------------------------------
message("Re-aggregating with htdv_simstudy_summary().")
summ <- htdv_simstudy_summary(res, alpha = 0.05)

utils::write.csv(summ, summ_file, row.names = FALSE)
message("Wrote: ", summ_file)

# ---- Re-run pre-registered benchmarks --------------------------------
check_B1 <- function(summ) {
  size_cells <- summ[summ$delta == 0, ]
  size_cells$size_ok <- abs(size_cells$reject_rate - 0.05) <= 0.02
  ignore <- size_cells$tail == "t2_1" & size_cells$layer %in% c("har", "boot")
  size_cells$size_ok[ignore] <- NA
  keys <- with(size_cells, paste(n, phi, tail, imb, sep = "|"))
  ok <- vapply(split(size_cells$size_ok, keys), function(v)
                 sum(v, na.rm = TRUE) >= min(3L, sum(!is.na(v))),
               logical(1L))
  list(pass_rate = mean(ok), n_cells = length(ok),
       failures = names(ok)[!ok])
}
check_B2 <- function(summ) {
  keys <- with(summ, paste(n, phi, tail, imb, layer, sep = "|"))
  sp <- split(summ[, c("delta", "reject_rate")], keys)
  mono <- vapply(sp, function(d) {
    d <- d[order(d$delta), ]
    all(diff(d$reject_rate) >= -0.05)
  }, logical(1L))
  list(pass_rate = mean(mono), n_cells = length(mono),
       failures = names(mono)[!mono])
}
check_B3 <- function(summ) {
  b <- summ[summ$layer == "bayes", ]
  ok <- b$coverage >= 0.93
  list(pass_rate = mean(ok, na.rm = TRUE), n_cells = nrow(b),
       failures = with(b[!ok & !is.na(ok), ],
                       paste(n, phi, tail, imb, delta, sep = "|")))
}

message("\n=== Pre-registered benchmark checks (v2, sign-corrected) ===")
B1 <- check_B1(summ); B2 <- check_B2(summ); B3 <- check_B3(summ)
message(sprintf("B1 (size control):       %.3f over %d cells",
                B1$pass_rate, B1$n_cells))
message(sprintf("B2 (power monotonic):    %.3f over %d cell-layers",
                B2$pass_rate, B2$n_cells))
message(sprintf("B3 (Bayes coverage>=.93):%.3f over %d cells",
                B3$pass_rate, B3$n_cells))

saveRDS(list(B1 = B1, B2 = B2, B3 = B3), bench_file)
message("Wrote: ", bench_file)

# ---- Layer-aware coverage diagnostics --------------------------------
message("\n=== Coverage by layer (sign-corrected, all delta) ===")
by_layer <- split(summ, summ$layer)
for (lyr in names(by_layer)) {
  d <- by_layer[[lyr]]
  message(sprintf("%-6s n=%d mean=%.4f sd=%.4f min=%.3f max=%.3f",
                  lyr, nrow(d), mean(d$coverage, na.rm = TRUE),
                  sd(d$coverage, na.rm = TRUE),
                  min(d$coverage, na.rm = TRUE),
                  max(d$coverage, na.rm = TRUE)))
}

message("\n=== Difficult cells (Bayesian diag_pass < 0.7) ===")
flag_difficult <- function(summ, threshold = 0.7) {
  bayes <- summ[summ$layer == "bayes", , drop = FALSE]
  flagged <- bayes[is.finite(bayes$diag_pass_rate) &
                   bayes$diag_pass_rate < threshold, , drop = FALSE]
  flagged[order(flagged$diag_pass_rate),
          c("n", "phi", "tail", "imb", "delta", "layer", "diag_pass_rate"),
          drop = FALSE]
}
warns <- if (exists("htdv_simstudy_warnings", mode = "function")) {
  htdv_simstudy_warnings(summ, threshold = 0.7)
} else {
  flag_difficult(summ, threshold = 0.7)
}
if (nrow(warns) == 0L) {
  message("None.")
} else {
  message(sprintf("Total: %d cell-layer rows", nrow(warns)))
  print(utils::head(warns, 20))
  if (nrow(warns) > 20L) message(sprintf("... and %d more.", nrow(warns) - 20L))
}

message("\nDone. Paste simstudy_summary_v2.csv back to update Section 12-bis.B.")
