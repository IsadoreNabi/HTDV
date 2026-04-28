# =====================================================================
# run_simstudy.R
# -----------------------------------------------------------------
# Executes the full-factorial Monte Carlo simulation study of
# Section 12-bis of the companion paper, using the exported function
# HTDV::htdv_simstudy().
#
# Target hardware: Linux workstation with >= 16 cores and >= 32 GB RAM.
# Expected wall-clock: ~22-28 hours with 16 cores.
#
# Outputs (written next to this script):
#   simstudy_raw.rds    -- full per-(cell, rep, layer) data frame
#   simstudy_summary.csv -- aggregated table (size, power, coverage, CI length)
#
# Run from R with:
#   source("run_simstudy.R")
# or from the shell with:
#   Rscript run_simstudy.R
# =====================================================================

suppressPackageStartupMessages({
  library(HTDV)
  library(rstan)
  library(parallel)
})

rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = "-march=native")

# ---- Configuration ---------------------------------------------------

CFG <- list(
  n_grid       = c(40L, 80L, 200L, 500L),
  phi_grid     = c(0, 0.3, 0.6, 0.85),
  tail_grid    = c("normal", "t5", "t3", "t2_1"),
  imb_grid     = c(1, 1.5, 3, 6),
  delta_grid   = c(0, 0.1, 0.25, 0.5),
  R            = 500L,
  seed         = 20260422L,
  n_cores      = 16L,
  layers       = c("har", "boot", "bayes"),
  bayes_chains = 2L,
  bayes_iter   = 600L,
  bayes_warmup = 300L,
  boot_R       = 999L,
  rope         = c(-0.1, 0.1),
  alpha        = 0.05
)

# Resolve script directory whether called via source() or Rscript.
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
inc_dir   <- file.path(out_dir, "simstudy_cells_dir")  # per-cell cache
summ_file <- file.path(out_dir, "simstudy_summary.csv")

# Logging is handled entirely by the shell:
#   nohup Rscript run_simstudy.R > sim_stdout.log 2>&1 &
# All message() output (master + mclapply workers) lands in sim_stdout.log,
# unbuffered via stderr. No R-level sink() here — that path was buffered
# and lost the worker progress lines.

message("=== htdv_simstudy run started at ", Sys.time(), " ===")
message("Grid: ",
        "n=", paste(CFG$n_grid, collapse = ","), "  ",
        "phi=", paste(CFG$phi_grid, collapse = ","), "  ",
        "tail=", paste(CFG$tail_grid, collapse = ","), "  ",
        "imb=", paste(CFG$imb_grid, collapse = ","), "  ",
        "delta=", paste(CFG$delta_grid, collapse = ","))
message("Replications per cell: ", CFG$R)
message("Layers: ", paste(CFG$layers, collapse = ", "))
message("Cores: ", CFG$n_cores)
message("Master seed: ", CFG$seed)
message("Per-cell cache directory: ", inc_dir)

# Fresh run (Camino C): wipe any prior per-cell cache so the run is one piece.
if (dir.exists(inc_dir)) {
  message("Wiping stale per-cell cache directory.")
  unlink(inc_dir, recursive = TRUE, force = TRUE)
}
dir.create(inc_dir, recursive = TRUE, showWarnings = FALSE)

t_start <- proc.time()

res <- htdv_simstudy(
  n_grid       = CFG$n_grid,
  phi_grid     = CFG$phi_grid,
  tail_grid    = CFG$tail_grid,
  imb_grid     = CFG$imb_grid,
  delta_grid   = CFG$delta_grid,
  R            = CFG$R,
  seed         = CFG$seed,
  n_cores      = CFG$n_cores,
  layers       = CFG$layers,
  bayes_chains = CFG$bayes_chains,
  bayes_iter   = CFG$bayes_iter,
  bayes_warmup = CFG$bayes_warmup,
  boot_R       = CFG$boot_R,
  rope         = CFG$rope,
  alpha        = CFG$alpha,
  out_dir      = inc_dir,
  progress     = TRUE
)

t_elapsed <- proc.time() - t_start
message(sprintf("=== htdv_simstudy finished in %.1f minutes ===",
                t_elapsed[["elapsed"]] / 60))

saveRDS(res, raw_file)
message("Wrote raw results to: ", raw_file)

summ <- htdv_simstudy_summary(res, alpha = CFG$alpha)
utils::write.csv(summ, summ_file, row.names = FALSE)
message("Wrote summary to: ", summ_file)

# ---- Pre-registered benchmark checks ---------------------------------

check_B1 <- function(summ) {
  # Empirical size |rej - 0.05| <= 0.02 on >= 3 of 4 layers per cell;
  # t2_1 exempts HAR and boot.
  size_cells <- summ[summ$delta == 0, ]
  size_cells$size_ok <- abs(size_cells$reject_rate - 0.05) <= 0.02
  ignore <- size_cells$tail == "t2_1" & size_cells$layer %in% c("har", "boot")
  size_cells$size_ok[ignore] <- NA
  keys <- with(size_cells,
               paste(n, phi, tail, imb, sep = "|"))
  ok <- vapply(split(size_cells$size_ok, keys), function(v)
                 sum(v, na.rm = TRUE) >= min(3L, sum(!is.na(v))),
               logical(1L))
  list(pass_rate = mean(ok), n_cells = length(ok), failures = names(ok)[!ok])
}

check_B2 <- function(summ) {
  # Monotone non-decreasing power in delta for every (cell, layer).
  keys <- with(summ, paste(n, phi, tail, imb, layer, sep = "|"))
  sp <- split(summ[, c("delta", "reject_rate")], keys)
  mono <- vapply(sp, function(d) {
    d <- d[order(d$delta), ]
    all(diff(d$reject_rate) >= -0.05)  # allow 5pp wiggle for MC noise
  }, logical(1L))
  list(pass_rate = mean(mono), n_cells = length(mono),
       failures = names(mono)[!mono])
}

check_B3 <- function(summ) {
  # Bayesian coverage >= 0.93 everywhere.
  b <- summ[summ$layer == "bayes", ]
  ok <- b$coverage >= 0.93
  list(pass_rate = mean(ok, na.rm = TRUE), n_cells = nrow(b),
       failures = with(b[!ok & !is.na(ok), ],
                       paste(n, phi, tail, imb, delta, sep = "|")))
}

message("\n=== Pre-registered benchmark checks ===")
B1 <- check_B1(summ)
message(sprintf("B1 (size control):      %.3f pass-rate over %d cells",
                B1$pass_rate, B1$n_cells))
B2 <- check_B2(summ)
message(sprintf("B2 (power monotonic):   %.3f pass-rate over %d cell-layers",
                B2$pass_rate, B2$n_cells))
B3 <- check_B3(summ)
message(sprintf("B3 (Bayes coverage>=0.93): %.3f pass-rate over %d cells",
                B3$pass_rate, B3$n_cells))

message("\n=== Diagnostic summary for Bayesian layer ===")
b <- res[res$layer == "bayes", ]
message(sprintf("diag_pass overall rate: %.3f",
                mean(b$diag_pass, na.rm = TRUE)))
message(sprintf("Mean Rhat_max:          %.4f",
                mean(b$rhat_max, na.rm = TRUE)))
message(sprintf("Mean min-ESS:           %.0f",
                mean(b$ess_min, na.rm = TRUE)))
message(sprintf("Total divergences:      %d",
                sum(b$divergences, na.rm = TRUE)))

saveRDS(list(B1 = B1, B2 = B2, B3 = B3),
        file.path(out_dir, "simstudy_benchmarks.rds"))
message("\nAll outputs written. You may now paste simstudy_summary.csv back.")
