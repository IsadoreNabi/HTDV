# =====================================================================
# run_benchmarks.R
# -----------------------------------------------------------------
# Empirical benchmark protocol for Section 12-ter of the companion
# paper. Three single-test datasets, all with public sources and
# pre-existing reference results from peer-reviewed literature.
#
#   E1. FRED-MD CPI inflation, post-1984 mean
#       File:    2026-03-MD.csv (McCracken-Ng monthly vintage)
#       Variable: CPIAUCSL (CPI All Urban Consumers, SA)
#       Test:    E[pi_t] over 1984-01..2026-02, where
#                pi_t = 1200 * (log CPI_t - log CPI_{t-1})
#       Reference: Stock & Watson (2007), JMCB. Post-1984 annualized
#                  monthly CPI inflation mean ~ 2.7% (with HAR-robust
#                  SE in the 0.3-0.5 range across vintages).
#
#   E2. Shiller log-CAPE mean reversion
#       File:    ie_data.xls (Shiller online data)
#       Variable: CAPE (cyclically-adjusted P/E ratio)
#       Test:    E[log CAPE_t] vs the rule-of-thumb fair value log(15).
#       Reference: Campbell & Shiller (1998). Long-run historical
#                  mean log CAPE ~ 2.85 (corresponding to CAPE ~ 17.3).
#
#   E3. US vs Canada 10-year government bond yields, mean equality
#       File:    us_canada_10y.csv (FRED multi-series download)
#       Variables: GS10 (US), IRLTLT01CAM156N (Canada)
#       Test:    E[GS10] - E[CA10] over 1990-01..2026-03.
#       Reference: long-run differential ~ -50 bps to +50 bps;
#                  no strong theoretical prior.
#
# Run from R with:
#   source("run_benchmarks.R")
# or from the shell with:
#   Rscript run_benchmarks.R
#
# Outputs (alongside this script):
#   benchmarks_E1_E2_E3.csv  one-row-per-dataset summary
#   benchmarks_full.rds      complete result object
# =====================================================================

suppressPackageStartupMessages({
  library(HTDV)
  library(rstan)
})
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(readxl))

rstan_options(auto_write = TRUE)
options(mc.cores = max(1L, parallel::detectCores() - 1L))

# ---- Resolve paths ---------------------------------------------------
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
out_dir          <- .resolve_script_dir()
fredmd_path      <- file.path(out_dir, "2026-03-MD.csv")
shiller_path     <- file.path(out_dir, "ie_data.xls")
us_canada_path   <- file.path(out_dir, "us_canada_10y.csv")
results_csv      <- file.path(out_dir, "benchmarks_E1_E2_E3.csv")
results_rds      <- file.path(out_dir, "benchmarks_full.rds")

for (p in c(fredmd_path, shiller_path, us_canada_path)) {
  if (!file.exists(p)) stop("Missing file: ", p)
}

# ---- Helpers ---------------------------------------------------------
har_ci_one_sample <- function(x, alpha = 0.05) {
  n <- length(x)
  m <- mean(x)
  lrv <- htdv_lrv(x, kernel = "qs", bandwidth = "andrews")$lrv
  se <- sqrt(lrv / n)
  z <- stats::qnorm(1 - alpha / 2)
  list(point = m, lrv = lrv, se = se,
       ci = c(m - z * se, m + z * se))
}

har_ci_two_sample <- function(x1, x2, alpha = 0.05) {
  n1 <- length(x1); n2 <- length(x2)
  m1 <- mean(x1);   m2 <- mean(x2)
  lrv1 <- htdv_lrv(x1, kernel = "qs", bandwidth = "andrews")$lrv
  lrv2 <- htdv_lrv(x2, kernel = "qs", bandwidth = "andrews")$lrv
  se <- sqrt(lrv1 / n1 + lrv2 / n2)
  est <- m2 - m1
  z <- stats::qnorm(1 - alpha / 2)
  list(point = est, lrv1 = lrv1, lrv2 = lrv2, se = se,
       ci = c(est - z * se, est + z * se))
}

bayes_one_sample <- function(x, chains = 4L, iter = 4000L,
                             warmup = 2000L, seed = 1L) {
  fit <- htdv_fit(x, model = "tac", chains = chains, iter = iter,
                  warmup = warmup, refresh = 0L, seed = seed)
  draws <- as.numeric(rstan::extract(fit$stanfit,
                                     pars = "theta",
                                     permuted = TRUE)$theta)
  diag <- htdv_diagnostics(fit)
  list(point = mean(draws),
       ci = as.numeric(stats::quantile(draws, c(0.025, 0.975))),
       diag_pass = diag$passed,
       rhat_max = max(diag$rhat, na.rm = TRUE),
       ess_min = min(c(diag$ess_bulk, diag$ess_tail), na.rm = TRUE),
       divergences = diag$divergences,
       fit = fit)
}

bayes_two_sample <- function(x1, x2, chains = 4L, iter = 4000L,
                             warmup = 2000L, seed = 1L) {
  stan_model <- HTDV:::.load_stan_model("simstudy_two_sample")
  fit <- suppressWarnings(
    rstan::sampling(stan_model,
                    data = list(N1 = length(x1), N2 = length(x2),
                                x1 = x1, x2 = x2),
                    chains = chains, iter = iter, warmup = warmup,
                    refresh = 0L, seed = seed,
                    control = list(adapt_delta = 0.95,
                                   max_treedepth = 12L))
  )
  draws <- as.numeric(rstan::extract(fit, pars = "delta",
                                     permuted = TRUE)$delta)
  sumr <- rstan::summary(fit,
                         pars = c("alpha1","alpha2","phi","sigma","delta"),
                         probs = c(0.025, 0.975))$summary
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divs <- sum(sapply(sp, function(p) sum(p[, "divergent__"])))
  list(point = mean(draws),
       ci = as.numeric(stats::quantile(draws, c(0.025, 0.975))),
       diag_pass = max(sumr[, "Rhat"], na.rm = TRUE) < 1.01 &&
                   min(sumr[, "n_eff"], na.rm = TRUE) > 200 && divs == 0L,
       rhat_max = max(sumr[, "Rhat"], na.rm = TRUE),
       ess_min = min(sumr[, "n_eff"], na.rm = TRUE),
       divergences = as.integer(divs),
       fit = fit)
}

boot_ci <- function(x, R = 1999L, alpha = 0.05, seed = 1L,
                    .stat = mean) {
  bb <- htdv_boot(x, .stat, R = R, type = "stationary",
                  block_length = "auto", level = 1 - alpha, seed = seed)
  list(point = bb$t0,
       ci = as.numeric(bb$ci_percentile),
       block_length = bb$block_length)
}

boot_ci_two_sample <- function(x1, x2, R = 1999L, alpha = 0.05,
                               seed = 1L) {
  set.seed(seed)
  b1 <- HTDV:::.ppw_block_length(x1)
  b2 <- HTDV:::.ppw_block_length(x2)
  diffs <- numeric(R)
  for (r in seq_len(R)) {
    xb1 <- HTDV:::.stationary_bootstrap_sample(x1, p = 1 / b1)
    xb2 <- HTDV:::.stationary_bootstrap_sample(x2, p = 1 / b2)
    diffs[r] <- mean(xb2) - mean(xb1)
  }
  q <- as.numeric(stats::quantile(diffs, probs = c(alpha / 2, 1 - alpha / 2),
                                  names = FALSE))
  list(point = mean(x2) - mean(x1),
       ci = q, block_length_x1 = b1, block_length_x2 = b2)
}

agreement_status <- function(ref, ci_lo, ci_hi) {
  if (!is.finite(ref) || !is.finite(ci_lo) || !is.finite(ci_hi)) return(NA_character_)
  if (ref >= ci_lo && ref <= ci_hi) "agreement"
  else if (ref < ci_lo) "htdv-strict (below ref)"
  else "htdv-strict (above ref)"
}

# =====================================================================
# E1: FRED-MD CPI inflation, post-1984
# =====================================================================
message("\n=== E1: FRED-MD CPI inflation (post-1984) ===")

fm_raw <- read.csv(fredmd_path, header = TRUE, na.strings = c("", "NA"),
                   stringsAsFactors = FALSE)
# Drop the second row labeled "Transform:"
fm <- fm_raw[!grepl("^Transform", fm_raw$sasdate, ignore.case = TRUE), ]
fm$date <- as.Date(fm$sasdate, format = "%m/%d/%Y")
cpi <- as.numeric(fm$CPIAUCSL)
keep <- !is.na(cpi) & !is.na(fm$date)
cpi <- cpi[keep]
fm_dates <- fm$date[keep]

# Annualized monthly CPI inflation
infl <- c(NA, 1200 * diff(log(cpi)))

post84 <- !is.na(infl) & fm_dates >= as.Date("1984-01-01")
x_e1 <- as.numeric(infl[post84])
n_e1 <- length(x_e1)
sample_start_e1 <- format(min(fm_dates[post84]), "%Y-%m")
sample_end_e1   <- format(max(fm_dates[post84]), "%Y-%m")
message(sprintf("  Sample: %s to %s | n=%d | mean=%.3f%% sd=%.3f%%",
                sample_start_e1, sample_end_e1, n_e1,
                mean(x_e1), stats::sd(x_e1)))

har_e1   <- har_ci_one_sample(x_e1)
bayes_e1 <- bayes_one_sample(x_e1, seed = 11L)
boot_e1  <- boot_ci(x_e1, seed = 11L)

ref_e1     <- 2.7
ref_e1_se  <- 0.4

E1 <- list(
  dataset = "FRED-MD CPI inflation post-1984",
  reference = "Stock & Watson (2007), JMCB",
  reference_value = ref_e1,
  reference_se = ref_e1_se,
  n = n_e1, sample = paste(sample_start_e1, "to", sample_end_e1),
  har = har_e1, bayes = bayes_e1, boot = boot_e1
)

# =====================================================================
# E2: Shiller log-CAPE mean reversion
# =====================================================================
message("\n=== E2: Shiller log-CAPE mean reversion ===")

# The Shiller xls has a multi-line preamble; the actual data table
# starts around row 8 in the "Data" sheet. Read with skip = 7 and look
# for the CAPE column by name pattern (it has been labeled "CAPE",
# "P/E10", "PE10", "P/E10 ratio" in different vintages).
sh <- as.data.frame(suppressWarnings(
  readxl::read_excel(shiller_path, sheet = "Data", skip = 7,
                     col_names = TRUE)
))
# Find date column
date_col <- grep("^Date", names(sh), ignore.case = TRUE, value = TRUE)[1]
if (is.na(date_col)) {
  date_col <- names(sh)[1]
  message("  No 'Date' column found; using the first column: ", date_col)
}
# Find CAPE column
cape_candidates <- grep("CAPE|P/?E.?10|PE.?10|Cyclically",
                        names(sh), ignore.case = TRUE, value = TRUE)
if (!length(cape_candidates))
  stop("Could not locate a CAPE-like column in the Shiller xls. ",
       "Available columns: ", paste(names(sh), collapse = ", "))
cape_col <- cape_candidates[1]
message("  Using column: ", cape_col)

cape_raw <- suppressWarnings(as.numeric(sh[[cape_col]]))
keep <- is.finite(cape_raw) & cape_raw > 0
cape <- cape_raw[keep]
log_cape <- log(cape)
n_e2 <- length(log_cape)
message(sprintf("  n=%d | mean log CAPE = %.4f | mean CAPE = %.2f",
                n_e2, mean(log_cape), exp(mean(log_cape))))

har_e2   <- har_ci_one_sample(log_cape)
bayes_e2 <- bayes_one_sample(log_cape, seed = 12L)
boot_e2  <- boot_ci(log_cape, seed = 12L)

ref_e2 <- 2.85   # Campbell-Shiller 1998: long-run historical mean log-CAPE

E2 <- list(
  dataset = "Shiller log-CAPE",
  reference = "Campbell & Shiller (1998)",
  reference_value = ref_e2,
  reference_se = NA_real_,
  n = n_e2, sample = "1881-01 to latest",
  har = har_e2, bayes = bayes_e2, boot = boot_e2
)

# =====================================================================
# E3: US vs Canada 10-year yield, mean equality
# =====================================================================
message("\n=== E3: US vs Canada 10-year yield, mean equality ===")

uc <- read.csv(us_canada_path, header = TRUE, na.strings = c("", ".", "NA"),
               stringsAsFactors = FALSE)
uc$date <- as.Date(uc$observation_date)
uc$gs10 <- suppressWarnings(as.numeric(uc$GS10))
uc$ca10 <- suppressWarnings(as.numeric(uc$IRLTLT01CAM156N))
common <- !is.na(uc$gs10) & !is.na(uc$ca10) & uc$date >= as.Date("1990-01-01")
x_us <- uc$gs10[common]
x_ca <- uc$ca10[common]
n_us <- length(x_us); n_ca <- length(x_ca)
sample_start_e3 <- format(min(uc$date[common]), "%Y-%m")
sample_end_e3   <- format(max(uc$date[common]), "%Y-%m")
message(sprintf("  Sample: %s to %s | n_US=%d n_CA=%d | mean_US=%.3f mean_CA=%.3f",
                sample_start_e3, sample_end_e3, n_us, n_ca,
                mean(x_us), mean(x_ca)))

har_e3   <- har_ci_two_sample(x_us, x_ca)
bayes_e3 <- bayes_two_sample(x_us, x_ca, seed = 13L)
boot_e3  <- boot_ci_two_sample(x_us, x_ca, seed = 13L)

# No theoretical reference; we report observed differential and let
# the empirical result speak for itself. We pin a "naive iid t-test"
# comparison as the would-be classical answer.
naive_se_e3 <- sqrt(stats::var(x_us) / n_us + stats::var(x_ca) / n_ca)
naive_diff_e3 <- mean(x_ca) - mean(x_us)
naive_ci_e3 <- naive_diff_e3 + c(-1, 1) * stats::qnorm(0.975) * naive_se_e3

E3 <- list(
  dataset = "US vs Canada 10y yields (mean differential CA - US)",
  reference = "Naive Welch t-test (iid baseline)",
  reference_value = naive_diff_e3,
  reference_se = naive_se_e3,
  reference_naive_ci = naive_ci_e3,
  n = n_us, sample = paste(sample_start_e3, "to", sample_end_e3),
  har = har_e3, bayes = bayes_e3, boot = boot_e3
)

# =====================================================================
# Aggregate and write
# =====================================================================
message("\n=== Aggregating ===")

row_for <- function(E) {
  data.frame(
    dataset = E$dataset,
    reference = E$reference,
    n = E$n,
    sample = E$sample,
    reference_value = E$reference_value,
    reference_se = E$reference_se,
    har_point = E$har$point,
    har_ci_lo = E$har$ci[1L],
    har_ci_hi = E$har$ci[2L],
    bayes_point = E$bayes$point,
    bayes_ci_lo = E$bayes$ci[1L],
    bayes_ci_hi = E$bayes$ci[2L],
    bayes_diag_pass = E$bayes$diag_pass,
    bayes_rhat_max = E$bayes$rhat_max,
    bayes_ess_min = E$bayes$ess_min,
    bayes_divergences = E$bayes$divergences,
    boot_point = E$boot$point,
    boot_ci_lo = E$boot$ci[1L],
    boot_ci_hi = E$boot$ci[2L],
    agreement_har = agreement_status(E$reference_value,
                                     E$har$ci[1L], E$har$ci[2L]),
    agreement_bayes = agreement_status(E$reference_value,
                                       E$bayes$ci[1L], E$bayes$ci[2L]),
    agreement_boot = agreement_status(E$reference_value,
                                      E$boot$ci[1L], E$boot$ci[2L]),
    stringsAsFactors = FALSE
  )
}

results_df <- rbind(row_for(E1), row_for(E2), row_for(E3))
utils::write.csv(results_df, results_csv, row.names = FALSE)
saveRDS(list(E1 = E1, E2 = E2, E3 = E3, summary = results_df),
        results_rds)

message("\nWrote: ", results_csv)
message("Wrote: ", results_rds)

message("\n=== Summary table ===")
for (i in seq_len(nrow(results_df))) {
  r <- results_df[i, ]
  message(sprintf("\n  %s", r$dataset))
  message(sprintf("    n=%d | sample: %s", r$n, r$sample))
  message(sprintf("    Reference (%s): %.4f%s",
                  r$reference, r$reference_value,
                  if (is.finite(r$reference_se))
                    sprintf(" (SE %.4f)", r$reference_se) else ""))
  message(sprintf("    HAR:    point=%.4f  CI=[%.4f, %.4f]  vs ref: %s",
                  r$har_point, r$har_ci_lo, r$har_ci_hi, r$agreement_har))
  message(sprintf("    Bayes:  point=%.4f  CI=[%.4f, %.4f]  vs ref: %s",
                  r$bayes_point, r$bayes_ci_lo, r$bayes_ci_hi,
                  r$agreement_bayes))
  message(sprintf("    Boot:   point=%.4f  CI=[%.4f, %.4f]  vs ref: %s",
                  r$boot_point, r$boot_ci_lo, r$boot_ci_hi,
                  r$agreement_boot))
  message(sprintf("    Bayes diagnostics: pass=%s  Rhat_max=%.4f  min-ESS=%.0f  divergences=%d",
                  r$bayes_diag_pass, r$bayes_rhat_max,
                  r$bayes_ess_min, r$bayes_divergences))
}
