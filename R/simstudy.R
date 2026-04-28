#' Factorial Monte Carlo Simulation Study for Dependent-Unbalanced Data
#'
#' Runs the pre-registered factorial Monte Carlo study of Section 12-bis of
#' the companion paper. The study crosses four factors (per-group sample
#' size, AR(1) coefficient, innovation tail, imbalance ratio) and evaluates
#' three inferential layers (HAR-Wald, stationary block bootstrap,
#' hierarchical Bayesian HMC) on each replication.
#'
#' @section Data-generating process:
#' Two independent AR(1) series \eqn{x^{(1)}_t} of length \eqn{n_1=n} and
#' \eqn{x^{(2)}_t} of length \eqn{n_2 = \max(2, \mathrm{round}(n / \mathrm{imb}))}
#' are generated. The innovations are drawn from \code{rnorm},
#' scaled Student-\eqn{t_{(5)}}, scaled Student-\eqn{t_{(3)}} or scaled
#' Student-\eqn{t_{(2.1)}} so that the nominal innovation variance is one
#' when the law has a finite second moment. Group one has mean zero; group
#' two has mean \eqn{\Delta\cdot\sigma_\infty} where
#' \eqn{\sigma_\infty = 1/\sqrt{1 - \phi^2}}.
#'
#' @section Inferential layers:
#' \itemize{
#'   \item \strong{HAR-Wald (\code{"har"}).} Welch-style two-sample statistic
#'     \eqn{z = (\bar x_1 - \bar x_2) / \sqrt{\hat\sigma_1^2/n_1 + \hat\sigma_2^2/n_2}}
#'     with long-run variances estimated via \code{\link{htdv_lrv}} (QS
#'     kernel, Andrews bandwidth) on each group's demeaned series.
#'     Rejection at asymptotic \eqn{\chi^2_1} critical value.
#'   \item \strong{Stationary block bootstrap (\code{"boot"}).} Each group
#'     resampled independently using
#'     \code{\link{htdv_boot}} with \code{block_length = "auto"} (Patton-
#'     Politis-White). Percentile interval on
#'     \eqn{\bar x_1^* - \bar x_2^*}.
#'   \item \strong{Hierarchical Bayesian (\code{"bayes"}).} Two-sample AR(1)
#'     Stan model (\code{simstudy_two_sample.stan}); posterior interval on
#'     \eqn{\delta = \alpha_2 - \alpha_1}. A run is flagged
#'     \code{diag_pass = TRUE} only if \eqn{\hat R < 1.01}, bulk and tail
#'     ESS above \eqn{200}, and zero divergences.
#' }
#'
#' @param n_grid Integer vector of primary sample sizes \eqn{n_1}.
#' @param phi_grid Numeric vector of AR(1) coefficients in \eqn{(-1,1)}.
#' @param tail_grid Character vector; any subset of \code{"normal"},
#'   \code{"t5"}, \code{"t3"}, \code{"t2_1"}.
#' @param imb_grid Numeric vector of imbalance ratios \eqn{n_1/n_2}.
#' @param delta_grid Numeric vector of location shifts in units of
#'   \eqn{\sigma_\infty}.
#' @param R Integer number of Monte Carlo replications per cell.
#' @param seed Integer master seed.
#' @param n_cores Integer number of workers for cell-level parallelism.
#'   On non-Unix platforms parallelism falls back to serial.
#' @param layers Character vector; any subset of \code{c("har","boot","bayes")}.
#' @param bayes_chains Number of HMC chains per Bayesian fit.
#' @param bayes_iter Total HMC iterations (warmup + sampling).
#' @param bayes_warmup HMC warmup iterations.
#' @param boot_R Bootstrap replicates per call.
#' @param rope Length-2 numeric vector; ROPE for the Bayesian decision on
#'   the raw (non-standardized) delta.
#' @param alpha Significance level (nominal size of tests; \eqn{1-\alpha} is
#'   the coverage target for intervals).
#' @param out_dir Optional directory for per-cell incremental RDS results.
#'   Each cell is written to its own \code{cell_<id>.rds} so that concurrent
#'   workers cannot collide on a shared file. Final aggregation reads the
#'   directory back, so successful cells are preserved even if some workers
#'   later die. If \code{NULL}, no per-cell artifacts are written and the
#'   only result is the in-memory return value.
#' @param progress Logical; print a one-line status per completed cell.
#'
#' @return A data frame with one row per (cell, replication, layer). Columns:
#' \code{cell_id, n, phi, tail, imb, n1, n2, delta, sigma_inf, layer,
#' replicate, reject, ci_lo, ci_hi, estimate, covered, ci_length,
#' rhat_max, ess_min, divergences, diag_pass, runtime_sec}.
#'
#' @references
#' Kiefer, N.M., & Vogelsang, T.J. (2005). Econometric Theory 21(6):
#' 1130-1164.
#' Patton, A., Politis, D.N., & White, H. (2009). Econometric Reviews
#' 28(4): 372-375.
#' Kruschke, J.K. (2018). Advances in Methods and Practices in
#' Psychological Science 1(2): 270-280.
#'
#' @examples
#' \donttest{
#' res <- htdv_simstudy(n_grid = c(40, 80),
#'                      phi_grid = c(0, 0.6),
#'                      tail_grid = c("normal"),
#'                      imb_grid = c(1, 3),
#'                      delta_grid = c(0, 0.25),
#'                      R = 5L, n_cores = 1L,
#'                      layers = c("har", "boot"))
#' head(htdv_simstudy_summary(res))
#' }
#'
#' @export
htdv_simstudy <- function(n_grid = c(40L, 80L, 200L, 500L),
                          phi_grid = c(0, 0.3, 0.6, 0.85),
                          tail_grid = c("normal", "t5", "t3", "t2_1"),
                          imb_grid = c(1, 1.5, 3, 6),
                          delta_grid = c(0, 0.1, 0.25, 0.5),
                          R = 500L,
                          seed = 20260422L,
                          n_cores = 1L,
                          layers = c("har", "boot", "bayes"),
                          bayes_chains = 2L,
                          bayes_iter = 600L,
                          bayes_warmup = 300L,
                          boot_R = 999L,
                          rope = c(-0.1, 0.1),
                          alpha = 0.05,
                          out_dir = NULL,
                          progress = TRUE) {

  layers <- match.arg(layers, several.ok = TRUE,
                      choices = c("har", "boot", "bayes"))
  n_grid <- as.integer(n_grid)
  R <- as.integer(R)
  if (any(n_grid < 10L)) stop("n_grid entries must be >= 10.", call. = FALSE)
  if (R < 2L) stop("R must be >= 2.", call. = FALSE)
  if (!is.numeric(rope) || length(rope) != 2L || rope[1L] >= rope[2L])
    stop("'rope' must be a length-2 increasing numeric.", call. = FALSE)

  stan_model <- if ("bayes" %in% layers) {
    .load_stan_model("simstudy_two_sample")
  } else NULL

  grid <- expand.grid(n = n_grid, phi = phi_grid, tail = tail_grid,
                      imb = imb_grid, delta = delta_grid,
                      stringsAsFactors = FALSE)
  grid$cell_id <- seq_len(nrow(grid))

  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (progress) {
    message(sprintf("htdv_simstudy: %d cells x %d replications x %d layers.",
                    nrow(grid), R, length(layers)))
  }

  cell_path <- function(i) {
    file.path(out_dir, sprintf("cell_%05d.rds", i))
  }

  run_one <- function(i) {
    if (!is.null(out_dir) && file.exists(cell_path(i))) {
      if (progress) {
        message(sprintf("  cell %d/%d | resumed from cache, skipped",
                        i, nrow(grid)))
      }
      return(invisible(NULL))
    }
    spec <- grid[i, , drop = FALSE]
    cell_seed <- .simstudy_cell_seed(seed, i)
    out <- .simstudy_run_cell(spec = spec, R = R, cell_seed = cell_seed,
                              layers = layers, stan_model = stan_model,
                              bayes_chains = bayes_chains,
                              bayes_iter = bayes_iter,
                              bayes_warmup = bayes_warmup,
                              boot_R = boot_R, rope = rope, alpha = alpha)
    if (!is.null(out_dir)) {
      .simstudy_append_cell(out, out_dir)
    }
    if (progress) {
      tot <- sum(out$runtime_sec, na.rm = TRUE)
      message(sprintf("  cell %d/%d | n=%d phi=%.2f tail=%s imb=%.2f delta=%.2f | %.1fs",
                      i, nrow(grid), spec$n, spec$phi, spec$tail,
                      spec$imb, spec$delta, tot))
    }
    if (is.null(out_dir)) out else invisible(NULL)
  }

  use_mc <- identical(.Platform$OS.type, "unix") && n_cores > 1L
  worker_results <- if (use_mc) {
    parallel::mclapply(seq_len(nrow(grid)), run_one,
                       mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    lapply(seq_len(nrow(grid)), run_one)
  }

  n_worker_errors <- sum(vapply(worker_results,
                                function(x) inherits(x, "try-error"),
                                logical(1L)))
  if (n_worker_errors > 0L) {
    warning(sprintf("htdv_simstudy: %d/%d workers errored. Surviving cells are preserved on disk.",
                    n_worker_errors, length(worker_results)),
            call. = FALSE)
  }

  if (!is.null(out_dir)) {
    files <- list.files(out_dir, pattern = "^cell_\\d+\\.rds$",
                        full.names = TRUE)
    if (!length(files)) {
      stop("No per-cell result files found in '", out_dir, "'.", call. = FALSE)
    }
    parts <- lapply(files, function(f) {
      tryCatch(readRDS(f), error = function(e) NULL)
    })
    ok <- vapply(parts, is.data.frame, logical(1L))
    if (sum(ok) < length(parts)) {
      warning(sprintf("htdv_simstudy: %d cell file(s) unreadable; skipped.",
                      sum(!ok)), call. = FALSE)
    }
    do.call(rbind, parts[ok])
  } else {
    df_ok <- vapply(worker_results, is.data.frame, logical(1L))
    if (!any(df_ok)) {
      stop("No worker returned usable data and out_dir was NULL.",
           call. = FALSE)
    }
    do.call(rbind, worker_results[df_ok])
  }
}

#' Flag Simulation-Study Cells in the Limit-of-Identification Zone
#'
#' Reads an aggregated \code{\link{htdv_simstudy_summary}} output and returns
#' the subset of cells where the Bayesian diagnostic-pass rate falls below a
#' user-set threshold. Such cells are typically located in the corner of the
#' design where strong autocorrelation meets a small sample (high \eqn{\phi}
#' and low \eqn{n}); the AR(1) likelihood there approaches its limit of
#' identification, MCMC diagnostics flag the resulting posterior, and
#' practitioners should treat the cell with extra caution—either by
#' increasing HMC iterations, tightening the priors on \eqn{\phi}, or
#' falling back to the conformal anchor whose validity does not require
#' identifiability.
#'
#' @param summ Data frame from \code{\link{htdv_simstudy_summary}}.
#' @param threshold Minimum acceptable Bayesian \code{diag_pass_rate}.
#'   Default \code{0.7}.
#'
#' @return A data frame with the columns \code{n, phi, tail, imb, delta,
#'   layer, diag_pass_rate}, restricted to Bayesian rows whose pass rate is
#'   below \code{threshold}, plus an integer column \code{flag_index}
#'   (sequential identifier of flagged cells, useful for cross-referencing
#'   in reports and in vignettes).
#'
#' @examples
#' \donttest{
#' res <- htdv_simstudy(n_grid = c(40, 200),
#'                      phi_grid = c(0, 0.85),
#'                      tail_grid = "normal",
#'                      imb_grid = c(1, 6),
#'                      delta_grid = c(0, 0.25),
#'                      R = 50L)
#' summ <- htdv_simstudy_summary(res)
#' htdv_simstudy_warnings(summ, threshold = 0.7)
#' }
#'
#' @export
htdv_simstudy_warnings <- function(summ, threshold = 0.7) {
  if (!is.data.frame(summ)) stop("'summ' must be a data frame.", call. = FALSE)
  required <- c("n", "phi", "tail", "imb", "delta", "layer", "diag_pass_rate")
  missing_cols <- setdiff(required, names(summ))
  if (length(missing_cols)) {
    stop("'summ' is missing columns: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      threshold <= 0 || threshold >= 1) {
    stop("'threshold' must be a scalar in (0, 1).", call. = FALSE)
  }
  bayes <- summ[summ$layer == "bayes", , drop = FALSE]
  flagged <- bayes[is.finite(bayes$diag_pass_rate) &
                   bayes$diag_pass_rate < threshold, required, drop = FALSE]
  flagged <- flagged[order(flagged$diag_pass_rate), , drop = FALSE]
  rownames(flagged) <- NULL
  flagged$flag_index <- seq_len(nrow(flagged))
  flagged
}

#' Aggregate a Simulation-Study Output
#'
#' Computes empirical size (or rejection rate), coverage, and mean interval
#' length per (cell x layer) from the raw output of \code{\link{htdv_simstudy}}.
#'
#' @param res Data frame returned by \code{htdv_simstudy}.
#' @param alpha Nominal size used in the study (must match what
#'   \code{htdv_simstudy} was called with).
#'
#' @return A data frame with one row per (cell, layer) and columns
#' \code{n, phi, tail, imb, delta, layer, n_reps, reject_rate, coverage,
#' mean_ci_length, mean_runtime, diag_pass_rate}.
#'
#' @examples
#' \donttest{
#' res <- htdv_simstudy(n_grid = c(40), phi_grid = c(0.3),
#'                      tail_grid = "normal", imb_grid = 1,
#'                      delta_grid = c(0, 0.3), R = 5,
#'                      layers = c("har"))
#' htdv_simstudy_summary(res)
#' }
#'
#' @export
htdv_simstudy_summary <- function(res, alpha = 0.05) {
  if (!is.data.frame(res)) stop("'res' must be a data frame.", call. = FALSE)
  by_cols <- c("n", "phi", "tail", "imb", "delta", "layer")
  key <- do.call(paste, c(res[by_cols], sep = "|"))
  split_res <- split(res, key)
  out <- lapply(split_res, function(d) {
    data.frame(
      n = d$n[1L], phi = d$phi[1L], tail = d$tail[1L],
      imb = d$imb[1L], delta = d$delta[1L], layer = d$layer[1L],
      n_reps = nrow(d),
      reject_rate = mean(d$reject, na.rm = TRUE),
      coverage = mean(d$covered, na.rm = TRUE),
      mean_ci_length = mean(d$ci_length, na.rm = TRUE),
      mean_runtime = mean(d$runtime_sec, na.rm = TRUE),
      diag_pass_rate = if (any(!is.na(d$diag_pass)))
        mean(d$diag_pass, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  agg <- do.call(rbind, out)
  rownames(agg) <- NULL
  agg
}

# --- Internal helpers --------------------------------------------------------

.simstudy_cell_seed <- function(master, i) {
  as.integer((as.numeric(master) + 1e4 * i) %% .Machine$integer.max)
}

.simstudy_rep_seed <- function(cell_seed, r) {
  as.integer((as.numeric(cell_seed) + r) %% .Machine$integer.max)
}

.simstudy_append_cell <- function(rows, dir) {
  cell_id <- rows$cell_id[1L]
  tmp <- tempfile(tmpdir = dir, pattern = sprintf("cell_%05d_", cell_id),
                  fileext = ".rds.tmp")
  saveRDS(rows, file = tmp)
  final <- file.path(dir, sprintf("cell_%05d.rds", cell_id))
  file.rename(tmp, final)
  invisible(TRUE)
}

.simstudy_gen_ar1 <- function(n, phi, tail, mean_shift = 0, burn = 200L) {
  m <- n + burn
  eps <- switch(tail,
    "normal" = stats::rnorm(m),
    "t5"     = stats::rt(m, df = 5) / sqrt(5 / 3),
    "t3"     = stats::rt(m, df = 3) / sqrt(3 / 1),
    "t2_1"   = stats::rt(m, df = 2.1) / sqrt(2.1 / 0.1),
    stop("unknown tail: ", tail, call. = FALSE)
  )
  x <- numeric(m)
  x[1L] <- eps[1L]
  for (t in 2:m) x[t] <- phi * x[t - 1L] + eps[t]
  x[(burn + 1L):m] + mean_shift
}

.simstudy_run_cell <- function(spec, R, cell_seed, layers, stan_model,
                               bayes_chains, bayes_iter, bayes_warmup,
                               boot_R, rope, alpha) {
  n1 <- spec$n
  n2 <- max(2L, as.integer(round(spec$n / spec$imb)))
  sigma_inf <- 1 / sqrt(max(1 - spec$phi^2, 1e-12))
  true_delta <- spec$delta * sigma_inf
  rows <- vector("list", R * length(layers))
  k <- 0L
  for (r in seq_len(R)) {
    rep_seed <- .simstudy_rep_seed(cell_seed, r)
    set.seed(rep_seed)
    x1 <- .simstudy_gen_ar1(n1, spec$phi, spec$tail, 0)
    x2 <- .simstudy_gen_ar1(n2, spec$phi, spec$tail, true_delta)

    for (lyr in layers) {
      t0 <- proc.time()[["elapsed"]]
      layer_out <- tryCatch(
        switch(lyr,
          "har"   = .simstudy_layer_har(x1, x2, alpha),
          "boot"  = .simstudy_layer_boot(x1, x2, boot_R, alpha, rep_seed),
          "bayes" = .simstudy_layer_bayes(x1, x2, stan_model,
                                          bayes_chains, bayes_iter,
                                          bayes_warmup, rope, alpha,
                                          rep_seed)
        ),
        error = function(e) list(reject = NA, ci_lo = NA_real_,
                                 ci_hi = NA_real_, estimate = NA_real_,
                                 rhat_max = NA_real_, ess_min = NA_real_,
                                 divergences = NA_integer_,
                                 diag_pass = NA, err = conditionMessage(e))
      )
      elapsed <- proc.time()[["elapsed"]] - t0
      covered <- if (is.finite(layer_out$ci_lo) && is.finite(layer_out$ci_hi))
        (true_delta >= layer_out$ci_lo) && (true_delta <= layer_out$ci_hi)
      else NA
      ci_length <- if (is.finite(layer_out$ci_lo) && is.finite(layer_out$ci_hi))
        layer_out$ci_hi - layer_out$ci_lo else NA_real_

      k <- k + 1L
      rows[[k]] <- data.frame(
        cell_id   = spec$cell_id,
        n         = n1,
        phi       = spec$phi,
        tail      = spec$tail,
        imb       = spec$imb,
        n1        = n1, n2 = n2,
        delta     = spec$delta,
        sigma_inf = sigma_inf,
        layer     = lyr,
        replicate = r,
        reject    = as.logical(layer_out$reject),
        ci_lo     = layer_out$ci_lo,
        ci_hi     = layer_out$ci_hi,
        estimate  = layer_out$estimate,
        covered   = covered,
        ci_length = ci_length,
        rhat_max    = if (!is.null(layer_out$rhat_max)) layer_out$rhat_max
                      else NA_real_,
        ess_min     = if (!is.null(layer_out$ess_min)) layer_out$ess_min
                      else NA_real_,
        divergences = if (!is.null(layer_out$divergences))
                      layer_out$divergences else NA_integer_,
        diag_pass   = if (!is.null(layer_out$diag_pass)) layer_out$diag_pass
                      else NA,
        runtime_sec = elapsed,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

.simstudy_layer_har <- function(x1, x2, alpha) {
  # Sign convention: estimand is delta = E[X2] - E[X1], matching the Stan
  # layer's parameter delta = alpha2 - alpha1. Coverage of CI is checked
  # against the same true_delta = +shift used by the data generator.
  n1 <- length(x1); n2 <- length(x2)
  m1 <- mean(x1); m2 <- mean(x2)
  lrv1 <- htdv_lrv(x1, kernel = "qs", bandwidth = "andrews")$lrv
  lrv2 <- htdv_lrv(x2, kernel = "qs", bandwidth = "andrews")$lrv
  se <- sqrt(lrv1 / n1 + lrv2 / n2)
  est <- m2 - m1
  z <- est / se
  z_crit <- stats::qnorm(1 - alpha / 2)
  list(reject = abs(z) > z_crit,
       ci_lo = est - z_crit * se,
       ci_hi = est + z_crit * se,
       estimate = est,
       rhat_max = NA_real_, ess_min = NA_real_,
       divergences = NA_integer_, diag_pass = NA)
}

.simstudy_layer_boot <- function(x1, x2, R, alpha, seed) {
  # Sign convention: see .simstudy_layer_har. delta = E[X2] - E[X1].
  n1 <- length(x1); n2 <- length(x2)
  b1 <- .ppw_block_length(x1)
  b2 <- .ppw_block_length(x2)
  set.seed(seed)
  diffs <- numeric(R)
  for (r in seq_len(R)) {
    xb1 <- .stationary_bootstrap_sample(x1, p = 1 / b1)
    xb2 <- .stationary_bootstrap_sample(x2, p = 1 / b2)
    diffs[r] <- mean(xb2) - mean(xb1)
  }
  est <- mean(x2) - mean(x1)
  q <- as.numeric(stats::quantile(diffs, probs = c(alpha / 2, 1 - alpha / 2),
                                  names = FALSE))
  list(reject = (q[1L] > 0) || (q[2L] < 0),
       ci_lo = q[1L], ci_hi = q[2L], estimate = est,
       rhat_max = NA_real_, ess_min = NA_real_,
       divergences = NA_integer_, diag_pass = NA)
}

.simstudy_layer_bayes <- function(x1, x2, stan_model, chains, iter, warmup,
                                  rope, alpha, seed) {
  sd <- list(N1 = length(x1), N2 = length(x2), x1 = x1, x2 = x2)
  fit <- suppressWarnings(
    rstan::sampling(stan_model, data = sd, chains = chains,
                    iter = iter, warmup = warmup, refresh = 0L,
                    seed = seed,
                    control = list(adapt_delta = 0.95,
                                   max_treedepth = 12L))
  )
  draws <- as.numeric(rstan::extract(fit, pars = "delta",
                                     permuted = TRUE)$delta)
  ci <- as.numeric(stats::quantile(draws,
                                   probs = c(alpha / 2, 1 - alpha / 2),
                                   names = FALSE))
  est <- mean(draws)
  rej_ci <- (ci[1L] > 0) || (ci[2L] < 0)
  sf_sum <- rstan::summary(fit,
                           pars = c("alpha1","alpha2","phi","sigma","delta"),
                           probs = c(0.025, 0.975))$summary
  rhat_max <- max(sf_sum[, "Rhat"], na.rm = TRUE)
  ess_min <- min(sf_sum[, "n_eff"], na.rm = TRUE)
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divs <- sum(sapply(sp, function(p) sum(p[, "divergent__"])))
  diag_pass <- is.finite(rhat_max) && rhat_max < 1.01 &&
               is.finite(ess_min) && ess_min > 200 && divs == 0L
  list(reject = rej_ci,
       ci_lo = ci[1L], ci_hi = ci[2L], estimate = est,
       rhat_max = rhat_max, ess_min = ess_min,
       divergences = as.integer(divs), diag_pass = diag_pass)
}
