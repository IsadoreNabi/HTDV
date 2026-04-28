#' MCMC Diagnostics for htdv_fit Objects
#'
#' Computes split-Rhat, bulk and tail ESS, E-BFMI, divergence count, and
#' treedepth saturation; all standard HMC diagnostics (Vehtari et al. 2021,
#' Betancourt 2016).
#'
#' @param fit An \code{htdv_fit} object.
#' @param rhat_threshold Convergence threshold for Rhat.
#' @param ess_threshold Minimum effective sample size.
#' @param bfmi_threshold Minimum acceptable energy-BFMI.
#'
#' @return A list with \code{rhat}, \code{ess_bulk}, \code{ess_tail},
#'   \code{bfmi}, \code{divergences}, \code{max_treedepth},
#'   \code{passed} (logical), and \code{failures} (character vector).
#'
#' @references
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Buerkner, P.-C.
#' (2021). Rank-normalization, folding, and localization: An improved Rhat
#' for assessing convergence of MCMC. Bayesian Analysis 16(2): 667-718.
#' Betancourt, M. (2016). Diagnosing Suboptimal Cotangent Disintegrations.
#' <doi:10.48550/arXiv.1604.00695>.
#'
#' @examples
#' \donttest{
#' x <- rnorm(50)
#' fit <- htdv_fit(x, model = "tac", chains = 2, iter = 500,
#'                 refresh = 0, seed = 1)
#' htdv_diagnostics(fit)$passed
#' }
#'
#' @export
htdv_diagnostics <- function(fit, rhat_threshold = 1.01,
                             ess_threshold = 400, bfmi_threshold = 0.3) {
  if (!inherits(fit, "htdv_fit")) {
    stop("'fit' must be an 'htdv_fit' object.", call. = FALSE)
  }
  sf <- fit$stanfit
  mon <- rstan::monitor(sf, print = FALSE, warmup = 0)
  rh <- mon[, "Rhat"]
  essb <- mon[, "Bulk_ESS"]
  esst <- mon[, "Tail_ESS"]
  sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
  div <- sum(sapply(sp, function(p) sum(p[, "divergent__"])))
  tdmax <- max(sapply(sp, function(p) max(p[, "treedepth__"])))
  energy_bfmi <- vapply(sp, function(p) {
    e <- p[, "energy__"]
    var_e <- stats::var(e)
    if (!is.finite(var_e) || var_e <= 0) return(NA_real_)
    sum(diff(e)^2) / (length(e) - 1L) / var_e
  }, numeric(1L))
  fails <- character(0L)
  if (any(rh > rhat_threshold, na.rm = TRUE)) {
    fails <- c(fails, sprintf("Rhat > %.3f on some parameter", rhat_threshold))
  }
  if (any(essb < ess_threshold, na.rm = TRUE)) {
    fails <- c(fails, "Bulk ESS below threshold")
  }
  if (any(esst < ess_threshold, na.rm = TRUE)) {
    fails <- c(fails, "Tail ESS below threshold")
  }
  if (div > 0L) {
    fails <- c(fails, sprintf("%d post-warmup divergences", div))
  }
  if (any(is.na(energy_bfmi)) || any(energy_bfmi < bfmi_threshold, na.rm = TRUE)) {
    fails <- c(fails, "E-BFMI below threshold")
  }
  list(rhat = rh, ess_bulk = essb, ess_tail = esst,
       bfmi = energy_bfmi, divergences = div,
       max_treedepth = tdmax,
       passed = length(fails) == 0L,
       failures = fails)
}

#' Posterior-Predictive Checks on Dependence Statistics
#'
#' Replicates autocorrelation, Ljung-Box, spectral discrepancy, and a
#' detrended-fluctuation Hurst exponent under the posterior and computes
#' posterior-predictive \code{p}-values.
#'
#' @param fit An \code{htdv_fit} object whose posterior-predictive replicates
#'   are stored as parameter \code{x_rep}.
#' @param x_obs The observed data vector.
#' @param lag_max Maximum lag for autocorrelation.
#'
#' @return A list with \code{p_acf}, \code{p_ljung_box}, \code{p_spectral},
#'   \code{p_hurst}.
#'
#' @examples
#' \donttest{
#' x <- rnorm(50)
#' fit <- htdv_fit(x, model = "tac", chains = 2, iter = 500,
#'                 refresh = 0, seed = 1)
#' htdv_ppc(fit, x)
#' }
#'
#' @export
htdv_ppc <- function(fit, x_obs, lag_max = 10L) {
  if (!inherits(fit, "htdv_fit")) {
    stop("'fit' must be an 'htdv_fit' object.", call. = FALSE)
  }
  .assert_numeric_vector(x_obs, "x_obs")
  sf <- fit$stanfit
  pars <- sf@model_pars
  if (!"x_rep" %in% pars) {
    stop("Posterior-predictive replicate 'x_rep' not found.", call. = FALSE)
  }
  xrep <- rstan::extract(sf, pars = "x_rep", permuted = TRUE)$x_rep
  t_obs_acf <- sum(stats::acf(x_obs, lag.max = lag_max, plot = FALSE)$acf[-1L]^2)
  t_obs_lb <- .ljung_box(x_obs, lag = max(1L, floor(log(length(x_obs)))))
  t_obs_sp <- .spectral_discrepancy(x_obs)
  t_obs_h <- .dfa_hurst(x_obs)
  t_rep_acf <- apply(xrep, 1L, function(xr) {
    sum(stats::acf(xr, lag.max = lag_max, plot = FALSE)$acf[-1L]^2)
  })
  t_rep_lb <- apply(xrep, 1L, function(xr) {
    .ljung_box(xr, lag = max(1L, floor(log(length(xr)))))
  })
  t_rep_sp <- apply(xrep, 1L, .spectral_discrepancy)
  t_rep_h <- apply(xrep, 1L, .dfa_hurst)
  list(p_acf = mean(t_rep_acf >= t_obs_acf),
       p_ljung_box = mean(t_rep_lb >= t_obs_lb),
       p_spectral = mean(t_rep_sp >= t_obs_sp),
       p_hurst = mean(t_rep_h >= t_obs_h))
}

.ljung_box <- function(x, lag = 10L) {
  n <- length(x)
  r <- stats::acf(x, lag.max = lag, plot = FALSE)$acf[-1L]
  n * (n + 2) * sum(r^2 / (n - seq_along(r)))
}

.spectral_discrepancy <- function(x) {
  sp <- stats::spec.pgram(x, plot = FALSE, taper = 0, detrend = FALSE,
                          demean = TRUE)
  mean(abs(sp$spec - mean(sp$spec))^2)
}

.dfa_hurst <- function(x) {
  n <- length(x)
  if (n < 16L) return(NA_real_)
  y <- cumsum(x - mean(x))
  scales <- unique(round(exp(seq(log(4), log(n / 4), length.out = 10))))
  scales <- scales[scales >= 4L & scales <= n / 4L]
  if (length(scales) < 2L) return(NA_real_)
  fs <- numeric(length(scales))
  for (i in seq_along(scales)) {
    s <- scales[i]
    k <- floor(n / s)
    if (k < 1L) { fs[i] <- NA_real_; next }
    rs <- numeric(k)
    for (j in seq_len(k)) {
      seg <- y[((j - 1L) * s + 1L):(j * s)]
      tseg <- seq_len(s)
      beta <- stats::.lm.fit(cbind(1, tseg), seg)$coefficients
      fit <- beta[1L] + beta[2L] * tseg
      rs[j] <- mean((seg - fit)^2)
    }
    fs[i] <- sqrt(mean(rs, na.rm = TRUE))
  }
  ok <- is.finite(fs) & fs > 0
  if (sum(ok) < 2L) return(NA_real_)
  coef <- stats::.lm.fit(cbind(1, log(scales[ok])), log(fs[ok]))$coefficients
  as.numeric(coef[2L])
}
