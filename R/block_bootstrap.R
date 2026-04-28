#' Block Bootstrap with Automatic Block Length
#'
#' Implements circular and stationary block bootstraps (Politis-Romano 1992,
#' 1994) with Patton-Politis-White (2009) automatic block length selection.
#'
#' @param x Numeric vector.
#' @param statistic Function \code{function(x)} returning a scalar or numeric
#'   vector.
#' @param R Number of bootstrap replicates.
#' @param type \code{"circular"} or \code{"stationary"}.
#' @param block_length Either \code{"auto"} (Patton-Politis-White) or a
#'   positive integer.
#' @param level Confidence level for intervals.
#' @param seed Optional integer seed.
#'
#' @return A list with \code{t0} (statistic on original data), \code{t_star}
#'   (bootstrap distribution), \code{ci_percentile}, \code{ci_basic},
#'   \code{ci_studentized}, and \code{block_length}.
#'
#' @references
#' Politis, D.N., & Romano, J.P. (1994). The Stationary Bootstrap. JASA
#' 89(428): 1303-1313.
#' Patton, A., Politis, D.N., & White, H. (2009). Correction to Automatic
#' Block-Length Selection for the Dependent Bootstrap. Econometric Reviews
#' 28(4): 372-375.
#'
#' @examples
#' x <- arima.sim(model = list(ar = 0.5), n = 200, rand.gen = rnorm)
#' out <- htdv_boot(as.numeric(x), mean, R = 500, type = "stationary", seed = 1)
#' out$ci_percentile
#'
#' @export
htdv_boot <- function(x, statistic, R = 1999L,
                      type = c("circular", "stationary"),
                      block_length = "auto",
                      level = 0.95, seed = NULL) {
  type <- match.arg(type)
  .assert_numeric_vector(x, "x")
  if (!is.function(statistic)) {
    stop("'statistic' must be a function of x.", call. = FALSE)
  }
  .assert_probability(level, "level")
  if (!is.null(seed)) {
    old_state <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    } else NULL
    on.exit({
      if (!is.null(old_state)) {
        assign(".Random.seed", old_state, envir = globalenv())
      }
    }, add = TRUE)
    set.seed(seed)
  }
  n <- length(x)
  b_opt <- if (identical(block_length, "auto")) {
    .ppw_block_length(x)
  } else {
    if (!is.numeric(block_length) || block_length < 1L ||
        block_length != floor(block_length)) {
      stop("'block_length' must be 'auto' or a positive integer.",
           call. = FALSE)
    }
    as.integer(block_length)
  }
  t0 <- statistic(x)
  t_star <- matrix(NA_real_, nrow = R, ncol = length(t0))
  for (r in seq_len(R)) {
    xb <- if (type == "circular") {
      .circular_bootstrap_sample(x, b = b_opt)
    } else {
      .stationary_bootstrap_sample(x, p = 1 / b_opt)
    }
    t_star[r, ] <- as.numeric(statistic(xb))
  }
  q_lo <- (1 - level) / 2
  q_hi <- 1 - q_lo
  ci_perc <- apply(t_star, 2L, stats::quantile,
                   probs = c(q_lo, q_hi), names = FALSE,
                   na.rm = TRUE)
  ci_basic <- rbind(2 * t0 - ci_perc[2L, ], 2 * t0 - ci_perc[1L, ])
  dimnames(ci_basic) <- dimnames(ci_perc)
  sd_star <- apply(t_star, 2L, stats::sd, na.rm = TRUE)
  z_lo <- stats::qnorm(q_lo)
  z_hi <- stats::qnorm(q_hi)
  ci_student <- rbind(t0 + z_lo * sd_star, t0 + z_hi * sd_star)
  dimnames(ci_student) <- dimnames(ci_perc)
  list(t0 = t0, t_star = t_star, ci_percentile = ci_perc,
       ci_basic = ci_basic, ci_studentized = ci_student,
       block_length = b_opt, type = type, level = level)
}

.circular_bootstrap_sample <- function(x, b) {
  n <- length(x)
  nb <- ceiling(n / b)
  starts <- sample.int(n, size = nb, replace = TRUE)
  idx <- integer(nb * b)
  for (i in seq_len(nb)) {
    s <- starts[i]
    block_idx <- ((s - 1L) + seq_len(b) - 1L) %% n + 1L
    idx[((i - 1L) * b + 1L):(i * b)] <- block_idx
  }
  x[idx[seq_len(n)]]
}

.stationary_bootstrap_sample <- function(x, p) {
  n <- length(x)
  out <- integer(n)
  out[1L] <- sample.int(n, size = 1L)
  for (t in 2:n) {
    if (stats::runif(1L) < p) {
      out[t] <- sample.int(n, size = 1L)
    } else {
      out[t] <- out[t - 1L] %% n + 1L
    }
  }
  x[out]
}

.ppw_block_length <- function(x) {
  n <- length(x)
  u <- .demean(x)
  max_lag <- floor(min(3 * sqrt(n), n / 3))
  acv <- .auto_cov(u, lag_max = max_lag)
  rho <- acv / acv[1L]
  thr <- 2 * sqrt(log10(n) / n)
  kn_m <- .find_significant_lag(rho[-1L], threshold = thr)
  m_hat <- min(2 * kn_m, max_lag)
  k_vec <- seq(-m_hat, m_hat)
  lam <- .flat_top_kernel(k_vec / m_hat)
  g_hat <- sum(lam * abs(k_vec) * acv[abs(k_vec) + 1L])
  d_hat <- 2 * (acv[1L] + 2 * sum(lam[k_vec > 0] * acv[k_vec[k_vec > 0] + 1L]))^2
  b <- max(1, round((2 * g_hat^2 / d_hat)^(1 / 3) * n^(1 / 3)))
  min(b, n %/% 2)
}

.find_significant_lag <- function(rho_nonzero, threshold) {
  K_limit <- max(5L, round(log10(length(rho_nonzero) + 1L)))
  for (k in seq_along(rho_nonzero)) {
    window_end <- min(k + K_limit, length(rho_nonzero))
    if (all(abs(rho_nonzero[k:window_end]) < threshold)) {
      return(k)
    }
  }
  length(rho_nonzero)
}

.flat_top_kernel <- function(z) {
  az <- abs(z)
  ifelse(az <= 0.5, 1,
         ifelse(az <= 1, 2 * (1 - az), 0))
}
