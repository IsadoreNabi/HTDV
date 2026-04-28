#' Long-Run Variance Estimation (HAC)
#'
#' Heteroskedasticity-and-autocorrelation-consistent long-run variance with
#' data-driven Andrews (1991) optimal bandwidth, Bartlett, Parzen, and
#' Quadratic-Spectral kernels.
#'
#' @param x Numeric vector of observations.
#' @param kernel Kernel choice: \code{"bartlett"}, \code{"parzen"}, or
#'   \code{"qs"} (Quadratic-Spectral).
#' @param bandwidth Either \code{"andrews"} (data-driven) or a positive
#'   numeric scalar.
#' @param prewhiten Logical; if \code{TRUE}, AR(1) prewhiten-recolor
#'   (Andrews-Monahan 1992).
#'
#' @return A list with components \code{lrv} (long-run variance estimate),
#'   \code{bandwidth}, \code{kernel}, and \code{ar_coef} (only if prewhitened).
#'
#' @references
#' Andrews, D.W.K. (1991). Heteroskedasticity and Autocorrelation Consistent
#' Covariance Matrix Estimation. Econometrica 59(3): 817-858.
#' Newey, W.K. & West, K.D. (1987). A Simple, Positive Semi-Definite
#' Heteroskedasticity and Autocorrelation Consistent Covariance Matrix.
#' Econometrica 55(3): 703-708.
#'
#' @examples
#' set_seed_user <- 42
#' x <- arima.sim(model = list(ar = 0.6), n = 200, rand.gen = rnorm)
#' out <- htdv_lrv(as.numeric(x), kernel = "qs", bandwidth = "andrews")
#' out$lrv
#'
#' @export
htdv_lrv <- function(x, kernel = c("bartlett", "parzen", "qs"),
                     bandwidth = "andrews", prewhiten = FALSE) {
  kernel <- match.arg(kernel)
  .assert_numeric_vector(x, "x")
  n <- length(x)
  u <- .demean(x)
  ar_coef <- NA_real_
  if (isTRUE(prewhiten)) {
    fit <- stats::ar(u, aic = FALSE, order.max = 1L, method = "ols")
    if (length(fit$ar) >= 1L) {
      ar_coef <- as.numeric(fit$ar[1L])
      u <- as.numeric(stats::na.omit(fit$resid))
      n <- length(u)
    }
  }
  h <- if (identical(bandwidth, "andrews")) {
    .andrews_bandwidth(u, kernel = kernel)
  } else {
    .assert_positive_scalar(bandwidth, "bandwidth")
    as.numeric(bandwidth)
  }
  acv <- .auto_cov(u, lag_max = n - 1L)
  lrv <- .kernel_sum(acv, h = h, kernel = kernel, n = n)
  if (isTRUE(prewhiten) && is.finite(ar_coef) && abs(ar_coef) < 0.999) {
    lrv <- lrv / (1 - ar_coef)^2
  }
  list(lrv = lrv, bandwidth = h, kernel = kernel, ar_coef = ar_coef)
}

.auto_cov <- function(u, lag_max) {
  n <- length(u)
  lag_max <- min(lag_max, n - 1L)
  out <- numeric(lag_max + 1L)
  for (k in 0:lag_max) {
    out[k + 1L] <- sum(u[1:(n - k)] * u[(1 + k):n]) / n
  }
  out
}

.kernel_weight <- function(z, kernel) {
  if (kernel == "bartlett") {
    pmax(0, 1 - abs(z))
  } else if (kernel == "parzen") {
    az <- abs(z)
    w <- ifelse(az <= 0.5, 1 - 6 * az^2 + 6 * az^3,
                ifelse(az <= 1, 2 * (1 - az)^3, 0))
    w
  } else if (kernel == "qs") {
    out <- numeric(length(z))
    nz <- abs(z) > .Machine$double.eps
    zz <- z[nz]
    arg <- 6 * pi * zz / 5
    out[nz] <- (25 / (12 * pi^2 * zz^2)) *
      (sin(arg) / arg - cos(arg))
    out[!nz] <- 1
    out
  } else {
    stop("Unknown kernel.", call. = FALSE)
  }
}

.kernel_sum <- function(acv, h, kernel, n) {
  lag_max <- length(acv) - 1L
  k_vec <- seq_len(lag_max)
  z <- k_vec / h
  w <- .kernel_weight(z, kernel)
  val <- acv[1L] + 2 * sum(w * acv[-1L])
  max(val, .Machine$double.eps)
}

.andrews_bandwidth <- function(u, kernel) {
  n <- length(u)
  fit <- stats::ar(u, aic = FALSE, order.max = 1L, method = "ols")
  rho <- if (length(fit$ar) >= 1L) as.numeric(fit$ar[1L]) else 0
  sig2 <- if (length(fit$var.pred) >= 1L) as.numeric(fit$var.pred) else .safe_var(u)
  if (!is.finite(rho) || abs(rho) >= 0.999) rho <- 0.99 * sign(rho)
  num <- 4 * rho^2 * sig2^2
  den <- (1 - rho)^6 * (1 + rho)^2
  alpha2 <- num / max(den, .Machine$double.eps)
  num4 <- 4 * rho^2 * sig2^2
  den4 <- (1 - rho)^8
  alpha4 <- num4 / max(den4, .Machine$double.eps)
  h <- switch(kernel,
    bartlett = 1.1447 * (alpha2 * n)^(1 / 3),
    parzen   = 2.6614 * (alpha4 * n)^(1 / 5),
    qs       = 1.3221 * (alpha4 * n)^(1 / 5)
  )
  max(h, 1)
}
