#' Whittle Likelihood for Stationary Series
#'
#' Computes the Whittle pseudo-log-likelihood at Fourier frequencies for a
#' user-supplied spectral density.
#'
#' @param x Numeric vector.
#' @param spec_fun Function \code{function(lambda, theta)} returning spectral
#'   density values at angular frequencies \code{lambda}.
#' @param theta Numeric parameter vector passed to \code{spec_fun}.
#' @param taper Logical; if \code{TRUE}, a 10-percent cosine taper is applied.
#'
#' @return A list with \code{loglik} (scalar), \code{frequencies},
#'   \code{periodogram}, and \code{fitted_spectrum}.
#'
#' @references
#' Whittle, P. (1953). Estimation and Information in Stationary Time Series.
#' Arkiv foer Matematik, 2(5), 423-434.
#' Dzhaparidze, K. (1986). Parameter Estimation and Hypothesis Testing in
#' Spectral Analysis of Stationary Time Series. Springer.
#'
#' @examples
#' x <- arima.sim(model = list(ar = 0.5), n = 200, rand.gen = rnorm)
#' ar1_spec <- function(lambda, theta) {
#'   phi <- theta[1]; sigma2 <- theta[2]
#'   sigma2 / (2 * pi * (1 - 2 * phi * cos(lambda) + phi^2))
#' }
#' htdv_whittle(as.numeric(x), ar1_spec, c(0.5, 1))$loglik
#'
#' @export
htdv_whittle <- function(x, spec_fun, theta, taper = FALSE) {
  .assert_numeric_vector(x, "x")
  if (!is.function(spec_fun)) {
    stop("'spec_fun' must be a function of (lambda, theta).", call. = FALSE)
  }
  n <- length(x)
  u <- .demean(x)
  if (isTRUE(taper)) u <- .cosine_taper(u, p = 0.1)
  pg <- .periodogram(u)
  lambda <- pg$lambda
  inp <- pg$intensity
  fv <- spec_fun(lambda, theta)
  if (any(!is.finite(fv)) || any(fv <= 0)) {
    return(list(loglik = -Inf, frequencies = lambda, periodogram = inp,
                fitted_spectrum = fv))
  }
  loglik <- -sum(log(fv) + inp / fv)
  list(loglik = loglik, frequencies = lambda, periodogram = inp,
       fitted_spectrum = fv)
}

.periodogram <- function(u) {
  n <- length(u)
  m <- floor(n / 2)
  j <- seq_len(m)
  lambda <- 2 * pi * j / n
  fft_u <- stats::fft(u)
  inp <- (1 / (2 * pi * n)) * Mod(fft_u[j + 1L])^2
  list(lambda = lambda, intensity = inp)
}

.cosine_taper <- function(u, p = 0.1) {
  n <- length(u)
  m <- floor(p * n / 2)
  if (m < 1) return(u)
  w <- rep(1, n)
  seq_edge <- seq_len(m)
  w[seq_edge] <- 0.5 * (1 - cos(pi * (seq_edge - 0.5) / m))
  w[(n - m + 1):n] <- rev(w[seq_edge])
  u * w
}
