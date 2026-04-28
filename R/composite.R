#' Composite Log-Likelihood for Dependent Data
#'
#' Pairwise or \code{k}-tuple composite log-likelihood (Varin, Reid and Firth
#' 2011) for strictly stationary sequences.
#'
#' @param x Numeric vector.
#' @param density_fun Function \code{function(block, theta)} returning the
#'   joint density of a \code{(k+1)}-block of consecutive observations.
#' @param theta Numeric parameter vector.
#' @param k Window width (pairwise if \code{k = 1}).
#' @param log Logical; if \code{TRUE}, \code{density_fun} already returns log
#'   densities.
#'
#' @return A list with \code{loglik} (scalar), \code{n_blocks}, and
#'   \code{block_values}.
#'
#' @references
#' Varin, C., Reid, N., & Firth, D. (2011). An overview of composite
#' likelihood methods. Statistica Sinica, 21(1), 5-42.
#'
#' @examples
#' x <- arima.sim(model = list(ar = 0.4), n = 200, rand.gen = rnorm)
#' gauss_ar1 <- function(block, theta) {
#'   phi <- theta[1]; sigma2 <- theta[2]
#'   a <- block[1]; b <- block[2]
#'   dnorm(a, 0, sqrt(sigma2 / (1 - phi^2)), log = TRUE) +
#'     dnorm(b, phi * a, sqrt(sigma2), log = TRUE)
#' }
#' htdv_composite(as.numeric(x), gauss_ar1, c(0.4, 1), k = 1, log = TRUE)$loglik
#'
#' @export
htdv_composite <- function(x, density_fun, theta, k = 1L, log = FALSE) {
  .assert_numeric_vector(x, "x")
  if (!is.function(density_fun)) {
    stop("'density_fun' must be a function of (block, theta).", call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || k < 1L || k != floor(k)) {
    stop("'k' must be a positive integer.", call. = FALSE)
  }
  n <- length(x)
  if (n < k + 1L) {
    stop("'x' is too short for the requested window 'k'.", call. = FALSE)
  }
  vals <- numeric(n - k)
  for (t in seq_len(n - k)) {
    block <- x[t:(t + k)]
    v <- density_fun(block, theta)
    vals[t] <- if (isTRUE(log)) v else log(max(v, .Machine$double.eps))
  }
  list(loglik = sum(vals), n_blocks = n - k, block_values = vals)
}
