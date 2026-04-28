#' Fixed-Bandwidth HAR Inference (Kiefer-Vogelsang 2005)
#'
#' Fixed-b HAR Wald test where the bandwidth is a fixed fraction of the sample
#' size. Critical values are obtained by simulation from the asymptotic
#' functional distribution.
#'
#' @param x Numeric vector.
#' @param theta0 Null value for the mean (default 0).
#' @param B Fixed bandwidth fraction in (0, 1].
#' @param kernel Kernel choice (see \code{\link{htdv_lrv}}).
#' @param sims Number of Monte Carlo draws for the critical-value simulation.
#' @param alpha Significance level.
#' @param seed Optional integer seed for reproducibility. NULL leaves the
#'   user's RNG state untouched.
#'
#' @return A list with \code{statistic}, \code{critical_value}, \code{p_value},
#'   \code{reject}, \code{bandwidth}, and \code{kernel}.
#'
#' @references
#' Kiefer, N.M., & Vogelsang, T.J. (2005). A New Asymptotic Theory for
#' Heteroskedasticity-Autocorrelation Robust Tests. Econometric Theory
#' 21(6): 1130-1164.
#'
#' @examples
#' x <- arima.sim(model = list(ar = 0.5), n = 200, rand.gen = rnorm)
#' htdv_fixedb(as.numeric(x), theta0 = 0, B = 0.2, sims = 1000, seed = 1)$p_value
#'
#' @export
htdv_fixedb <- function(x, theta0 = 0,
                        B = 0.1,
                        kernel = c("bartlett", "parzen", "qs"),
                        sims = 5000L, alpha = 0.05,
                        seed = NULL) {
  kernel <- match.arg(kernel)
  .assert_numeric_vector(x, "x")
  if (!is.numeric(B) || length(B) != 1L || B <= 0 || B > 1) {
    stop("'B' must be a scalar in (0, 1].", call. = FALSE)
  }
  .assert_probability(alpha, "alpha")
  n <- length(x)
  b <- max(floor(B * n), 1L)
  u <- x - theta0
  mu_hat <- mean(x)
  resid <- x - mu_hat
  lrv <- .kernel_sum(.auto_cov(resid, lag_max = n - 1L),
                     h = b, kernel = kernel, n = n)
  stat <- n * (mu_hat - theta0)^2 / lrv
  null_draws <- .simulate_fixedb_null(n = n, B = B, kernel = kernel,
                                      sims = sims, seed = seed)
  cv <- as.numeric(stats::quantile(null_draws, probs = 1 - alpha,
                                   names = FALSE))
  pv <- mean(null_draws >= stat)
  list(statistic = stat, critical_value = cv, p_value = pv,
       reject = stat > cv, bandwidth = b, kernel = kernel)
}

.simulate_fixedb_null <- function(n, B, kernel, sims, seed = NULL) {
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
  b <- max(floor(B * n), 1L)
  out <- numeric(sims)
  for (s in seq_len(sims)) {
    z <- stats::rnorm(n)
    mu <- mean(z)
    r <- z - mu
    lrv <- .kernel_sum(.auto_cov(r, lag_max = n - 1L),
                       h = b, kernel = kernel, n = n)
    out[s] <- n * mu^2 / lrv
  }
  out
}
