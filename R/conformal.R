#' Adaptive Conformal Inference for Dependent Data
#'
#' Online adaptive conformal prediction intervals with gradient update on the
#' miscoverage rate (Gibbs-Candes 2021). Provides long-run marginal coverage
#' without exchangeability.
#'
#' @param x Numeric vector of sequential observations.
#' @param predictor Function \code{function(history)} returning a scalar
#'   one-step point prediction.
#' @param residual_fun Function of residuals to calibrate the score; default
#'   absolute residual.
#' @param alpha_target Target miscoverage rate.
#' @param lambda Step size for the gradient update.
#' @param burn_in Integer; the first \code{burn_in} observations are used only
#'   to warm up the predictor.
#'
#' @return A list with \code{intervals} (matrix of lower/upper),
#'   \code{coverage} (empirical), \code{alpha_path}, and
#'   \code{point_predictions}.
#'
#' @references
#' Gibbs, I., & Candes, E. (2021). Adaptive Conformal Inference under
#' Distribution Shift. NeurIPS 34: 1660-1672.
#'
#' @examples
#' x <- arima.sim(model = list(ar = 0.6), n = 200, rand.gen = rnorm)
#' pred <- function(hist) if (length(hist) >= 1) hist[length(hist)] else 0
#' out <- htdv_conformal(as.numeric(x), pred, alpha_target = 0.1,
#'                       lambda = 0.05, burn_in = 20)
#' out$coverage
#'
#' @export
htdv_conformal <- function(x, predictor,
                           residual_fun = function(e) abs(e),
                           alpha_target = 0.1,
                           lambda = 0.05,
                           burn_in = 20L) {
  .assert_numeric_vector(x, "x")
  if (!is.function(predictor)) {
    stop("'predictor' must be a function of the history vector.",
         call. = FALSE)
  }
  .assert_probability(alpha_target, "alpha_target")
  if (!is.numeric(lambda) || length(lambda) != 1L || lambda <= 0 || lambda >= 1) {
    stop("'lambda' must be a scalar in (0, 1).", call. = FALSE)
  }
  n <- length(x)
  if (burn_in < 1L || burn_in >= n) {
    stop("'burn_in' must be an integer in [1, length(x)-1].", call. = FALSE)
  }
  alpha_t <- alpha_target
  scores <- numeric(0L)
  yhat <- numeric(n)
  intervals <- matrix(NA_real_, nrow = n, ncol = 2L)
  alpha_path <- numeric(n)
  colnames(intervals) <- c("lower", "upper")
  for (t in seq_len(n)) {
    hist_t <- x[seq_len(t - 1L)]
    yhat[t] <- if (t > 1L) predictor(hist_t) else 0
    alpha_path[t] <- alpha_t
    if (t > burn_in && length(scores) > 0L) {
      q <- as.numeric(stats::quantile(scores,
                                      probs = min(1, max(0, 1 - alpha_t)),
                                      names = FALSE, na.rm = TRUE))
      intervals[t, ] <- c(yhat[t] - q, yhat[t] + q)
      covered <- (x[t] >= intervals[t, 1L]) && (x[t] <= intervals[t, 2L])
      alpha_t <- alpha_t + lambda * (alpha_target - (1 - as.integer(covered)))
      alpha_t <- min(0.99, max(0.01, alpha_t))
    }
    if (t > 1L) {
      scores <- c(scores, as.numeric(residual_fun(x[t] - yhat[t])))
    }
  }
  valid <- !is.na(intervals[, 1L])
  cov_emp <- if (any(valid)) {
    mean((x[valid] >= intervals[valid, 1L]) &
         (x[valid] <= intervals[valid, 2L]))
  } else NA_real_
  list(intervals = intervals, coverage = cov_emp,
       alpha_path = alpha_path, point_predictions = yhat,
       alpha_target = alpha_target, lambda = lambda)
}
