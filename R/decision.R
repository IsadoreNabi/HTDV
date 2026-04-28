#' Region-of-Practical-Equivalence (ROPE) Decision
#'
#' Posterior decision under Kruschke's (2018) Region of Practical Equivalence.
#'
#' @param draws Numeric vector of posterior draws of the parameter of
#'   interest.
#' @param rope Length-2 numeric vector \code{c(lower, upper)} defining the
#'   practical-equivalence region.
#' @param level Credibility level for the highest-density interval.
#'
#' @return A list with \code{decision} (\code{"accept"}, \code{"reject"},
#'   \code{"undecided"}), \code{hdi}, \code{rope}, and \code{prob_in_rope}.
#'
#' @references
#' Kruschke, J.K. (2018). Rejecting or Accepting Parameter Values in Bayesian
#' Estimation. Advances in Methods and Practices in Psychological Science
#' 1(2): 270-280.
#'
#' @examples
#' draws <- rnorm(4000, mean = 0.05, sd = 0.1)
#' htdv_rope(draws, rope = c(-0.1, 0.1), level = 0.95)$decision
#'
#' @export
htdv_rope <- function(draws, rope, level = 0.95) {
  .assert_numeric_vector(draws, "draws", min_len = 100L)
  if (!is.numeric(rope) || length(rope) != 2L || rope[1L] >= rope[2L]) {
    stop("'rope' must be a length-2 numeric vector with rope[1] < rope[2].",
         call. = FALSE)
  }
  .assert_probability(level, "level")
  hdi <- .hdi(draws, level = level)
  inside <- (hdi[1L] >= rope[1L]) && (hdi[2L] <= rope[2L])
  outside <- (hdi[2L] < rope[1L]) || (hdi[1L] > rope[2L])
  decision <- if (inside) "accept" else if (outside) "reject" else "undecided"
  prob_in <- mean(draws >= rope[1L] & draws <= rope[2L])
  list(decision = decision, hdi = hdi, rope = rope,
       level = level, prob_in_rope = prob_in)
}

.hdi <- function(draws, level = 0.95) {
  s <- sort(draws)
  n <- length(s)
  k <- max(1L, floor(level * n))
  if (k >= n) return(c(s[1L], s[n]))
  widths <- s[(k + 1L):n] - s[1L:(n - k)]
  i <- which.min(widths)
  c(s[i], s[i + k])
}

#' Bridge-Sampling Bayes Factor
#'
#' Computes the Bayes factor between two fitted \code{\link{htdv_fit}}
#' objects via bridge sampling. Requires the \code{bridgesampling} package.
#'
#' @param fit1 Model in the numerator.
#' @param fit0 Model in the denominator.
#' @param ... Additional arguments passed to
#'   \code{bridgesampling::bridge_sampler}.
#'
#' @return A list with \code{bf10}, \code{log_bf10}, \code{logml1},
#'   \code{logml0}, and \code{error_percentage}.
#'
#' @references
#' Meng, X.-L., & Wong, W.H. (1996). Simulating Ratios of Normalizing
#' Constants. Statistica Sinica 6(4): 831-860.
#' Gronau, Q.F. et al. (2017). A Tutorial on Bridge Sampling. JMP 81: 80-97.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("bridgesampling", quietly = TRUE)) {
#'   x <- rnorm(50)
#'   f1 <- htdv_fit(x, model = "tac", chains = 2, iter = 500,
#'                  refresh = 0, seed = 1)
#'   f0 <- htdv_fit(x, model = "tac", chains = 2, iter = 500,
#'                  refresh = 0, seed = 2)
#'   htdv_bf(f1, f0)$bf10
#' }
#' }
#'
#' @export
htdv_bf <- function(fit1, fit0, ...) {
  .require_suggested("bridgesampling")
  if (!inherits(fit1, "htdv_fit") || !inherits(fit0, "htdv_fit")) {
    stop("'fit1' and 'fit0' must be 'htdv_fit' objects.", call. = FALSE)
  }
  bs1 <- bridgesampling::bridge_sampler(fit1$stanfit, silent = TRUE, ...)
  bs0 <- bridgesampling::bridge_sampler(fit0$stanfit, silent = TRUE, ...)
  log_bf <- bs1$logml - bs0$logml
  list(bf10 = exp(log_bf), log_bf10 = log_bf,
       logml1 = bs1$logml, logml0 = bs0$logml,
       error_percentage = c(bridgesampling::error_measures(bs1)$percentage,
                            bridgesampling::error_measures(bs0)$percentage))
}

#' WAIC and Leave-Future-Out Cross-Validation
#'
#' Computes the Widely Applicable Information Criterion (Watanabe 2010) and
#' leave-future-out expected log predictive density (Buerkner-Gabry-Vehtari
#' 2020) from a fitted \code{\link{htdv_fit}} object.
#'
#' @param fit An \code{htdv_fit} object with pointwise log-likelihood stored
#'   as parameter \code{log_lik}.
#' @param L Block-prefix length; leave-future-out starts at time \code{L+1}.
#'   Default \code{floor(n/5)}.
#' @param k_threshold Pareto-k refresh threshold; when exceeded, importance
#'   sampling is refreshed.
#'
#' @return A list with \code{waic}, \code{elpd_lfo}, \code{elpd_se},
#'   \code{pareto_k}, and \code{refresh_times}.
#'
#' @references
#' Watanabe, S. (2010). Asymptotic Equivalence of Bayes Cross Validation and
#' Widely Applicable Information Criterion. JMLR 11: 3571-3594.
#' Buerkner, P.-C., Gabry, J., & Vehtari, A. (2020). Approximate
#' Leave-Future-Out Cross-Validation for Bayesian Time Series Models. Journal
#' of Statistical Computation and Simulation 90(14): 2499-2523.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("loo", quietly = TRUE)) {
#'   x <- rnorm(50)
#'   fit <- htdv_fit(x, model = "tac", chains = 2, iter = 500,
#'                   refresh = 0, seed = 1)
#'   htdv_waic_lfo(fit)$waic
#' }
#' }
#'
#' @export
htdv_waic_lfo <- function(fit, L = NULL, k_threshold = 0.7) {
  .require_suggested("loo")
  if (!inherits(fit, "htdv_fit")) {
    stop("'fit' must be an 'htdv_fit' object.", call. = FALSE)
  }
  ll <- .extract_log_lik(fit)
  if (is.null(ll)) {
    stop("Log-likelihood 'log_lik' not found in the Stan fit.", call. = FALSE)
  }
  n <- ncol(ll)
  if (is.null(L)) L <- max(1L, floor(n / 5L))
  w <- loo::waic(ll)
  elpd_lfo <- NA_real_
  elpd_se <- NA_real_
  pareto_k <- rep(NA_real_, n)
  refreshes <- integer(0L)
  if (L < n) {
    elpds <- numeric(n - L)
    for (t in (L + 1L):n) {
      lw <- .log_weights_lfo(ll, t, L)
      pk <- .pareto_k_fast(lw)
      pareto_k[t] <- pk
      if (!is.finite(pk) || pk > k_threshold) {
        refreshes <- c(refreshes, t)
        lw <- rep(0, length(lw))
      }
      elpds[t - L] <- .log_mean_exp(ll[, t] + lw - .log_sum_exp(lw))
    }
    elpd_lfo <- sum(elpds)
    elpd_se <- sqrt(length(elpds)) * stats::sd(elpds)
  }
  list(waic = w, elpd_lfo = elpd_lfo, elpd_se = elpd_se,
       pareto_k = pareto_k, refresh_times = refreshes)
}

.extract_log_lik <- function(fit) {
  sf <- fit$stanfit
  pars <- sf@model_pars
  if (!"log_lik" %in% pars) return(NULL)
  as.matrix(rstan::extract(sf, pars = "log_lik", permuted = TRUE)$log_lik)
}

.log_weights_lfo <- function(log_lik, t, L) {
  if (t - L <= 1L) return(numeric(nrow(log_lik)))
  rowSums(log_lik[, (L + 1L):(t - 1L), drop = FALSE])
}

.pareto_k_fast <- function(lw) {
  w <- exp(lw - max(lw))
  n <- length(w)
  tail_n <- max(5L, floor(n / 5L))
  s <- sort(w, decreasing = TRUE)[seq_len(tail_n)]
  if (any(s <= 0)) return(NA_real_)
  .gpd_shape(s)
}

.gpd_shape <- function(x) {
  m <- length(x)
  if (m < 5L) return(NA_real_)
  x <- sort(x)
  prior <- 3
  n_grid <- 30L
  bs <- 1 / x[m] + (1 - sqrt(m / (seq_len(n_grid) - 0.5))) /
    (prior * x[floor(m / 4 + 0.5)])
  ks <- sapply(bs, function(b) {
    k <- mean(log1p(-b * x))
    m * (log(-b / k) - k - 1)
  })
  w <- 1 / sapply(bs, function(b) {
    sum(exp(ks - max(ks)) * (1 / (1e-9 + abs(b - bs))))
  })
  b_hat <- sum(bs * w) / sum(w)
  k_hat <- mean(log1p(-b_hat * x))
  k_hat
}

.log_mean_exp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

.log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

#' Predictive Stacking of Bayesian Models
#'
#' Computes stacking weights (Yao, Vehtari, Simpson, Gelman 2018) across a
#' list of fitted models by LOO-log-score maximization.
#'
#' @param fits Named list of \code{htdv_fit} objects.
#' @param method Stacking objective. Only \code{"log_score"} supported in this
#'   release.
#'
#' @return A list with \code{weights}, \code{model_names}, and
#'   \code{elpd_per_model}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("loo", quietly = TRUE)) {
#'   x <- rnorm(50)
#'   f1 <- htdv_fit(x, model = "tac", chains = 2, iter = 500,
#'                  refresh = 0, seed = 1)
#'   f2 <- htdv_fit(x, model = "wsc", chains = 2, iter = 500,
#'                  refresh = 0, seed = 2)
#'   htdv_stack(list(tac = f1, wsc = f2))$weights
#' }
#' }
#'
#' @export
htdv_stack <- function(fits, method = c("log_score")) {
  method <- match.arg(method)
  .require_suggested("loo")
  if (!is.list(fits) || length(fits) < 2L) {
    stop("'fits' must be a list with at least two 'htdv_fit' objects.",
         call. = FALSE)
  }
  ll_list <- lapply(fits, .extract_log_lik)
  if (any(vapply(ll_list, is.null, logical(1L)))) {
    stop("All fits must contain 'log_lik'.", call. = FALSE)
  }
  ncols <- vapply(ll_list, ncol, integer(1L))
  if (length(unique(ncols)) != 1L) {
    stop("All fits must have 'log_lik' with identical number of observations; ",
         "got ", paste(ncols, collapse = ", "),
         ". Whittle and time-domain likelihoods cannot be stacked directly.",
         call. = FALSE)
  }
  loo_list <- lapply(ll_list, function(ll) loo::loo(ll, cores = 1L))
  lpd_point <- vapply(loo_list,
                      function(l) as.numeric(l$pointwise[, "elpd_loo"]),
                      numeric(ncols[[1L]]))
  if (!is.matrix(lpd_point)) lpd_point <- as.matrix(lpd_point)
  w <- .stacking_weights(lpd_point)
  list(weights = stats::setNames(w, names(fits)),
       model_names = names(fits),
       elpd_per_model = vapply(loo_list,
                               function(l) l$estimates["elpd_loo", 1L],
                               numeric(1L)))
}

.stacking_weights <- function(lpd) {
  K <- ncol(lpd)
  n <- nrow(lpd)
  neg_obj <- function(z) {
    w <- .simplex(z)
    -sum(log(pmax(lpd %*% w, .Machine$double.eps)))
  }
  start <- rep(0, K - 1L)
  opt <- stats::optim(start, neg_obj, method = "BFGS")
  .simplex(opt$par)
}

.simplex <- function(z) {
  e <- c(exp(z), 1)
  e / sum(e)
}
