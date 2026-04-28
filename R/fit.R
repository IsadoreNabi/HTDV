#' Hierarchical Bayesian Fit for Dependent Variables
#'
#' Fits the hierarchical Bayesian model described in the companion paper via
#' Hamiltonian Monte Carlo (NUTS) using \pkg{rstan}. Nuisance parameters
#' \code{sigma_inf}, \code{gamma_mix}, and (for WSC) \code{b} carry priors and
#' are integrated over.
#'
#' @param x Numeric vector of observations.
#' @param model One of \code{"tac"}, \code{"wsc"}, \code{"mpc"},
#'   \code{"whittle"}, \code{"composite"}.
#' @param chains Number of HMC chains.
#' @param iter Total iterations per chain (warmup plus sampling).
#' @param warmup Number of warmup iterations.
#' @param refresh Stan progress refresh rate (0 suppresses output).
#' @param adapt_delta HMC adaptation target acceptance.
#' @param max_treedepth Maximum NUTS treedepth.
#' @param seed Optional integer seed forwarded to \pkg{rstan}.
#' @param ... Additional arguments forwarded to \code{rstan::sampling}.
#'
#' @return An object of class \code{htdv_fit}, a list with components
#'   \code{stanfit}, \code{model}, \code{data}, \code{n}, \code{call}.
#'
#' @references
#' Hoffman, M.D., & Gelman, A. (2014). The No-U-Turn Sampler. JMLR 15(1):
#' 1593-1623.
#' Betancourt, M. (2016). Diagnosing Suboptimal Cotangent Disintegrations in
#' HMC. <doi:10.48550/arXiv.1604.00695>.
#'
#' @examples
#' \donttest{
#' x <- rnorm(50)
#' fit <- htdv_fit(x, model = "tac", chains = 2, iter = 500,
#'                 refresh = 0, seed = 1)
#' class(fit)
#' }
#'
#' @export
htdv_fit <- function(x,
                     model = c("tac", "wsc", "mpc", "whittle", "composite"),
                     chains = 4L, iter = 4000L, warmup = iter %/% 2L,
                     refresh = max(iter %/% 10L, 1L),
                     adapt_delta = 0.99, max_treedepth = 12L,
                     seed = NULL, ...) {
  model <- match.arg(model)
  .assert_numeric_vector(x, "x")
  stan_name <- switch(model,
    tac       = "tac_hierarchical",
    wsc       = "wsc_hierarchical",
    mpc       = "mpc_hierarchical",
    whittle   = "whittle_hierarchical",
    composite = "composite_hierarchical"
  )
  sm <- .load_stan_model(stan_name)
  stan_data <- .build_stan_data(x = x, model = model)
  control <- list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  args <- list(object = sm, data = stan_data, chains = chains,
               iter = iter, warmup = warmup, refresh = refresh,
               control = control, ...)
  if (!is.null(seed)) args$seed <- seed
  stanfit <- do.call(rstan::sampling, args)
  out <- list(stanfit = stanfit, model = model, data = stan_data,
              n = length(x), call = match.call())
  class(out) <- "htdv_fit"
  out
}

.build_stan_data <- function(x, model) {
  n <- length(x)
  lrv <- htdv_lrv(x, kernel = "qs", bandwidth = "andrews")$lrv
  base <- list(N = n, x = as.numeric(x),
               sigma_inf_scale = max(2 * sqrt(lrv), .Machine$double.eps),
               gamma_shape = 2, gamma_rate = 1,
               mu0 = mean(x),
               tau0 = max(10 * stats::sd(x), .Machine$double.eps))
  if (model == "wsc") {
    base$b_lower <- max(2L, floor(0.5 * sqrt(n)))
    base$b_upper <- max(base$b_lower + 1L, ceiling(2 * sqrt(n)))
  }
  base
}

#' @rdname htdv_fit
#' @export
htdv_tac <- function(x, ...) htdv_fit(x = x, model = "tac", ...)

#' @rdname htdv_fit
#' @export
htdv_wsc <- function(x, ...) htdv_fit(x = x, model = "wsc", ...)

#' @rdname htdv_fit
#' @export
htdv_mpc <- function(x, ...) htdv_fit(x = x, model = "mpc", ...)

#' Berger-Robust Envelope Posterior
#'
#' Builds an envelope posterior over a list of \code{htdv_fit} objects by
#' entropy-maximizing mixture weights subject to predictive consistency
#' (Berger 1994).
#'
#' @param fits List of \code{htdv_fit} objects on the same data.
#' @param target Name of the Stan parameter over which to build the envelope.
#'
#' @return A list with \code{draws} (envelope draws for the target parameter),
#'   \code{weights}, \code{intervals}, and \code{component_intervals}.
#'
#' @references
#' Berger, J.O. (1994). An Overview of Robust Bayesian Analysis. Test 3(1):
#' 5-124.
#'
#' @examples
#' \donttest{
#' x <- rnorm(50)
#' f1 <- htdv_fit(x, model = "tac", chains = 2, iter = 500,
#'                refresh = 0, seed = 1)
#' f2 <- htdv_fit(x, model = "wsc", chains = 2, iter = 500,
#'                refresh = 0, seed = 2)
#' htdv_envelope(list(f1, f2), target = "theta")$intervals
#' }
#'
#' @export
htdv_envelope <- function(fits, target = "theta") {
  if (!is.list(fits) || !all(vapply(fits, inherits, logical(1L), "htdv_fit"))) {
    stop("'fits' must be a list of 'htdv_fit' objects.", call. = FALSE)
  }
  draws_list <- lapply(fits, function(f) {
    d <- rstan::extract(f$stanfit, pars = target, permuted = TRUE)[[1L]]
    as.numeric(d)
  })
  lengths <- vapply(draws_list, length, integer(1L))
  if (length(unique(lengths)) > 1L) {
    lo <- min(lengths)
    draws_list <- lapply(draws_list, function(d) d[seq_len(lo)])
  }
  ent <- vapply(draws_list, function(d) {
    sd_d <- stats::sd(d)
    if (!is.finite(sd_d) || sd_d <= 0) 0 else log(sd_d)
  }, numeric(1L))
  w <- exp(ent - max(ent))
  w <- w / sum(w)
  k <- length(draws_list)
  n_env <- min(vapply(draws_list, length, integer(1L)))
  pick <- sample.int(k, size = n_env, replace = TRUE, prob = w)
  envelope <- vapply(seq_len(n_env), function(i) draws_list[[pick[i]]][i],
                     numeric(1L))
  comp_ci <- t(vapply(draws_list, function(d) {
    as.numeric(stats::quantile(d, probs = c(0.025, 0.975), names = FALSE))
  }, numeric(2L)))
  env_ci <- as.numeric(stats::quantile(envelope, probs = c(0.025, 0.975),
                                       names = FALSE))
  list(draws = envelope, weights = w,
       intervals = env_ci, component_intervals = comp_ci)
}

#' @export
print.htdv_fit <- function(x, ...) {
  cat("HTDV fit\n")
  cat("  model:     ", x$model, "\n", sep = "")
  cat("  n:         ", x$n, "\n", sep = "")
  cat("  chains:    ", x$stanfit@sim$chains, "\n", sep = "")
  cat("  iters:     ", x$stanfit@sim$iter, "\n", sep = "")
  invisible(x)
}

#' @export
summary.htdv_fit <- function(object, pars = NULL, probs = c(0.025, 0.5, 0.975),
                             ...) {
  s <- if (is.null(pars)) {
    rstan::summary(object$stanfit, probs = probs)
  } else {
    rstan::summary(object$stanfit, pars = pars, probs = probs)
  }
  s$summary
}

#' @export
plot.htdv_fit <- function(x, pars = NULL, ...) {
  if (is.null(pars)) pars <- "theta"
  rstan::traceplot(x$stanfit, pars = pars, ...)
}

#' Print-Friendly Summary for an HTDV Fit
#'
#' @param fit An \code{htdv_fit} object.
#' @param rope Optional two-element vector for ROPE decision on the target.
#' @param target Parameter name to summarize.
#'
#' @return A data frame with posterior mean, sd, 2.5\%, 50\%, 97.5\%, Rhat,
#'   bulk and tail ESS for \code{target}, plus the ROPE decision if requested.
#'
#' @examples
#' \donttest{
#' x <- rnorm(50)
#' fit <- htdv_fit(x, model = "tac", chains = 2, iter = 500,
#'                 refresh = 0, seed = 1)
#' htdv_summary(fit, target = "theta")
#' }
#'
#' @export
htdv_summary <- function(fit, rope = NULL, target = "theta") {
  if (!inherits(fit, "htdv_fit")) {
    stop("'fit' must be an 'htdv_fit' object.", call. = FALSE)
  }
  s <- rstan::summary(fit$stanfit, pars = target,
                      probs = c(0.025, 0.5, 0.975))$summary
  df <- as.data.frame(s)
  if (!is.null(rope)) {
    draws <- as.numeric(rstan::extract(fit$stanfit, pars = target,
                                       permuted = TRUE)[[1L]])
    dec <- htdv_rope(draws, rope = rope)
    df$rope_decision <- dec$decision
  }
  df
}
