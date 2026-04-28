.assert_numeric_vector <- function(x, name = "x", min_len = 2L) {
  if (!is.numeric(x) || !is.null(dim(x))) {
    stop(sprintf("'%s' must be a numeric vector.", name), call. = FALSE)
  }
  if (length(x) < min_len) {
    stop(sprintf("'%s' must have length >= %d.", name, min_len), call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(sprintf("'%s' must not contain NA/NaN/Inf.", name), call. = FALSE)
  }
  invisible(TRUE)
}

.assert_positive_scalar <- function(x, name = "x") {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(sprintf("'%s' must be a positive finite scalar.", name), call. = FALSE)
  }
  invisible(TRUE)
}

.assert_probability <- function(p, name = "p") {
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p <= 0 || p >= 1) {
    stop(sprintf("'%s' must be a scalar in (0, 1).", name), call. = FALSE)
  }
  invisible(TRUE)
}

.has_package <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

.require_suggested <- function(pkg) {
  if (!.has_package(pkg)) {
    stop(sprintf(
      "Package '%s' is required for this function. Please install it.",
      pkg
    ), call. = FALSE)
  }
  invisible(TRUE)
}

.with_restored_options <- function(expr) {
  old <- options()
  on.exit(options(old), add = TRUE)
  force(expr)
}

.load_stan_model <- function(name) {
  cache <- get("stan_cache", envir = .htdv_env)
  if (!is.null(cache[[name]])) return(cache[[name]])
  stan_file <- system.file("stan", paste0(name, ".stan"), package = "HTDV")
  if (!nzchar(stan_file)) {
    stop(sprintf("Stan file '%s.stan' not found in HTDV/inst/stan.", name),
         call. = FALSE)
  }
  model <- rstan::stan_model(file = stan_file, auto_write = FALSE)
  cache[[name]] <- model
  assign("stan_cache", cache, envir = .htdv_env)
  model
}

.demean <- function(x) x - mean(x)

.safe_var <- function(x) {
  v <- stats::var(x)
  if (!is.finite(v) || v <= 0) 1e-12 else v
}
