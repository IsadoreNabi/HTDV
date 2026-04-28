test_that("htdv_whittle returns finite log-likelihood for AR(1)", {
  set.seed(501)
  x <- stats::arima.sim(model = list(ar = 0.5), n = 200,
                        rand.gen = stats::rnorm)
  ar1_spec <- function(lambda, theta) {
    phi <- theta[1]
    sigma2 <- theta[2]
    sigma2 / (2 * pi * (1 - 2 * phi * cos(lambda) + phi^2))
  }
  out <- htdv_whittle(as.numeric(x), ar1_spec, c(0.5, 1))
  expect_true(is.finite(out$loglik))
  expect_equal(length(out$frequencies), length(out$periodogram))
})

test_that("htdv_whittle rejects non-positive spectrum", {
  x <- stats::rnorm(50)
  bad_spec <- function(lambda, theta) rep(-1, length(lambda))
  out <- htdv_whittle(x, bad_spec, theta = 1)
  expect_true(is.infinite(out$loglik))
})

test_that("htdv_composite evaluates pairwise Gaussian correctly", {
  set.seed(502)
  x <- stats::arima.sim(model = list(ar = 0.3), n = 150,
                        rand.gen = stats::rnorm)
  gauss_ar1 <- function(block, theta) {
    phi <- theta[1]
    sigma2 <- theta[2]
    a <- block[1]
    b <- block[2]
    stats::dnorm(a, 0, sqrt(sigma2 / (1 - phi^2)), log = TRUE) +
      stats::dnorm(b, phi * a, sqrt(sigma2), log = TRUE)
  }
  out <- htdv_composite(as.numeric(x), gauss_ar1, c(0.3, 1),
                        k = 1L, log = TRUE)
  expect_true(is.finite(out$loglik))
  expect_equal(out$n_blocks, length(x) - 1L)
})
