test_that("htdv_fit returns an htdv_fit object for a minimal TAC run", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  set.seed(601)
  x <- stats::rnorm(40)
  fit <- suppressWarnings(
    htdv_fit(x, model = "tac", chains = 1L, iter = 200L,
             warmup = 100L, refresh = 0L, seed = 1L))
  expect_s3_class(fit, "htdv_fit")
  expect_equal(fit$model, "tac")
  expect_equal(fit$n, length(x))
})

test_that("htdv_envelope combines draws across models", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  set.seed(602)
  x <- stats::rnorm(40)
  f1 <- suppressWarnings(
    htdv_fit(x, model = "tac", chains = 1L, iter = 200L,
             warmup = 100L, refresh = 0L, seed = 2L))
  f2 <- suppressWarnings(
    htdv_fit(x, model = "composite", chains = 1L, iter = 200L,
             warmup = 100L, refresh = 0L, seed = 3L))
  env <- htdv_envelope(list(f1, f2), target = "theta")
  expect_true(length(env$draws) > 0L)
  expect_length(env$intervals, 2L)
})
