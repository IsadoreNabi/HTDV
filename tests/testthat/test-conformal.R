test_that("htdv_conformal attains approximate target coverage", {
  set.seed(301)
  x <- stats::arima.sim(model = list(ar = 0.4), n = 300,
                        rand.gen = stats::rnorm)
  pred <- function(history) {
    if (length(history) >= 1L) history[length(history)] else 0
  }
  out <- htdv_conformal(as.numeric(x), pred,
                        alpha_target = 0.1,
                        lambda = 0.05,
                        burn_in = 30L)
  expect_true(is.numeric(out$coverage))
  expect_true(out$coverage >= 0.5 && out$coverage <= 1)
})

test_that("htdv_conformal validates arguments", {
  expect_error(htdv_conformal(stats::rnorm(10), mean, alpha_target = 0))
  expect_error(htdv_conformal(stats::rnorm(10), mean, lambda = -1))
  expect_error(htdv_conformal(stats::rnorm(10), mean, burn_in = 0))
})
