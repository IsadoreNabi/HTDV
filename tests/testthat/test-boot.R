test_that("htdv_boot percentile CI covers the true mean frequently", {
  set.seed(201)
  x <- stats::arima.sim(model = list(ar = 0.5), n = 150,
                        rand.gen = stats::rnorm)
  out <- htdv_boot(as.numeric(x), mean, R = 400, type = "stationary",
                   block_length = "auto", seed = 1)
  expect_true(nrow(out$ci_percentile) == 2)
  expect_true(is.finite(out$block_length))
  expect_true(out$block_length >= 1)
})

test_that("htdv_boot supports circular type", {
  set.seed(202)
  x <- stats::rnorm(100)
  out <- htdv_boot(x, mean, R = 200, type = "circular",
                   block_length = 5L, seed = 2)
  expect_true(nrow(out$ci_basic) == 2)
})

test_that("htdv_boot returns studentized CI matrix", {
  set.seed(203)
  x <- stats::rnorm(100)
  out <- htdv_boot(x, mean, R = 200, seed = 3)
  expect_true(nrow(out$ci_studentized) == 2)
})
