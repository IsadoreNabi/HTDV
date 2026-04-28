test_that("htdv_fixedb returns a valid p-value on iid data", {
  set.seed(101)
  x <- stats::rnorm(150)
  out <- htdv_fixedb(x, theta0 = 0, B = 0.1,
                     sims = 500L, seed = 1)
  expect_true(out$p_value >= 0 && out$p_value <= 1)
  expect_true(is.finite(out$statistic))
})

test_that("htdv_fixedb tends to reject when H0 is violated", {
  set.seed(102)
  x <- stats::rnorm(200, mean = 2)
  out <- htdv_fixedb(x, theta0 = 0, B = 0.1,
                     sims = 500L, seed = 2)
  expect_true(out$reject)
})

test_that("htdv_fixedb checks argument validity", {
  expect_error(htdv_fixedb(stats::rnorm(50), B = 0))
  expect_error(htdv_fixedb(stats::rnorm(50), B = 1.5))
  expect_error(htdv_fixedb(stats::rnorm(50), alpha = -0.1))
})
