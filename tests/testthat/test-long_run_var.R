test_that("htdv_lrv returns a positive LRV for iid data", {
  set.seed(17)
  x <- stats::rnorm(200)
  out <- htdv_lrv(x, kernel = "qs", bandwidth = "andrews")
  expect_true(is.numeric(out$lrv))
  expect_gt(out$lrv, 0)
  expect_true(is.finite(out$lrv))
})

test_that("htdv_lrv inflates LRV for positively correlated data", {
  set.seed(18)
  x_ar <- stats::arima.sim(model = list(ar = 0.8), n = 500,
                           rand.gen = stats::rnorm)
  out_ar <- htdv_lrv(as.numeric(x_ar), kernel = "bartlett",
                     bandwidth = "andrews")
  set.seed(19)
  out_iid <- htdv_lrv(stats::rnorm(500), kernel = "bartlett",
                      bandwidth = "andrews")
  expect_gt(out_ar$lrv, out_iid$lrv)
})

test_that("htdv_lrv rejects invalid inputs", {
  expect_error(htdv_lrv(c(1, 2, NA)))
  expect_error(htdv_lrv(1))
  expect_error(htdv_lrv(stats::rnorm(50), bandwidth = -1))
})

test_that("htdv_lrv supports all three kernels", {
  set.seed(20)
  x <- stats::rnorm(100)
  for (k in c("bartlett", "parzen", "qs")) {
    out <- htdv_lrv(x, kernel = k, bandwidth = "andrews")
    expect_true(out$lrv > 0)
  }
})
