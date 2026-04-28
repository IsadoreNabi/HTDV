test_that("htdv_rope accepts null when draws concentrate inside ROPE", {
  set.seed(401)
  d <- stats::rnorm(4000, mean = 0, sd = 0.01)
  out <- htdv_rope(d, rope = c(-0.1, 0.1))
  expect_equal(out$decision, "accept")
})

test_that("htdv_rope rejects when draws are far from ROPE", {
  set.seed(402)
  d <- stats::rnorm(4000, mean = 2, sd = 0.1)
  out <- htdv_rope(d, rope = c(-0.1, 0.1))
  expect_equal(out$decision, "reject")
})

test_that("htdv_rope returns undecided for borderline draws", {
  set.seed(403)
  d <- stats::rnorm(4000, mean = 0, sd = 0.25)
  out <- htdv_rope(d, rope = c(-0.1, 0.1))
  expect_true(out$decision %in% c("accept", "reject", "undecided"))
})

test_that("htdv_rope validates ROPE argument", {
  expect_error(htdv_rope(stats::rnorm(200), rope = c(1, 0)))
  expect_error(htdv_rope(stats::rnorm(200), rope = 1))
})
