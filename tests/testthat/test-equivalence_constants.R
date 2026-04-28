test_that("htdv_equivalence_constants returns positive K_TAC and K_MPC", {
  out <- htdv_equivalence_constants(gamma = 2, q = 6, n = 100)
  expect_true(out$K_TAC > 0)
  expect_true(out$K_MPC > 0)
  expect_true(out$c_U > 0)
})

test_that("higher gamma reduces K_TAC", {
  a <- htdv_equivalence_constants(gamma = 2, q = 6, n = 100)
  b <- htdv_equivalence_constants(gamma = 4, q = 6, n = 100)
  expect_lt(b$K_TAC, a$K_TAC)
})

test_that("larger n tightens c_L", {
  small <- htdv_equivalence_constants(gamma = 4, q = 6, n = 200)
  big   <- htdv_equivalence_constants(gamma = 4, q = 6, n = 50000)
  expect_gt(big$c_L, small$c_L)
  expect_true(big$c_L > 0)
})

test_that("htdv_equivalence_constants validates inputs", {
  expect_error(htdv_equivalence_constants(gamma = 0.5))
  expect_error(htdv_equivalence_constants(gamma = 2, q = 1))
  expect_error(htdv_equivalence_constants(gamma = 1.1, q = 2.01))
})
