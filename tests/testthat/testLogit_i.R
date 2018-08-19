library(Rrcgev)
testthat::context("Logit_i shares")

testthat::test_that("logit class works as expected", {
  mean_u <- rep(2.0, 10)

  shares <- exp(mean_u) / (1 + sum(exp(mean_u)))

  test_shares <- Rrcgev:::logit_i_test_shares(mean_u)

  inc_val <- log(1 + sum(exp(mean_u)))

  test_inc_val <- Rrcgev:::logit_i_test_val(mean_u)

  jacobian <- shares %*% t(shares)
  diag(jacobian) <- shares * (1 - shares)

  test_jacobian <- Rrcgev:::logit_i_test_jac(mean_u)

  testthat::expect_equal(test_jacobian, jacobian)
  testthat::expect_equal(test_inc_val, inc_val)
  testthat::expect_equal(test_shares, shares)
})
