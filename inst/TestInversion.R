library(Rrcgev)
library(microbenchmark)

set.seed(1234)
my_num_draws <- 1000
my_dim_sigma <- 1
my_rand_mat <- matrix(rnorm(my_num_draws*my_dim_sigma, 0, 1), my_dim_sigma, my_num_draws)
i_mat <- diag(my_dim_sigma)

my_which_nln_prmtrs <- as.vector(i_mat) == 1

test <- UnobsvHetero__new(my_rand_mat, my_which_nln_prmtrs,my_num_draws, my_dim_sigma)

test_mat <- UnobsvHetero__getNln_parmtrs_mat(test, rep(1.0, 5))

test_mat2 <- UnobsvHetero__getNui_mat(test, test_mat)

J <- 10

alpha <- -4
beta <- 2

prods_chars <- sample(c(1,0), J, replace = TRUE)
prices <- rep(1, J) + rnorm(J, mean = 0, sd = 0.25)

my_mu <- t(t(prices)) %*% (alpha + test_mat2)
my_delta <- beta * prods_chars

my_index_i <- matrix(my_delta, J, my_num_draws, byrow = FALSE) + my_mu
max_index <- apply(my_index_i, MARGIN = 2, FUN = max)
my_index_i2 <- my_index_i - matrix(max_index, J, my_num_draws, byrow = TRUE)
exp_index_i <- exp(my_index_i2)

s_i <- exp_index_i / matrix((colSums(exp_index_i) + exp(0-max_index)), J, my_num_draws, byrow = TRUE)
my_shares = rowMeans(s_i)

share_test <- DemandLogit__new(my_mu, 0.0)

index_i_test <- DemandLogit__shares(share_test, my_delta)

share_bench <- microbenchmark(
  index_i_test <- DemandLogit__shares(share_test, my_delta)
)

abs(my_shares - index_i_test) <= .Machine$double.eps

test_delta_hat_fp <- DemandLogit__invert_fp(share_test, index_i_test, delta_0 = rep(1.0, J), tol = 1.0e-15, max_count = 1000)

test_delta_hat_nr <- DemandLogit__invert_nr(share_test, index_i_test, delta_0 = rep(1.0, J), tol = 1.0e-15, max_count = 1000)

test_delta_hat_nrb <- DemandLogit__invert_nrb(share_test, index_i_test, delta_0 = rep(1.0, J), tol = 1.0e-15, max_count = 1000)

test_delta_hat_nrl <- DemandLogit__invert_nrl(share_test, index_i_test, delta_0 = rep(1.0, J), tol = 1.0e-15, max_count = 1000)

test_delta_hat_fpj <- DemandLogit__invert_fpj(share_test, index_i_test, delta_0 = rep(1.0, J), tol = 1.0e-15, max_count = 1000)


fp_bench <- microbenchmark(
  test_delta_hat_fp <- DemandLogit__invert_fp(share_test, index_i_test, delta_0 = rep(1.0, J), tol = 1.0e-15, max_count = 1000)
)

fpj_bench <- microbenchmark(
  test_delta_hat_fpj <- DemandLogit__invert_fpj(share_test, index_i_test, delta_0 = rep(1.0, J), tol = 1.0e-15, max_count = 1000)
)

nr_bench <- microbenchmark(
  test_delta_hat_nr <- DemandLogit__invert_nr(share_test, index_i_test, delta_0 = rep(1.0, J), tol = 1.0e-15, max_count = 1000)
)

nrb_bench <- microbenchmark(
  test_delta_hat_nrb <- DemandLogit__invert_nrb(share_test, index_i_test, delta_0 = rep(1.0, J), tol = 1.0e-15, max_count = 1000)
)

nrl_bench <- microbenchmark(
  test_delta_hat_nrl <- DemandLogit__invert_nrl(share_test, index_i_test, delta_0 = rep(1.0, J), tol = 1.0e-15, max_count = 1000)
)

# close_delta <- DemandLogit__invert_fp(share_test, index_i_test, delta_0 = rep(1.0, J), tol = 1.0e-8, max_count = 1000)
#
# cfp_bench <- microbenchmark(
#   test_delta_hat_fp <- DemandLogit__invert_fp(share_test, index_i_test, delta_0 = close_delta, tol = 1.0e-15, max_count = 1000)
# )
#
# cnr_bench <- microbenchmark(
#   test_delta_hat_nr <- DemandLogit__invert_nr(share_test, index_i_test, delta_0 = close_delta, tol = 1.0e-15, max_count = 1000)
# )
#
# cnrb_bench <- microbenchmark(
#   test_delta_hat_nrb <- DemandLogit__invert_nrb(share_test, index_i_test, delta_0 = close_delta, tol = 1.0e-15, max_count = 1000)
# )
