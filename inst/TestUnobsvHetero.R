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



