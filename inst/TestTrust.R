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

index_i_test <- DemandLogit__compute(share_test, my_delta)

share_bench <- microbenchmark(
  index_i_test <- DemandLogit__compute(share_test, my_delta)
)

share_tests <- DemandLogit__getShares(share_test)
inc_val_tests <- DemandLogit__getIncValue(index_i_test)

abs(my_shares - share_tests) <= .Machine$double.eps

f_maker <- function(s_0, mu_mat, u_opt_out){
  dl_ptr <- DemandLogit__new(mu_mat, u_opt_out)
  f_ <- function(x){
    res <- list()
    ## compute shares
    DemandLogit__compute(dl_ptr, x)

    s <- DemandLogit__getShares(dl_ptr)
    V <- DemandLogit__getIncValue(dl_ptr)

    res$value <- V - t(s_0) %*% x
    res$gradient <- s - s_0
    temp <- s %*% t(s)
    diag(temp) <- s*(1-s)
    res$hessian <- temp
    return(res)
  }
  return(f_)
}

my_f_ <- f_maker(s_0 = my_shares, mu_mat = my_mu, u_opt_out = 0.0)

my_f_(rep(1.0,10))
my_f_(my_delta)

r_test <- microbenchmark::microbenchmark(
  trust_test_result <- trust::trust(my_f_, parinit = rep(1.0,10), rinit = 0.1, rmax = 10.0, fterm = .Machine$double.eps,  blather = TRUE)
)

f_maker_cpp <- function(s_0, mu_mat, u_opt_out){
  dlr_ptr <- DemandLogitRoot__new(s_0, mu_mat, u_opt_out)
  f_ <- function(x){
    DemandLogitRoot__compute_obj(dlr_ptr, x)
    res <- list()
    res$value <- DemandLogitRoot__getValue(dlr_ptr)
    res$gradient <- DemandLogitRoot__getGradient(dlr_ptr)
    res$hessian <- DemandLogitRoot__getHessian(dlr_ptr)
    return(res)
  }
  return(f_)
}

my_f_cpp <- f_maker_cpp(s_0 = my_shares, mu_mat = my_mu, u_opt_out = 0.0)

test <- my_f_(rep(1.0,10))
test_cpp <- my_f_cpp(rep(1.0,10))

abs(test$value - test_cpp$value) <= .Machine$double.eps
abs(test$gradient - test_cpp$gradient) <= .Machine$double.eps
abs(test$hessian - test$hessian) <= .Machine$double.eps

cpp_test <- microbenchmark::microbenchmark(
  trust_test_result <- trust::trust(my_f_cpp, parinit = rep(1.0,10), rinit = 0.1, rmax = 10.0, fterm = .Machine$double.eps,  blather = TRUE)
)

invert_test <- DemandInvertTrust_Logit__new(delta_init = rep(1.0, 10), s_0 = my_shares, mu_=my_mu, u_opt_out_ = 0.0,
                             max_radius = 10.0, thresh = c(0.125,0.25,0.75),
                             scale = c(0.25,2.0), tol = 1e-3, max_count = 1000)

test_deltas <- DemandInvertTrust_Logit__root_cauchy(invert_test)
