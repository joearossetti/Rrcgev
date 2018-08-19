library(Rrcgev)

vec <- 1:10
start_vec <- c(0, 5)
end_vec <- c(4,9)
sum1 <- sum(vec[(start_vec[1]+1):(end_vec[1]+1)])
sum2 <- sum(vec[(start_vec[2]+1):(end_vec[2]+1)])

vec_sub(x=vec, start_vec = start_vec, end_vec = end_vec, M = 2)
