

## ----hyper-parameters--------------------------------------------------------------------------------------------------------------------
poly_degree = 20
a = 0.01
b = 100
burnin = 2000
approach = "SIRGPAM"
approach_abrv = "sv"
thinning = 1
mcmc_sample = 2000
setwd("~/Documents/BCAR/")
data_path = "./Data/MNIST_3_8/"
subtask = 8

## ---------- Sample size  ---------- ##
sample_sizes = c(50, 150, 1000)
computation_times = rep(NA, length(sample_sizes))

## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------
for(s in 1:length(sample_sizes)){
sample_size = sample_sizes[s]

training_prefix = sprintf("training_%d", sample_size)
library(BayesGPfit)
library(RcppTN)
library(mgcv)
library(ggplot2)



Rcpp::sourceCpp("./Scripts/BasisExpansion_GP.cpp")
Rcpp::sourceCpp("./Scripts/Zikai/update_pj_mat.cpp")




mnist_data_testing = readRDS(file.path(data_path, "testing.rds"))
testing_X = matrix(NA,
                   nrow = nrow(mnist_data_testing$mnist_data_testing_X),
                   ncol = 28 * 28)
testing_y = as.integer(mnist_data_testing$mnist_data_testing_Y == subtask)

for (s in 1:nrow(mnist_data_testing$mnist_data_testing_X)) {
  testing_X[s, ] = c(mnist_data_testing$mnist_data_testing_X[s, , ])
}
testing_X = testing_X / 255



i = 1

cat("---------------- i = ", i, "----------------\n")
## ----load data---------------------------------------------------------------------------------------------------------------------------
mnist_data_training = readRDS(file.path(
  data_path,
  training_prefix,
  sprintf("%s_%i.rds", training_prefix, i)
))

training_y = as.integer(mnist_data_training$mnist_data_training_Y_sub == subtask)
training_X = matrix(
  NA,
  nrow = nrow(mnist_data_training$mnist_data_training_X_sub),
  ncol = 28 * 28
)
for (s in 1:nrow(mnist_data_training$mnist_data_training_X_sub)) {
  training_X[s, ] = c(mnist_data_training$mnist_data_training_X_sub[s, , ])
}
training_X = training_X / 255


## ----find non-zero-column----------------------------------------------------------------------------------------------------------------
#all elements in training_X >= 0
non_zero_column = as.integer(apply(training_X, 2, function(x)
  (sum(x) > 0)))
exclude_column = 1 - non_zero_column



training_theta_1 = matrix(rnorm((poly_degree + 1) * length(non_zero_column), 0, 1)
                          ,
                          (poly_degree + 1),
                          length(non_zero_column))
training_alpha_1 = rnorm(1, 0, 1)
training_z_1 = rep(training_alpha_1, nrow(training_X))
training_delta_1 = rbinom(ncol(training_X), 1, 0.5)
training_a_1 = 1 / rgamma(1, 0.5, 1)


for (j in 1:length(non_zero_column)) {
  temp <- as.matrix(training_X[, j])
  B <- GP.eigen.funcs.fast(temp,
                           poly_degree = poly_degree,
                           a = a,
                           b = b)
  training_z_1 <- training_z_1 + training_delta_1[j] * B %*% training_theta_1[, j]
}


prior_selection_prob = apply(training_X, 2, var)
partition = matrix(prior_selection_prob, 28, 28)
partition[prior_selection_prob > 0.05] = 1
partition[prior_selection_prob <= 0.05] = 2


pj_mat = update_pj_mat(matrix(training_delta_1, 28, 28),
                       partition,
                       radius_partition = rep(1, length(unique(partition))))
pj_mat[pj_mat == 1] <- 0.95
pj_mat[pj_mat == 0] <- 0.05


## ----------------------------------------------------------------------------------------------------------------------------------------
time_taken = system.time({
model_1 = Basis_Expansion(
  y = training_y,
  X = training_X,
  partition = partition,
  initial_alpha = training_alpha_1,
  initial_a = training_a_1,
  initial_delta = training_delta_1,
  initial_z = training_z_1 ,
  initial_theta = training_theta_1,
  poly_degree = poly_degree,
  poly_a = a,
  poly_b = b,
  oc_theta = T,
  oc_alpha = T,
  oc_delta = T,
  oc_z = T,
  oc_sigma = T,
  oc_a = T,
  burnin = burnin,
  thinning = thinning,
  mcmc_sample = mcmc_sample,
  threshold = 0.6,
  initial_sigma_sq = 1,
  initial_prob_mat = pj_mat,
  radius_candidats = c(1, 2),
  verbose = 500,
  radius_update_thinning = 100,
  excluded_vox = exclude_column,
  prior_max = 1,
  prior_min = 0
)
})



## ----training_result_model_1-------------------------------------------------------------------------------------------------------------
predict_training_p_model_1 <- rep(model_1$post_mean$alpha, nrow(training_X))
predict_training_y_model_1 <- as.numeric(rep(0, nrow(training_X)))
# delta_post = as.numeric(rep(0,ncol(training_X)))
# delta_post[non_zero_column == 1] = model_1$post_mean$delta
delta_post = model_1$post_mean$delta
# theta_post = array(0, dim = c(dim(model_1$post_mean$theta)[1:2], ncol(training_X)))
# theta_post[,,non_zero_column == 1] = model_1$post_mean$theta
theta_post = model_1$post_mean$theta
for (j in 1:ncol(training_X)) {
  # training_X_scale[,j] = training_X[,non_zero_column[j]]
  if (non_zero_column[j] == 1) {
    temp <- as.matrix(training_X[, j])
    B <- GP.eigen.funcs.fast(temp,
                             poly_degree = poly_degree,
                             a = a,
                             b = b)
    predict_training_p_model_1 <-  predict_training_p_model_1 + delta_post[j] *
      B %*% theta_post[, , j]
  }
  predict_training_y_model_1 <- rbinom(nrow(training_X), 1, pnorm(predict_training_p_model_1))
}

print(mean(training_y == predict_training_y_model_1))
computation_times[s] = time_taken[3]
}

saveRDS(computation_times, file.path(data_path, "computation_times_sasbcar.rds"))