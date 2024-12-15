subject_id = 179
poly_degree = 20
a = 0.010000
b = 10
n_chains = 1
burnin = 2000
thinning = 1
mcmc_sample = 2000
setwd(".")
subject = sprintf("K%d", subject_id)

library(BayesGPfit)
library(RcppTN)
library(mgcv)
library(glmnet)
library(loo)
# poly_degree = 10
# a = 0.01
# b = 2
# burnin = 1000    
# thinning = 1
# mcmc_sample = 1000

## ----hyper-parameters--------------------------------------------------------------------------------------------------------------------

kernel_params_str = sprintf("poly%d_b%d", poly_degree, b)


loo_best = Inf


Rcpp::sourceCpp("./scripts/BasisExpansion_GP_BCI_analysis.cpp")
# Rcpp::sourceCpp("/Volumes/easystore/Dropbox/U-Mich/Research/Scalar-on-image regression via GP Additive Model/Scripts/Zikai/update_pj_mat.cpp")



dat = readRDS(sprintf("./BCI/data/%s.rds",subject))




training_X = dat$X_train
training_y = dat$y_train

training_X = training_X/max(abs(training_X))

exclude_column = as.integer(apply(training_X, 2, function(x) all(x == 0) ))
non_zero_column = 1-exclude_column


training_theta_1 = matrix(rnorm((poly_degree + 1) * length(non_zero_column), 0, 1),
                          (poly_degree + 1),
                          length(non_zero_column))




fit_glmnnet = glmnet::glmnet(training_X, training_y, family = "binomial",alpha = 1,)
cv.fit=cv.glmnet(training_X,training_y, alpha = 1, lambda = NULL, family = "binomial", nfolds = 5, type.measure = "mse")
beta_coef = as.numeric(coef(fit_glmnnet, s = cv.fit$lambda.min))
training_alpha_1 = beta_coef[1]

training_z_1 = rep(0, nrow(training_X))
training_z_1[training_y == 1] = 1
training_z_1[training_y == 0] = 0
training_delta_1 = as.integer(beta_coef != 0)[-1]

training_a_1 = 1 / rgamma(1, 0.5, 1)


prior_selection_prob = apply(training_X, 2, var)
partition = prior_selection_prob
partition[prior_selection_prob > quantile(prior_selection_prob,0.5)] = 1
partition[prior_selection_prob <= quantile(prior_selection_prob,0.5)] = 2


pj_mat = rep(0.5, ncol(training_X))


channels = rep(1:16, each = 26)
timepoints = rep(1:26, 16)
chanel2ind = matrix(1:(26*16), nrow = 26, ncol = 16)
order_bci = c('F3', 'Fz', 'F4', 'T7', 'C3', 'Cz', 'C4', 'T8', 'CP3', 'CP4', 'P3', 'Pz', 'P4', 'PO7', 'PO8', 'Oz')
dist_mat = read.csv("./BCI/distance_mat.csv")[,-1]
rownames(dist_mat) = colnames(dist_mat)

model_list = list()
## ----------------------------------------------------------------------------------------------------------------------------------------
for (s in 1:n_chains){
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
    channels2tp = chanel2ind-1,
    mcmc_sample = mcmc_sample,
    threshold = 0.6,
    initial_sigma_sq = 1,
    radius_candidats = c(2,4,6),
    verbose = 50,
    radius_update_thinning = 100,
    excluded_vox = exclude_column,
    channels = channels,
    timepoints = timepoints,
    n_timepoints = 26,
    n_channels = 16,
    prior_max = 1,
    prior_min = 0,
    initial_prob_vec = pj_mat,
    time_length_candidates = c(1,2),
    dist_mat = as.matrix(dist_mat)
  )
  model_best = model_1
}


filename = sprintf("./BCI_result/%s_mcmc.rds", subject)
if(!dir.exists(file.path("./BCI_result", kernel_params_str))) 
  dir.create(file.path("./BCI_result/", kernel_params_str), recursive = T)
saveRDS(model_best, sprintf("./BCI_result/%s/%s_mcmc.rds", kernel_params_str, subject))



SIRGPAM_prediction = function(mod, testing_X, delta_post = NULL) {
  require(BayesGPfit)

  delta_post = mod$post_mean$delta

  predict_testing_y_model_1 <- rep(mod$post_mean$alpha, nrow(testing_X))

  for (j in 1:ncol(testing_X)) {
    cat(j, " ")
    temp <- as.matrix(testing_X[, j])
    B <- GP.eigen.funcs.fast(temp,
                             poly_degree = poly_degree,
                             a = a,
                             b = b)

    theta_post = mod$post_mean$theta

    predict_testing_y_model_1 = predict_testing_y_model_1 + delta_post[j] *
      B %*% theta_post[, , j]

  }


  prob_predict <- pnorm(predict_testing_y_model_1)
  predict_testing_y_model_1 <- as.integer(pnorm(predict_testing_y_model_1)>0.7)





  return(list(prob_predict = prob_predict,
              label_predict = predict_testing_y_model_1) )
}



testing_X = dat$X_test  
testing_y = dat$y_test
testing_X = testing_X/max(abs(dat$X_train))
test = SIRGPAM_prediction(model_best, testing_X = testing_X)
# test_pred = apply(test, 2, function(x) as.integer(mean(x)>=0.5) )
print(mean(test$label_predict == testing_y))
print(caret::confusionMatrix(factor(test$label_predict), factor(testing_y)))

if(!dir.exists(file.path("./BCI_prediction/", kernel_params_str))) dir.create(file.path("./BCI_prediction/", kernel_params_str), recursive = T)
saveRDS(test, sprintf("./BCI_prediction/%s/K%d_prediction.rds", kernel_params_str, subject_id))
