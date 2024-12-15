# TLOTO

library(TLOHO)
library(fields)
setwd("~/Documents/BCAR/")
sample_sizes = c(50, 150, 1000)
computation_times = rep(NA, length(sample_sizes))

for(t in 1:length(sample_sizes)){
  sample_size = sample_sizes[t]
  mnist_data_training = readRDS(sprintf("./Data/MNIST_3_8/training_%d/training_%d_1.rds", sample_size, sample_size))
  training_y = as.integer(mnist_data_training$mnist_data_training_Y_sub == 7)
  training_X = matrix(NA, nrow = nrow(mnist_data_training$mnist_data_training_X_sub), ncol = 28*28)
  for(s in 1:nrow(mnist_data_training$mnist_data_training_X_sub)){
    training_X[s,] = c(mnist_data_training$mnist_data_training_X_sub[s,,])
  }
  training_X = training_X/255
  
  
  mnist_data_testing = readRDS("./Data/MNIST_3_8/testing.rds")
  testing_X = matrix(NA, nrow = nrow(mnist_data_testing$mnist_data_testing_X), ncol = 28*28)
  testing_y = as.integer(mnist_data_testing$mnist_data_testing_Y == 7)
  
  for(s in 1:nrow(mnist_data_testing$mnist_data_testing_X)){
    testing_X[s,] = c(mnist_data_testing$mnist_data_testing_X[s,,])
  }
  testing_X = testing_X/255
  
  out_dir = "./Data/MNIST_3_8/TLOTO"
  if(!dir.exists(out_dir))  dir.create(out_dir)
  
  

  i = 1
    ## ----load data---------------------------------------------------------------------------------------------------------------------------
    mnist_data_training = readRDS(sprintf("./Data/MNIST_3_8/training_%d/training_%d_1.rds", sample_size, sample_size))
    training_y = as.integer(mnist_data_training$mnist_data_training_Y_sub == 8)
    training_X = matrix(NA, nrow = nrow(mnist_data_training$mnist_data_training_X_sub), ncol = 28*28)
    for(s in 1:nrow(mnist_data_training$mnist_data_training_X_sub)){
      training_X[s,] = c(mnist_data_training$mnist_data_training_X_sub[s,,])
    }
    training_X = training_X/255
    
    
    ## ----find non-zero-column----------------------------------------------------------------------------------------------------------------
    #all elements in training_X >= 0
    non_zero_column = as.integer(apply(training_X, 2, function(x) (sum(x)>0)))
    exclude_column = 1- non_zero_column
    
    
    non_zero_column = as.logical(as.integer(apply(training_X, 2, sum) > 0))
    time_taken = system.time({
      graph0 = igraph::make_lattice(c(28,28))# construct lattice graph
      graph1 = igraph::delete_vertices(graph0, v = which(!non_zero_column))
      fit <- tloho_lm(training_y,
                      training_X[,non_zero_column],
                      graph0 =  graph1,
                      nsave = 2000,
                      nburn = 2000,
                      nthin = 1) # fit
      
    })

    
    
    beta_est = rep(0, 784)
    beta_est[non_zero_column] = fit$median_beta_est
    image.plot(matrix(beta_est, 28, 28))
    
    
    y_pred = pnorm(testing_X %*% beta_est)
    y_pred = as.integer(y_pred>0.5)
    print(mean(y_pred == testing_y))
    computation_times[t] = time_taken
    # saveRDS(list(model = fit, y_pred = y_pred, testing_accuracy = mean(y_pred == testing_y)), 
            # file.path(out_dir, sprintf("model_%d_tloto_n50_2.rds", i)))
    
  
}

saveRDS(computation_times, file = "./Data/MNIST_3_8/TLOTO/computation_times_tloho.rds")
