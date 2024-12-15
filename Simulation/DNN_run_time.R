
library(keras)
setwd("~/Documents/BCAR/")
dataset = "MNIST_3_8"
task_num = 8
out_dir = sprintf("./Data/%s/DNN", dataset)
if(!dir.exists(out_dir))  dir.create(out_dir, recursive = T)

sample_sizes = c(50, 150, 1000)
computation_times_dnn = c(NA, NA, NA)
mnist_data_testing = readRDS(sprintf("./Data/%s/testing.rds", dataset))
mnist_data_testing_X = mnist_data_testing$mnist_data_testing_X
mnist_data_testing_Y = as.integer(mnist_data_testing$mnist_data_testing_Y == task_num)

i = 1
for (t in 1:length(sample_sizes)){
  sample_size = sample_sizes[t]
  model <- keras_model_sequential() %>%
    layer_flatten(input_shape = c(28, 28)) %>%
    layer_dense(units = 21, activation = 'relu') %>%
    layer_dense(units = 2, activation = 'sigmoid')
  
  model %>% compile(
    optimizer = 'adam',
    loss = 'sparse_categorical_crossentropy',
    metrics = c('accuracy')
  )
  
  mnist_data_training = readRDS(sprintf(".//Data/%s/training_%d/training_%d_%d.rds", dataset, sample_size, sample_size, i))
  mnist_data_training_X = mnist_data_training$mnist_data_training_X_sub
  mnist_data_training_Y = as.integer(mnist_data_training$mnist_data_training_Y_sub == task_num)
  
  time_taken = system.time({
    model %>% fit(mnist_data_training_X, mnist_data_training_Y, epochs = 5, verbose = 0)
  })
  
  computation_times_dnn[t] = time_taken
}
saveRDS(computation_times_dnn, file = "./Data/MNIST_3_8/TLOTO/computation_times_dnn.rds")
