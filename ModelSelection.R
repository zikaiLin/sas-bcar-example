# Selection of best model

MNIST_task = "3_8"
approach = "SIRGPAM"
abbrv_approach = "sv"
task_1 = 8
n_models = 3
dataset = "MNIST"
sample_size = 150


data_path = sprintf("./Data/%s_%s/", dataset, MNIST_task)
res_path = sprintf("./Data/%s_%s/%s/N%d",
                   dataset,
                   MNIST_task,
                   approach,
                   sample_size)
output_dir = file.path(res_path, "best_mod")
if(!dir.exists(output_dir)) dir.create(output_dir)


## ----load data---------------------------------------------------------------------------------------------------------------------------
for(i in 1:10){
  
  mnist_data_training = readRDS(file.path(data_path, sprintf("training_%d_%i.rds", sample_size,i)))
  training_y = as.integer(mnist_data_training$mnist_data_training_Y_sub == task_1)
  training_X = matrix(NA, nrow = nrow(mnist_data_training$mnist_data_training_X_sub), ncol = 28*28)
  
  
  for(s in 1:nrow(mnist_data_training$mnist_data_training_X_sub)){
    training_X[s,] = c(mnist_data_training$mnist_data_training_X_sub[s,,])
  }
  
  training_X = training_X/255
  
  mods = lapply(1:n_models, function(k){
    readRDS(file = file.path(res_path,sprintf("model_%d_%s_n%d_%d.rds", i, abbrv_approach, sample_size, k)))
  })
  
  mods_best_DIC = DIC_compare(mods,x_train = training_X, y_train = training_y)
  
  saveRDS(mods_best_DIC, file = file.path(output_dir, sprintf("model_%d_%s_n%d.rds", i, abbrv_approach, sample_size)))
}



