library(readxl)
library(tidyverse)

# number of trials from 1:9
subject_id = 179

# process example data
# data$pred_sgpam = c(data$Swlda) c(gp_pred$prob_predict)
character_accuracy = function(data, pred, k){
  data$pred_sgpam = pred
  data%>%
    rename(Id = `Character Id`,
           predict = pred_sgpam)->
    data
  
  # find true characters
  data %>%
    filter(true == 1) %>%
    group_by(Id) %>%
    summarise(true =  paste0(sort(unique(Code)), collapse = '-'))->
    true_characters
  # find predicted characters
  data%>%
    group_by(Id)%>%
    filter(row_number() / 12 <= k)%>%
    group_by(Id, Code)%>%
    summarise(predict = sum(predict))%>%
    ungroup()%>%
    mutate(row = Code <= 6)%>%
    group_by(Id, row)%>%
    filter(predict == max(predict))%>%
    ungroup()%>%
    group_by(Id)%>%
    summarise(predict =  paste0(sort(unique(Code)), collapse = '-'))->
    predicted_characters
  
  true_characters%>%
    left_join(predicted_characters, by = "Id")%>%
    summarise(accuracy = mean(true == predict))%>%
    pull(accuracy) ->
    accuracy
  
  accuracy
}



char_acc = data_frame(k  = 1:9,
                      char_acc_swlda = rep(NA, 9),
                      char_acc_gp = rep(NA, 9),
                      char_acc_dnn = rep(NA, 9),
                      char_acc_tloho = rep(NA, 9))

for(k in 1:9){
  data <- read_excel(sprintf("./BCI/K%d.xlsx", subject_id), sheet = "y_test")
  gp_pred <- readRDS(sprintf("./BCI_results/BCI_prediction_best/K%d_prediction.rds", subject_id))
  dnn_pred <- readRDS(sprintf("./BCI_results/DNN/model_BCI_DNN_n.rds"))$pred
  tloho_pred <- readRDS(sprintf("./BCI_results/TLOHO/ypred_tloto.rds"))
  
  # Accuracy of Swlda
  char_acc_swlda = character_accuracy(data, c(data$Swlda), k = k)
  
  # Accuracy of gp_pred
  char_acc_gp = character_accuracy(data, c(gp_pred$prob_predict), k = k)
  
  # Accuracy of gp_pred
  char_acc_tloho = character_accuracy(data, c(tloho_pred), k = k)
  
  # Accuracy of gp_pred
  char_acc_dnn = character_accuracy(data, c(dnn_pred), k = k)
  
  
  saveRDS(list(char_acc_swlda = char_acc_swlda,
               char_acc_gp = char_acc_gp,
               char_acc_tloho = char_acc_tloho,
               char_acc_dnn = char_acc_dnn),
          sprintf("./BCI_results/BCI_char_acc/K%d_char_acc.rds", subject_id))
  
  char_acc[k, "char_acc_swlda"] = char_acc_swlda
  char_acc[k, "char_acc_gp"] = char_acc_gp
  char_acc[k, "char_acc_tloho"] = char_acc_tloho
  char_acc[k, "char_acc_dnn"] = char_acc_dnn
  
}


library(ggpubr)
dat_to_plot = data.frame(
  K = factor(as.integer(rep(1:9, 4))),
  Accuracy = c(char_acc$char_acc_swlda, char_acc$char_acc_gp, char_acc$char_acc_tloho, char_acc$char_acc_dnn),
  Method = factor(rep(c("Swalda", "SAS-BCAR", "TLOHO", "DNN"), each = 9))
)
ggpubr::ggline(
  dat_to_plot,
  x = "K",
  y = "Accuracy",
  group = "Method",
  color = "Method",
  numeric.x.axis = F,
  palette = "aaas",
  title = "Character Prediction Accuracy",
  xlab = "Number of Sequences",
  ylab = "Character-based Prediction Accuracy",
  size = 1,
  point.size = 2,
  label.rectangle = T)
