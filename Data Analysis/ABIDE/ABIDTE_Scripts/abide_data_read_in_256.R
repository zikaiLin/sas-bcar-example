# Load required packages
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(tibble)
library(neurobase)
library(caret)
library(oro.nifti)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("./Scripts/find_supervoxel_neighbours.cpp")
source("./Scripts/compute_super_voxels_aal.R")

set.seed(20241111)

Phenotypic_V1_0b <- read.csv("./ABIDE/Phenotypic_V1_0b.csv", as.is = T)
# 1. Prepare phenotypic data
Phenotypic_V1_0b <- Phenotypic_V1_0b %>%
  mutate(
    SITE_ID = as.character(SITE_ID),
    SUB_ID = as.character(SUB_ID),
    ID = sprintf("%07d", as.numeric(SUB_ID)),  # Ensure SUB_ID is numeric
    SITE_ID = toupper(SITE_ID)  # Convert to uppercase for consistency
  )

# 2. List and parse image files
falff_dir <- "./ABIDE_Data/falff/"
image_files <- list.files(
  path = falff_dir,
  pattern = "_falff\\.nii\\.gz$",
  full.names = TRUE
)
image_filenames <- basename(image_files)

image_metadata <- tibble(
  filepath = image_files,
  filename = image_filenames
) %>%
  mutate(
    base_name = str_remove(filename, "_falff\\.nii\\.gz$"),
    SUB_ID = str_extract(base_name, "\\d{7}$"),
    SITE_ID = str_remove(base_name, "_\\d{7}$"),
    SITE_ID = toupper(SITE_ID)
  ) %>%
  select(filepath, filename, SITE_ID, SUB_ID)

# 3. Merge image metadata with phenotypic data
matched_data <- image_metadata %>%
  inner_join(
    Phenotypic_V1_0b,
    by = c("SITE_ID" = "SITE_ID", "SUB_ID" = "ID")
  )

cat("Number of matched subjects:", nrow(matched_data), "\n")

# 4. Read and vectorize NIfTI images
vectorize_nifti <- function(filepath) {
  img <- readnii(filepath)
  as.numeric(img)
}

non_vectorized_data <- function(filepath) {
  img <- readnii(filepath)
  img
}

predictor_list <- matched_data$filepath %>%
  map(~ vectorize_nifti(.x))

predictor_list_3d <- matched_data$filepath %>%
  map(~ non_vectorized_data(.x))

predictor_matrix <- do.call(rbind, predictor_list)
colnames(predictor_matrix) <- paste0("Vox_", seq_len(ncol(predictor_matrix)))

cat("Predictor matrix dimensions (Subjects x Voxels):", dim(predictor_matrix), "\n")
cat("Predictor matrix dimensions (3D) (Subjects x Voxels):", length(predictor_list_3d), "\n")

# 4.2 Re-resolute the imaging data based on 

load("./ABIDE_Data/AAL_partitions.RData")

# Get the unique elements of the array
unique_values <- sort(unique(c(AAL_256)))

# Create a mapping from original values to 0:256
remap <- setNames(0:(length(unique_values) - 1), unique_values)

# Remap the array using the mapping
AAL_256_Remap <- array(remap[as.character(AAL_256)], dim = dim(AAL_256))

AAL_256_neighbors <-  findSupervoxelNeighborsOptimized(AAL_256_Remap, dims = dim(AAL_256_Remap), unique_labels = unique(c(AAL_256_Remap)))
names(AAL_256_neighbors) <- unique(c(AAL_256_Remap))
predictors_super_voxels_aal256 <- lapply(predictor_list_3d,
                                         calculate_super_voxel_means,
                                         AAL_512 = AAL_256_Remap,
                                         unique_labels = unique(c(AAL_256_Remap)))

saveRDS(predictors_super_voxels_aal256, "./ABIDE_Data/predictors_super_voxels_aal256.rds")
saveRDS(AAL_256_neighbors, "./ABIDE_Data/AAL_256_neighbors.rds")

# Remove the first element of the list

# predictor_matrix_IFGtriangR <- predictor_matrix[, indices_IFGtriangR_flat]
# predictor_matrix_PoCGR <- predictor_matrix[, indices_PoCGR_flat]

# map 3d indices 
# indices_min_IFGtriangR <- apply(indices_IFGtriangR, 2, min)
# indices_max_IFGtriangR <- apply(indices_IFGtriangR, 2, max)
# indices_min_PoCGR <- apply(indices_PoCGR, 2, min)
# indices_max_PoCGR <- apply(indices_PoCGR, 2, max)

# # subset a cube bounding the ROI
# predictor_3d_IFGtriangR <- lapply(predictor_list_3d, function(x){
#   x[!indices_IFGtriangR] <- 0
#   x[indices_min_IFGtriangR[1]:indices_max_IFGtriangR[1],
#     indices_min_IFGtriangR[2]:indices_max_IFGtriangR[2],
#     indices_min_IFGtriangR[3]:indices_max_IFGtriangR[3]]
# }) 
# 
# predictor_3d_PoCGR <- lapply(predictor_list_3d, function(x){
#   x[!indices_PoCGR] <- 0
#   x[indices_min_PoCGR[1]:indices_max_PoCGR[1],
#     indices_min_PoCGR[2]:indices_max_PoCGR[2],
#     indices_min_PoCGR[3]:indices_max_PoCGR[3]]
# })

# 
# predictor_3d_IFGtriangR_cube = abind::abind(predictor_3d_IFGtriangR, along = 0)
# predictor_3d_PoCGR_cube = abind::abind(predictor_3d_PoCGR, along = 0)

# cat("Predictor matrix dimensions (3D) (Subjects x Voxels):", dim(predictor_3d_PoCGR_cube), "\n")
# cat("Predictor matrix dimensions (3D) (Subjects x Voxels):", dim(predictor_3d_IFGtriangR_cube), "\n")

predictors_super_voxels_aal256_mat <- abind::abind(predictors_super_voxels_aal256, along = 0)
if(!dir.exists("./processed_data")) dir.create("./processed_data", recursive = T)
saveRDS(predictors_super_voxels_aal256_mat, "./processed_data/predictors_super_voxels_aal256_mat.rds")

# 5. Extract outcome variable
outcome <- matched_data$DX_GROUP

# 6. Create covariate table
covariates <- matched_data %>%
  select(AGE_AT_SCAN, SEX, FIQ, HANDEDNESS_CATEGORY) %>%
  mutate(
    SEX = factor(SEX, levels = c(1, 2), labels = c("Male", "Female")),
    HANDEDNESS_CATEGORY = factor(HANDEDNESS_CATEGORY)
  )

# 7. Combine into final dataset
phenotype_data <- covariates %>%
  mutate(DX_GROUP = outcome)


final_dataset_super <- list(
  predictors = predictors_super_voxels_aal256_mat,
  outcome = outcome,
  covariates = covariates
)


# 8. Save the dataset
output_dir <- "./processed_data/"
dir.create(output_dir, showWarnings = FALSE, recursive = T)


saveRDS(final_dataset_super, file = file.path(output_dir, "final_dataset_super_256.rds"))
saveRDS(phenotype_data, file = file.path(output_dir, "phenotype_data_256.rds"))



# 9. Split the dataset into training and testing

# Set a seed for reproducibility
set.seed(20241111)

# Define the proportion of the data to use for training
train_index <- createDataPartition(final_dataset_super$outcome, p = 0.9, list = FALSE)

# Split data for IFGtriangR
train_dat <- list(
  predictors = final_dataset_super$predictors[train_index, ],
  outcome = final_dataset_super$outcome[train_index],
  covariates = final_dataset_super$covariates[train_index, ]
)

test_dat <- list(
  predictors = final_dataset_super$predictors[-train_index, ],
  outcome = final_dataset_super$outcome[-train_index],
  covariates = final_dataset_super$covariates[-train_index, ]
)





# Save training and testing datasets
saveRDS(train_dat, file = file.path(output_dir, "train_dat_256.rds"))
saveRDS(test_dat, file = file.path(output_dir, "test_dat_256.rds"))


cat("Training and testing datasets for have been saved.\n")



# ----------------------- DNN ----------------------- #
num_super_vox = length(unique(c(AAL_256_Remap)))
predictor_list_3d_super = predictor_list_3d
for (i in 1:length(predictor_list_3d)){
  if (i %% 100 == 0) cat(i, "\n")
  img = predictor_list_3d[[i]]
  for (j in 1:num_super_vox){
    img[AAL_256_Remap == j] = predictors_super_voxels_aal256_mat[i,j]
  }
  predictor_list_3d_super[[i]] = img
}


predictor_list_3d_super_mat = abind::abind(predictor_list_3d_super, along = 0)

# Split data for IFGtriangR
train_dat_nn <- list(
  predictors = predictor_list_3d_super_mat[train_index, , ,],
  outcome = final_dataset_super$outcome[train_index],
  covariates = final_dataset_super$covariates[train_index, ]
)

test_dat_nn <- list(
  predictors = predictor_list_3d_super_mat[-train_index,,, ],
  outcome = final_dataset_super$outcome[-train_index],
  covariates = final_dataset_super$covariates[-train_index, ]
)

saveRDS(predictor_list_3d_super, "./processed_data/predictor_list_3d_super.rds")
saveRDS(train_dat_nn, "./processed_data/train_dat_nn.rds")
saveRDS(test_dat_nn, "./processed_data/test_dat_nn.rds")
