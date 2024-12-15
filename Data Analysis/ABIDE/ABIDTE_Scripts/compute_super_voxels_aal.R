calculate_super_voxel_means <- function(A, AAL_512, unique_labels) {
  # Ensure the dimensions of A and AAL_512 match
  if (!all(dim(A) == dim(AAL_512))) {
    stop("Dimensions of A and AAL_512 must match.")
  }
  
  # Flatten the 3D arrays into vectors
  A_flat <- as.vector(A)
  AAL_512_flat <- as.vector(AAL_512)
  
  # Get unique super voxel labels
  # unique_labels <- unique(AAL_512_flat)
  
  # Initialize a vector to store the mean values for each super voxel label
  mean_values <- numeric(length(unique_labels))
  
  # Calculate the mean value for each super voxel label
  for (i in seq_along(unique_labels)) {
    label <- unique_labels[i]
    # Find the indices of voxels with the current label
    voxel_indices <- which(AAL_512_flat == label)
    # Compute the mean of the corresponding values in A
    mean_values[i] <- mean(A_flat[voxel_indices], na.rm = TRUE)
  }
  
  # Return a named vector with the unique labels and their mean values
  names(mean_values) <- unique_labels
  return(mean_values)
}
