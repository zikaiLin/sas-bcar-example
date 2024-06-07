# Create partition of an image

create_partition <- function(dim_x, dim_y, n_part = 4){
  x_nseq = dim_x %/% n_part
  y_nseq = dim_y %/% n_part
  
  x_seq = c(floor(seq(1,dim_x-x_nseq, length.out = n_part)), dim_x)
  y_seq = c(floor(seq(1,dim_y-y_nseq, length.out = n_part)), dim_y)
  
  partition_matrix = matrix(NA, dim_x, dim_y)
  count = 1
  for(i in 2:length(x_seq)-1){
    x_interval = c(x_seq[i], x_seq[i+1])
    for(j in 2:length(y_seq)-1){
      y_interval = c(y_seq[j], y_seq[j+1])
      partition_matrix[x_interval[1]: x_interval[2], y_interval[1]:y_interval[2]] =  count
      count = count +1
    }
  }
  return(partition_matrix)
}