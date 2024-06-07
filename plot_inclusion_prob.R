plot_mnist <- function(img,layout = NULL){
  par(mar=c(0,0,0,0))
  if(!is.matrix(img)){
    
    image(1:28, 1:28, matrix(img, nrow=28)[ , 28:1], 
          col = gray(seq(0, 1, 0.05)), xlab = "", ylab="",axes = FALSE,asp=1)
  } else{
    n = nrow(img)
    m = round(sqrt(n))
    k = ceiling(n/m)
    if(is.null(layout)){
      layout = c(m,k)
    }
    par(mfrow=layout)
    for(i in 1:n){
      image(1:28, 1:28, matrix(img[i,], nrow=28)[ , 28:1], 
            col = gray(seq(0, 1, 0.05)), xlab = "", ylab="",axes = FALSE,asp=1)
    }
  }
  par(mfrow=c(1,1),mar=c(5, 4, 4, 2) + 0.1)
}


