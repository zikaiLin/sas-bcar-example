plot_non_linear_GP = function(mod,
                              x,
                              dims_to_plot = 1,
                              poly_degree = 20,
                              a = 1,
                              b = 1000) {
  library(BayesGPfit)
  library(ggplot2)
  x = as.matrix(x[,dims_to_plot])
  grids = GP.generate.grids(d = 1,
                            num_grids = 50,
                            grids_lim = c(0 , 1))
  basis = GP.eigen.funcs.fast(x, poly_degree = poly_degree, a, b)
 

  par(mar = c(5, 5, 5, 5))
  x_order = order(x)
  to_plot = data.frame(x = x[x_order], y = c(basis %*% mod$post_mean$theta[, , dims_to_plot])[x_order])
  plt = ggplot(data = to_plot, mapping = aes(x = x, y = y)) + 
    geom_point(size = 6, colour = "red") +  # Scatter plot
    # geom_smooth(method = "loess", colour = "blue", se = FALSE) +  # Non-linear trend line without confidence interval
    xlab("Signal Intensity") +
    ylab("Function Values") +  # Mathematical expression for y-axis label
    xlim(range(to_plot$x, na.rm = TRUE)) +  # Adjust for NA values if present
    ylim(range(to_plot$y, na.rm = TRUE)) +  # Adjust for NA values if present
    theme_bw(base_size = 24) +  # Use base_size for base theme text size
    theme(
      text = element_text(size = 36),  # Enlarge text size for all text elements
      axis.title = element_text(size = 28),  # Enlarge axis labels
      axis.text = element_text(size = 28)  # Enlarge axis tick text
    )
  plt

}


K179_mcmc <- readRDS("./BCI_results/BCI_result_best/K179_mcmc.rds")
K179 <- readRDS("./BCI/data/K179.rds")
which(K179_mcmc$post_mean$delta == 1)

for(i in c(2,14,30,37,63,104,138,158,163,185,203,269,295,329,358,380)){
  plt = plot_non_linear_GP(
    mod = K179_mcmc,
    x = K179$X_test/max(abs(K179$X_test)),
    dims_to_plot = i,
    poly_degree = K179_mcmc$best_kernel$poly,
    a = 0.01,
    b = K179_mcmc$best_kernel$b
  )
  ggsave(plt,device = "png", filename = sprintf("./BCI_results/plot_K179_non_linear//dim_%d.png",i), width = 11.3, height = 5)
}

