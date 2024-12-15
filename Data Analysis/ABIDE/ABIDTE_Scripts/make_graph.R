# Build a graph

chanel2ind = matrix(1:(26*16), nrow = 26, ncol = 16)
dist_mat = as.matrix(read.csv("/Volumes/easystore/Dropbox/U-Mich/Research/Scalar-on-image regression via GP Additive Model/BCI/distance_mat.csv")[,-1])

library(igraph)
dist_thres = 2
time_thres = 1
adjacency_mat = matrix(NA, 416,416)

for(node in 1:415){
    chanel_i = which(chanel2ind == node, arr.ind = T)[2]
    tp_i = which(chanel2ind == node, arr.ind = T)[1]
    
    for (j in (node+1):416){
      chanel_j = which(chanel2ind == j, arr.ind = T)[2]
      tp_j = which(chanel2ind == j, arr.ind = T)[1]
      if (abs(tp_j - tp_i) <= time_thres & dist_mat[chanel_i, chanel_j] <= dist_thres){
        adjacency_mat[node, j] = adjacency_mat[j, node] = 1
      }
    }
    
}

graph_BCI = igraph::graph.adjacency(adjacency_mat, mode = "undirected", weighted = NULL)

plot(graph_BCI, layout = layout.circle(graph_BCI))
E(graph = graph_BCI)

saveRDS(graph_BCI, file = "./BCI_graph.rds")
