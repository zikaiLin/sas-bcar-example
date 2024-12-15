library(RNifti)
aal116 <- RNifti::readNifti("./ABIDE/aal116_3mm.nii.gz")
aal116_table <- read.table("./ABIDE/AAL90/Node_AAL116.node")
Rcpp::sourceCpp("./Scripts/Zikai/update_pj_mat_3d.cpp")

# 14   50.33  30.16  14.17  3  3 IFGtriang.R
# 58   41.43 -25.49  52.55  1  3      PoCG.R

sum(aal116 == 14)

indices_IFGtriangR <- which(aal116 == 14, arr.ind = TRUE)
indices_PoCGR <- which(aal116 == 58, arr.ind = TRUE)

indices_IFGtriangR_flat <- which(aal116 == 14)
indices_PoCGR_flat <- which(aal116 == 58)


neighbor_list_IFGtriangR <- build_neighbor_list(indices_IFGtriangR)
neighbor_list_PoCGR <- build_neighbor_list(indices_PoCGR)



saveRDS(neighbor_list_IFGtriangR, "./ABIDE/neighbor_list_IFGtriangR.rds")
saveRDS(neighbor_list_PoCGR, "./ABIDE/neighbor_list_PoCGR.rds")
saveRDS(indices_IFGtriangR, "./ABIDE/indices_IFGtriangR.rds")
saveRDS(indices_PoCGR, "./ABIDE/indices_PoCGR.rds")
saveRDS(indices_IFGtriangR_flat, "./ABIDE/indices_IFGtriangR_flat.rds")
saveRDS(indices_PoCGR_flat, "./ABIDE/indices_PoCGR_flat.rds")


delta_vec = rbinom(length(neighbor_list_PoCGR), 1, 0.5)


test_avg_inclusion = compute_average_inclusion(delta_vec, neighbor_list_PoCGR)
test_avg_inclusion[1]
print(delta_vec[neighbor_list_PoCGR[[1]]])
