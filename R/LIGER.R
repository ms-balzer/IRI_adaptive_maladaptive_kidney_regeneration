library(liger)
library(Seurat)
set.seed(123)

setwd('~/LIGER/')

# Stage I: Preprocessing and Normalization
matrix_list <- read10X(sample.dirs =c("/~/data/cellranger_count/norm1/outs/filtered_feature_bc_matrix",
                                      "/~/data/cellranger_count/norm2/outs/filtered_feature_bc_matrix",
                                      "/~/data/cellranger_count/norm3/outs/filtered_feature_bc_matrix",
                                      "/~/data/cellranger_count/norm4/outs/filtered_feature_bc_matrix",
                                      "/~/data/cellranger_count/norm5/outs/filtered_feature_bc_matrix",
                                      "/~/data/cellranger_count/Scl/outs/filtered_feature_bc_matrix",
                                      "/~/data/cellranger_count/IRF1/outs/filtered_feature_bc_matrix",
                                      "/~/data/cellranger_count/IRF2/outs/filtered_feature_bc_matrix",
                                      "/~/data/cellranger_count/IRF3/outs/filtered_feature_bc_matrix",
                                      "/~/data/cellranger_count/IRF4/outs/filtered_feature_bc_matrix",
                                      "/~/data/cellranger_count/IRF5/outs/filtered_feature_bc_matrix",
                                      "/~/data/cellranger_count/IRF6/outs/filtered_feature_bc_matrix",
                                      "/~/data/cellranger_count/IRF7/outs/filtered_feature_bc_matrix",
                                      "/~/data/cellranger_count/IRF8/outs/filtered_feature_bc_matrix",
                                      "/~/data/cellranger_count/IRF9/outs/filtered_feature_bc_matrix",
                                      "/~/data/cellranger_count/IRF10/outs/filtered_feature_bc_matrix",
                                      "/~/data/cellranger_count/IRF11/outs/filtered_feature_bc_matrix",
                                      "/~/data/cellranger_count/IRF12/outs/filtered_feature_bc_matrix"), 
                       sample.names = c("norm1", "norm2", "norm3", "norm4", "norm5", "Scl", 
                                        "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", 
                                        "IRF7", "IRF8", "IRF9", "IRF10", "IRF11", "IRF12"), 
                       merge = F)

colnames(matrix_list$norm1$`Gene Expression`) <- paste0("norm1_",colnames(matrix_list$norm1$`Gene Expression`))
colnames(matrix_list$norm2$`Gene Expression`) <- paste0("norm2_",colnames(matrix_list$norm2$`Gene Expression`))
colnames(matrix_list$norm3$`Gene Expression`) <- paste0("norm3_",colnames(matrix_list$norm3$`Gene Expression`))
colnames(matrix_list$norm4$`Gene Expression`) <- paste0("norm4_",colnames(matrix_list$norm4$`Gene Expression`))
colnames(matrix_list$norm5$`Gene Expression`) <- paste0("norm5_",colnames(matrix_list$norm5$`Gene Expression`))
colnames(matrix_list$Scl$`Gene Expression`) <- paste0("Scl_",colnames(matrix_list$Scl$`Gene Expression`))
colnames(matrix_list$IRF1$`Gene Expression`) <- paste0("IRF1_",colnames(matrix_list$IRF1$`Gene Expression`))
colnames(matrix_list$IRF2$`Gene Expression`) <- paste0("IRF2_",colnames(matrix_list$IRF2$`Gene Expression`))
colnames(matrix_list$IRF3$`Gene Expression`) <- paste0("IRF3_",colnames(matrix_list$IRF3$`Gene Expression`))
colnames(matrix_list$IRF4$`Gene Expression`) <- paste0("IRF4_",colnames(matrix_list$IRF4$`Gene Expression`))
colnames(matrix_list$IRF5$`Gene Expression`) <- paste0("IRF5_",colnames(matrix_list$IRF5$`Gene Expression`))
colnames(matrix_list$IRF6$`Gene Expression`) <- paste0("IRF6_",colnames(matrix_list$IRF6$`Gene Expression`))
colnames(matrix_list$IRF7$`Gene Expression`) <- paste0("IRF7_",colnames(matrix_list$IRF7$`Gene Expression`))
colnames(matrix_list$IRF8$`Gene Expression`) <- paste0("IRF8_",colnames(matrix_list$IRF8$`Gene Expression`))
colnames(matrix_list$IRF9$`Gene Expression`) <- paste0("IRF9_",colnames(matrix_list$IRF9$`Gene Expression`))
colnames(matrix_list$IRF10$`Gene Expression`) <- paste0("IRF10_",colnames(matrix_list$IRF10$`Gene Expression`))
colnames(matrix_list$IRF11$`Gene Expression`) <- paste0("IRF11_",colnames(matrix_list$IRF11$`Gene Expression`))
colnames(matrix_list$IRF12$`Gene Expression`) <- paste0("IRF12_",colnames(matrix_list$IRF12$`Gene Expression`))

Co_mm_IRF <- createLiger(list(norm1 = matrix_list$norm1$`Gene Expression`,
                              norm2 = matrix_list$norm2$`Gene Expression`,
                              norm3 = matrix_list$norm3$`Gene Expression`,
                              norm4 = matrix_list$norm4$`Gene Expression`,
                              norm5 = matrix_list$norm5$`Gene Expression`,
                              Scl = matrix_list$Scl$`Gene Expression`,
                              IRF1 = matrix_list$IRF1$`Gene Expression`,
                              IRF2 = matrix_list$IRF2$`Gene Expression`,
                              IRF3 = matrix_list$IRF3$`Gene Expression`,
                              IRF4 = matrix_list$IRF4$`Gene Expression`,
                              IRF5 = matrix_list$IRF5$`Gene Expression`,
                              IRF6 = matrix_list$IRF6$`Gene Expression`,
                              IRF7 = matrix_list$IRF7$`Gene Expression`,
                              IRF8 = matrix_list$IRF8$`Gene Expression`,
                              IRF9 = matrix_list$IRF9$`Gene Expression`,
                              IRF10 = matrix_list$IRF10$`Gene Expression`,
                              IRF11 = matrix_list$IRF11$`Gene Expression`,
                              IRF12 = matrix_list$IRF12$`Gene Expression`))
Co_mm_IRF #162607 cells

#Filter high quality cells = subsetLiger with only those cells that are high quality in Seurat
seur_Co_mm_IRF <- readRDS("/~/data/seurat/Co_mm_IRF/Co_mm_IRF_1.rds") #determined by Seurat standard pipeline
seur_Co_mm_IRF #113579 cells
cells.to.subset <- colnames(seur_Co_mm_IRF)
Co_mm_IRF_mtfilter <- subsetLiger(Co_mm_IRF, cells.use = cells.to.subset, remove.missing = T)
#subset 113579 cells of 162607 total cells.

Co_mm_IRF_mtfilter <- normalize(Co_mm_IRF_mtfilter)
Co_mm_IRF_mtfilter <- selectGenes(Co_mm_IRF_mtfilter)
print(length(Co_mm_IRF_mtfilter@var.genes)) # Print out the number of variable genes for LIGER analysis: 8183
Co_mm_IRF_mtfilter <- scaleNotCenter(Co_mm_IRF_mtfilter)



# Stage II: Joint Matrix Factorization (CPU heavy)
suggested.k <- suggestK(Co_mm_IRF_mtfilter)
suggested.l <- suggestLambda(Co_mm_IRF_mtfilter, k = suggested.k)
Co_mm_IRF_mtfilter <- optimizeALS(Co_mm_IRF_mtfilter,
                                        k = suggested.k,
                                        lambda = suggested.l,
                                        thresh = 1e-06,
                                        max.iters = 30,
                                        nrep = 1,
                                        rand.seed = 1)
Co_mm_IRF_mtfilter <- runTSNE(Co_mm_IRF_mtfilter, use.raw = T)



# Stage III: Quantile Normalization and Joint Clustering
Co_mm_IRF_mtfilter <- quantile_norm(Co_mm_IRF_mtfilter)
Co_mm_IRF_mtfilter <- louvainCluster(Co_mm_IRF_mtfilter, resolution = 0.25)



# Stage IV: Visualization and Downstream Analysis
Co_mm_IRF_mtfilter <- runTSNE(Co_mm_IRF_mtfilter)
Co_mm_IRF_mtfilter <- runUMAP(Co_mm_IRF_mtfilter, distance = 'cosine', n_neighbors = 30, min_dist = 0.3)



# Stage V: Liger to Seurat export
object <- ligerToSeurat(
  Co_mm_IRF_mtfilter,
  nms = NULL,
  renormalize = T,
  use.liger.genes = T,
  by.dataset = F)

saveRDS(object, file="/~/data/LIGER/object_res0.25.rds")


#================= EXIT =================
q(save = "no", status = 0, runLast = TRUE)
