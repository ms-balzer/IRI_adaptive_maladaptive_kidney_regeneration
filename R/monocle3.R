rm(list=ls())
library(dplyr)
library(Seurat)
library(data.table)
library(cowplot)
library(ggplot2)
library(monocle3)
library(htmlwidgets)
set.seed(123)

setwd('~/monocle3/')


#=========================================================================
#=========================================================================
#================= Convert Seurat object to monocle CDS ==================
#=========================================================================
#=========================================================================

### STEP 1: Pick random sample of cells with equal nCells from a number of interesting clusters 
#read in clustered information
seurat <- readRDS("~/Seurat/object_res0.25_freshPTsubset_final_res0.25_inclIRIscore.rds") #this object includes Phase and IRI score classification
seurat #19614 features across 28385 samples within 1 assay 


#now subset to only IRI_short_1d, IRI_short_3d, IRI_short_14d, Control
Idents(seurat) <- seurat$exp.time
table(Idents(seurat))
seurat@reductions$umap@assay.used <- "RNA"
PTlongControl <- subset(seurat, idents = c("IRI_long_1d", "IRI_long_3d", "IRI_long_14d", "Control"))
PTlongControl #19614 features across 23639 samples within 1 assay 

#evenly subsample those timepoints
table(Idents(PTlongControl))
cells.to.sample <- min(table(Idents(PTlongControl)))
cells.to.sample #497

# Sample from other clusters as many cells as there are cells in cluster12
# For reproducibility, set a random seed
set.seed(cells.to.sample)
sampled.cells_0 <- sample(x = WhichCells(seurat, idents = "IRI_long_1d"), size = cells.to.sample, replace = F)
sampled.cells_1 <- sample(x = WhichCells(seurat, idents = "IRI_long_3d"), size = cells.to.sample, replace = F)
sampled.cells_2 <- sample(x = WhichCells(seurat, idents = "IRI_long_14d"), size = cells.to.sample, replace = F)
sampled.cells_3 <- sample(x = WhichCells(seurat, idents = "Control"), size = cells.to.sample, replace = F)

# Create a vector of cells on which variable genes will be computed
cells.to.subset <- c(sampled.cells_0,
                     sampled.cells_1,
                     sampled.cells_2, 
                     sampled.cells_3)
#subset to PTtrajectoryM
PTlongControl@reductions$umap@assay.used <- "RNA"
PTtrajectoryM <- subset(PTlongControl, cells = cells.to.subset)
PTtrajectoryM #19614 features across 1988 samples within 1 assay 
seurat <- PTtrajectoryM
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat)) #seurat object does not have PCA reduction, need to run before converting to CDS


### Building the necessary parts for a basic cds
# part one, gene annotations
gene_annotation <- as.data.frame(rownames(seurat@reductions[["pca"]]@feature.loadings), 
                                 row.names = rownames(seurat@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

# part two, cell information
cell_metadata <- as.data.frame(seurat@assays[["RNA"]]@counts@Dimnames[[2]], 
                               row.names = seurat@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

# part three, counts sparse matrix
New_matrix <- seurat@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(seurat@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

### Construct the basic cds object
cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)

### Construct and assign the made up partition
recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

### Assign the cluster info
list_cluster <- seurat@meta.data$LIGER_clusters
names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]
cds_from_seurat@clusters@listData[["UMAP"]][["LIGER_clusters"]] <- list_cluster

### Assign clusters
cds_from_seurat@colData@listData$orig.ident <- seurat$orig.ident
cds_from_seurat@colData@listData$time.pod <- seurat$time.pod
cds_from_seurat@colData@listData$exp.cond <- seurat$exp.cond
cds_from_seurat@colData@listData$exp.time <- seurat$exp.time
cds_from_seurat@colData@listData$LIGER_clusters <- seurat$LIGER_clusters
cds_from_seurat@colData@listData$Phase <- seurat$Phase
cds_from_seurat@colData@listData$S.Score <- seurat$S.Score
cds_from_seurat@colData@listData$G2M.Score <- seurat$G2M.Score
cds_from_seurat@colData@listData$IRI_label <- seurat$IRI_label
cds_from_seurat@colData@listData$IRI.Score <- seurat$IRI.Score



cds <- cds_from_seurat
rm(cds_from_seurat)



#==================================================
#==================================================
#================= PREPROCESSING ==================
#==================================================
#==================================================

#================= Step 1: Normalize & pre-process the data =================
# CPU intensive, time-consuming.
cds <- preprocess_cds(cds, num_dim = 100, cores=8, method = "PCA")
cds <- preprocess_cds(cds, num_dim = 100, cores=8, method = "LSI")



#================= Step 2: Reduce dimensionality and visualize the cells (UMAP by default) =================
cds <- reduce_dimension(cds, reduction_method = "UMAP", umap.fast_sgd = FALSE, preprocess_method = "PCA", cores=1)
cds <- reduce_dimension(cds, reduction_method = "tSNE", preprocess_method = "PCA", cores=1)
cds <- reduce_dimension(cds, reduction_method = "PCA", preprocess_method = "PCA", cores=1)
cds <- reduce_dimension(cds, reduction_method = "LSI", preprocess_method = "LSI", cores=1)



#==================================================================================================================
#==================================================================================================================
#================= MANUALLY COPY PCA COORDINATES TO UMAP SLOT TO PROCEED WITH MONOCLE3 FUNCTIONS ==================
#==================================================================================================================
#==================================================================================================================
### Assign the cluster info
list_cluster <- seurat@meta.data$LIGER_clusters
names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]
cds@clusters@listData[["UMAP"]][["LIGER_clusters"]] <- list_cluster

### Assign clusters
cds@colData@listData$orig.ident <- seurat$orig.ident


### COPY PCA OVER UMAP SLOT
cds@int_colData@listData$reducedDims@listData$UMAP <- as.matrix(cds@int_colData@listData$reducedDims@listData$PCA[,1:2]) #(obviously only use the first two PCs as embeddings for fake UMAP!)
colnames(cds@int_colData@listData$reducedDims@listData$UMAP)



#================= Cluster the cells =================
res=3e-3
cds = cluster_cells(cds, resolution=res, reduction_method = "UMAP", k=29) #need to run with UMAP before learning graph!

#by clusters
pdf(file = "/home/balzermi/data/monocle3/LtoS_Co_mm_IRF_noKsp_mtfilter_freshPTsubset_PTtrajectoryM_pcatoumap/plot_cells_by_clusters_PCAtoUMAP.pdf")
plot_cells(cds, color_cells_by="cluster", reduction_method = "UMAP")
dev.off()

#by partitions
pdf(file = "/home/balzermi/data/monocle3/LtoS_Co_mm_IRF_noKsp_mtfilter_freshPTsubset_PTtrajectoryM_pcatoumap/plot_cells_by_partition_PCAtoUMAP.pdf")
plot_cells(cds, color_cells_by="partition", group_cells_by="partition", reduction_method = "UMAP")
dev.off()



#================= Find marker genes expressed by each cluster ================= SKIPPED THIS
#by cluster
marker_test_res <- top_markers(cds, group_cells_by="cluster", 
                               reference_cells=1000, cores=8)

# The data frame marker_test_res contains a number of metrics for how specifically expressed 
# each gene is in each cluster We could group the cells according to cluster, partition, 
# or any categorical variable in colData(cds). You can rank the table according to one or more 
# of the specificity metrics and take the top gene for each cluster. For example, pseudo_R2 is
# one such measure. We can rank markers according to pseudo_R2 like this:
top50_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(50, pseudo_R2)

top1_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top1_specific_marker_ids <- unique(top1_specific_markers %>% pull(gene_id))


# It's often informative to look at more than one marker, which you can do just by changing the 
# first argument to top_n():
top3_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top3_specific_marker_ids <- unique(top3_specific_markers %>% pull(gene_id))



#=========================================================================
#=========================================================================
#================= Constructing single-cell trajectories =================
#=========================================================================
#=========================================================================

#================= Learn the trajectory graph =================
# Next, we will fit a principal graph within each partition using the learn_graph() function:
cds <- learn_graph(cds, verbose = T) #this worked because standard k=25 and above, k was set to 29 at the cluster_cells stage!

#================= Order the cells in pseudotime =================
get_earliest_principal_node <- function(cds, assigned_cell_type="IRI_long_1d"){
  cell_ids <- which(colData(cds)[, "assigned_cell_type"] == assigned_cell_type)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))



#================= Finding genes that change as a function of pseudotime =================
NP_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=8)




#================= EXIT =================
q(save = "no", status = 0, runLast = TRUE)




