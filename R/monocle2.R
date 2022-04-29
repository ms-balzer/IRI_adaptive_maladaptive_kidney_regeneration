library(dplyr)
library(Seurat)
library(data.table)
library(cowplot)
library(ggplot2)
library(monocle)
library(htmlwidgets)
library(OneR)
set.seed(123)

setwd('~/monocle2/')

#===================================================================
#===================================================================
#========= LOADING DATA, CREATING CDS & ADDING METADATA ============
#===================================================================
#===================================================================

### STEP 1: Pick random sample of cells with equal nCells from a number of interesting clusters 
#read in clustered information
seurat <- readRDS("~/Seurat/PT1.7kbarcodes.rds") #this object includes Phase and IRI score classification
Idents(seurat) <- seurat$LIGER_clusters
table(Idents(seurat))

### STEP 3: create CDS object
## Building the necessary parts for a basic cds
# part one, gene annotations
gene_annotation <- as.data.frame(seurat@assays[["RNA"]]@counts@Dimnames[[1]], 
                                 row.names = seurat@assays[["RNA"]]@counts@Dimnames[[1]])
gene_annotation <- AnnotatedDataFrame(gene_annotation)
colnames(gene_annotation) <- "gene_short_name"

# part two, cell information
cell_metadata <- as.data.frame(seurat@assays[["RNA"]]@counts@Dimnames[[2]], 
                               row.names = seurat@assays[["RNA"]]@counts@Dimnames[[2]])
cell_metadata <- AnnotatedDataFrame(cell_metadata)
colnames(cell_metadata) <- "barcode"

# part three, counts sparse matrix (from fresh sample)
New_matrix <- seurat@assays[["RNA"]]@counts
expression_matrix <- New_matrix

## Construct the basic cds object
cds_from_seurat <- newCellDataSet(as(expression_matrix, "sparseMatrix"),
                                  phenoData = cell_metadata,
                                  featureData = gene_annotation,
                                  expressionFamily=negbinomial.size())
cds <- cds_from_seurat

### STEP 4: ADD METADATA
cds@phenoData@data$orig.ident <- seurat$orig.ident
cds@phenoData@data$time.pod <- seurat$time.pod
cds@phenoData@data$exp.cond <- seurat$exp.cond
cds@phenoData@data$exp.time <- seurat$exp.time
cds@phenoData@data$LIGER_clusters <- seurat$LIGER_clusters
cds@phenoData@data$Phase <- seurat$Phase
cds@phenoData@data$S.Score <- seurat$S.Score
cds@phenoData@data$G2M.Score <- seurat$G2M.Score
cds@phenoData@data$IRI_label <- seurat$IRI_label
cds@phenoData@data$IRI.Score <- seurat$IRI.Score



#=================================================
#=================================================
#================= PREPROCESSING =================
#=================================================
#=================================================
#Estimate size factors and dispersions
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

#Filtering low-quality cells
cds <- detectGenes(cds, min_expr = 0.1)
cds #still 1650 cells
print(head(fData(cds)))



#==================================================================
#==================================================================
#================= Classifying and Counting Cells =================
#==================================================================
#==================================================================

disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
cds <- reduceDimension(cds, max_components = 2, num_dim = 15, #need to adjust num_dim depending on elbow plot
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 11)



#=========================================================================
#=========================================================================
#================= Constructing Single Cell Trajectories =================
#=========================================================================
#=========================================================================
#========= Trajectory step 1: choose genes that define a cell's progress
#we use the fully unsupervised approach and do not do feature selection here:



#========= Trajectory step 2: reduce data dimensionality
cds <- reduceDimension(cds, 
                       max_components = 2,
                       reduction_method = 'DDRTree', 
                       verbose = T)



#========= Trajectory step 3: order cells along the trajectory
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = NULL, num_paths = NULL, reverse = T) #Need to reverse trajectory in this case



#====================================================================
#====================================================================
#================= Differential Expression Analysis =================
#====================================================================
#====================================================================

#================= Step 1: Basic Differential Analysis =================
marker_genes <- row.names(fData(cds))



#================= Step 2: Finding Genes that Distinguish Cell Type or State =================
#DEG across LIGER_clusters
diff_test_res_LIGER_clusters <- differentialGeneTest(cds,
                                                     fullModelFormulaStr = "~LIGER_clusters")
diff_test_res_LIGER_clusters[,c("gene_short_name", "pval", "qval")]



#================= Step 3: Finding Genes that Change as a Function of Pseudotime =================
to_be_tested <- row.names(subset(fData(cds),
                                 gene_short_name %in% c("Il1b")))
cds_subset <- cds[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res[,c("gene_short_name", "pval", "qval")]



#================= Step 4: Clustering Genes by Pseudotemporal Expression Pattern =================
marker_genes <- row.names(fData(cds))
diff_test_res_Pseudotime <- differentialGeneTest(cds[marker_genes,],
                                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names_Pseudotime <- row.names(subset(diff_test_res_Pseudotime, qval < 0.1))



#=========================================================================
#=========================================================================
#================= Analyzing Branches in Single-Cell Trajectories ========
#=========================================================================
#=========================================================================
BEAM_res <- BEAM(cds, branch_point = 1, cores = 8)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
















