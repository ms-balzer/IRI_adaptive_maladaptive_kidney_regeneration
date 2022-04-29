
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
set.seed(123)

setwd("~/VelocytoR/")

#==================================================================
#==============Estimating RNA Velocity using Seurat================
#==================================================================
### Load filtered Seurat object
seurat_object <- readRDS(file = "/Seurat/object_res0.25_freshPTsubset_final_res0.25.rds")
Idents(seurat_object) <- seurat_object$GMMclusters
table(Idents(seurat_object))
seurat_object@reductions$umap@assay.used <- "RNA"
PTtrajectory <- subset(seurat_object, idents = "notsubset", invert=T) #subset all GMM clusters except "notsubset"
PTtrajectory$GMMclusters <- factor(PTtrajectory$GMMclusters, levels = c("1", "2", "3", "4", "5"))
table(PTtrajectory$GMMclusters)
PTtrajectory@reductions$gmm <- PTtrajectory@reductions$umap #copy umap to new GMM slot
dr_coords <- readRDS('~/SCENIC/int/dr_coords.Rds') #load GMM
PTtrajectory@reductions$gmm@cell.embeddings <- dr_coords #overwrite new slot with GMM cell.embeddings
rm(seurat_object)

### Load loom file
ldat <- ReadVelocity(file = "~/loompy/Co_mm_IRI_noKsp.loom")
for (i in names(x = ldat)) {
  ### Store assay in a new variable
  assay <- ldat[[i]]
  
  ### Rename cell names in loom file to match cell names in Seurat object
  colnames(assay) <- gsub(':', '_', colnames(assay))
  colnames(assay) <- gsub('comb_', '', colnames(assay))
  colnames(assay) <- gsub('x', '', colnames(assay))
  
  ### Subset to filtered cells in Seurat object
  assay <- assay[,colnames(PTtrajectory)]
  
  ### Add assay to Seurat object
  PTtrajectory[[i]] <- CreateAssayObject(counts = assay)
}

### Proceed with Velocity in R
Idents(PTtrajectory) <- PTtrajectory$GMMclusters
PTtrajectory@reductions$pca <- PTtrajectory@reductions$gmm
rd2 <- readRDS("~/Slingshot/rd2.rds") #load rd2
PTtrajectory@reductions$pca@cell.embeddings <- rd2 #overwrite new slot with GMM cell.embeddings
#this is exactly the same as the gmm! PROCEED!

#RunVelocity
PTtrajectory5 <- RunVelocity(object = PTtrajectory, deltaT = 1, kCells = 25, fit.quantile = 0.02)

