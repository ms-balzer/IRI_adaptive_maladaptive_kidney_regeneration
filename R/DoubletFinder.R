library(dplyr)
library(Seurat)
library(data.table)
library(DoubletFinder)
set.seed(123)


object <- readRDS(file="/~/LIGER/object_res0.25.rds")

#========= Remove Doublets ===========
## Pre-process Seurat object (standard) ---------------------------------------------------------------------------------------
seu_kidney <- NormalizeData(object)
seu_kidney <- FindVariableFeatures(seu_kidney, selection.method = "vst", nfeatures = 2000)
seu_kidney <- ScaleData(seu_kidney)
seu_kidney <- RunPCA(seu_kidney)
seu_kidney <- RunUMAP(seu_kidney, dims = 1:10)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:10, sct = FALSE, num.cores=8)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
pdf(file = "/~/plots_DoubletFinder.pdf", width=5, height=5)
bcmvn_kidney <- find.pK(sweep.stats_kidney)
dev.off()
bcmvn_kidney$pK[which.max(bcmvn_kidney$BCmetric)] #which pK returns the maximum BCmetric?

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seu_kidney@meta.data$LIGER_clusters3
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.091*length(colnames(seu_kidney)))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_10336", sct = FALSE)
saveRDS(seu_kidney, "/~/seu_kidney_DoubletFinder.rds")

