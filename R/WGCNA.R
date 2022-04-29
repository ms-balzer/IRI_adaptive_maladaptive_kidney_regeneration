library(dplyr)
library(Seurat)
library(data.table)
library(cowplot)
library(ggplot2)
library(WGCNA)
enableWGCNAThreads()
library(flashClust)
set.seed(123)

setwd('~/WGCNA/')

object <- readRDS("~/Seurat/object_final_res0.25_inclIRIscore.rds")
object@reductions$umap@assay.used <- "RNA"
object$allcells <- "allcells"
table(object$allcells, object$exp.time)



################################################################################
# Metacell aggregation
################################################################################
library(monocle3)
library(cicero)
seurat_list <- list()
k = 25
celltypes <- unique(object$allcells)
celltypes <- celltypes[celltypes != 'PER.END']
for (cur_celltype in celltypes){
  condition_list <- list()
  for(condition in unique(object$exp.time)){
    print(paste(cur_celltype, condition))
    cur_seurat <- subset(object, allcells == cur_celltype & exp.time == condition)
    expr_matrix <- GetAssayData(cur_seurat, slot='data')
    genes <- data.frame(as.character(rownames(expr_matrix)))
    rownames(genes) <- rownames(expr_matrix)
    genes <- as.data.frame(cbind(genes,genes))
    colnames(genes) <- c("GeneSymbol", "gene_short_name")
    cds <- new_cell_data_set(
      expr_matrix,
      cell_metadata=cur_seurat@meta.data,
      gene_metadata=genes)
    reducedDims(cds)$UMAP <- cur_seurat@reductions$umap@cell.embeddings
    umap_coords <- reducedDims(cds)$UMAP
    metacell_cds <- make_cicero_cds(cds, reduced_coordinates=umap_coords, k=k, size_factor_normalize=FALSE)
    metacell_seurat <- CreateSeuratObject(counts = exprs(metacell_cds)/k)
    metacell_seurat$allcells <- cur_celltype
    metacell_seurat$exp.time <- condition
    metacell_seurat <- RenameCells(metacell_seurat, new.names=paste0(cur_celltype, '_', condition, '_', seq(1:ncol(metacell_seurat))))
    condition_list[[condition]] <- metacell_seurat
  }
  seurat_list[[cur_celltype]] <- merge(condition_list[[1]], y=condition_list[2:length(condition_list)])
}

metacell_seurat <- seurat_list[[1]]
metacell_seurat <- FindVariableFeatures(metacell_seurat, nfeatures=3000)
metacell_seurat <- ScaleData(metacell_seurat, features = VariableFeatures(metacell_seurat))
metacell_seurat <- RunPCA(metacell_seurat, features=VariableFeatures(metacell_seurat))
metacell_seurat <- RunUMAP(metacell_seurat, reduction='pca', dims=1:25, min.dist = 0.9, n.neighbors = 100L, spread = 5)
metacell_seurat <- metacell_seurat %>% 
  FindNeighbors(reduction = "pca", dims = 1:15) %>% 
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0), save.SNN=TRUE) %>% 
  identity()
metacell_seurat@meta.data$exp.time <- factor(metacell_seurat@meta.data$exp.time, 
                                             levels = c("Control", "IRI_short_1d", "IRI_short_3d", "IRI_short_14d", "IRI_long_1d", "IRI_long_3d", "IRI_long_14d"))



####Code by Vivek Swarup, PhD UC Irvine #############
###contact vswarup@uci.edu################
#####################################################
##Single-nuclei Consensus Weighted Gene Co-expression Network Analysis (scWGCNA)
##ROSMAP data was downloaded and processed as mentioned in Morabito et al., https://doi.org/10.1101/695221
## Mathys et al data was downloaded from Synapse (syn18485175; doi:10.7303/syn18485175)
## Our snRNA-seq and Bulk-tissue RNA-seq data is available at https://www.synapse.org/#!Synapse:syn22079621/
#####################################################

##Get the metacells for Proximal Tubule (PT) from the scRNA-seq data
targets.PT=metacell_seurat@meta.data
group=factor(targets.PT$exp.time,c("Control", "IRI_short_1d", "IRI_short_3d", "IRI_short_14d", "IRI_long_1d", "IRI_long_3d", "IRI_long_14d"))
datExpr.Cluster <- as.data.frame(GetAssayData(metacell_seurat, assay='RNA', slot='data')[VariableFeatures(metacell_seurat),])
datExpr.Cluster=as.data.frame(t(datExpr.Cluster))
rm(metacell_seurat)
gnS=colnames(datExpr.Cluster)
cat(length(gnS),'\n')

datExpr.Cluster=datExpr.Cluster[,match(gnS,colnames(datExpr.Cluster))]

# Make a multi-Expression data list containing all the data
nSets=1
setLabels=c("Cluster.scSeq.PT")
shortLabels=setLabels
multiExpr=vector(mode="list",length=nSets)
multiExpr[[1]] = list(data=as.data.frame(datExpr.Cluster)) #PTfreshsubset meta_cell
names(multiExpr[[1]]$data)=colnames(datExpr.Cluster)
rownames(multiExpr[[1]]$data)=rownames(datExpr.Cluster)
checkSets(multiExpr) # check data size
multiMeta=list(Cluster.NucSeq.PT=list(data=targets.PT))


## Network Construction
# Choose a set of soft-thresholding powers
powers = c(seq(1,10,by=1), seq(12,30, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 100,networkType="signed",corFnc="bicor")[[2]]);
saveRDS(powerTables, "powerTables_all_by_exp.time.rds")

# Plot the results:
pdf("1_Power_all_by_exp.time.pdf", height=10, width=18)
colors = c("blue", "red","black")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "mean connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off()
save(list=ls(),file="ClusterPT_all_by_exp.time.rda")



## Consensus WGCNA
softPower=18 ## power was chosen based on figure above
net=blockwiseConsensusModules(multiExpr, blocks = NULL,
                              randomSeed = 12345,
                              corType = "pearson",
                              power = softPower,
                              consensusQuantile = 0.3,
                              networkType = "signed",
                              TOMType = "unsigned",
                              TOMDenom = "min",
                              scaleTOMs = TRUE, scaleQuantile = 0.8,
                              sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                              useDiskCache = TRUE, chunkSize = NULL,
                              deepSplit = 4,
                              pamStage=FALSE,
                              detectCutHeight = 0.995, minModuleSize = 50,
                              mergeCutHeight = 0.2,
                              saveConsensusTOMs = TRUE,
                              consensusTOMFilePattern = "ConsensusTOM-block.%b_all_by_exp.time.rda")
save(list=ls(),file="ClusterPT_all_by_exp.time_2.rda")

consMEs = net$multiMEs;
moduleLabels = net$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
table(moduleColors)
consTree = net$dendrograms[[1]];
##UCI
MEs=moduleEigengenes(multiExpr[[1]]$data, colors = moduleColors, nPC=1)$eigengenes
MEs=orderMEs(MEs)
datExpr.UCI=multiExpr[[1]]$data
meInfo<-data.frame(rownames(datExpr.UCI), MEs)
colnames(meInfo)[1]= "SampleID"
KMEs<-signedKME(datExpr.UCI, MEs,outputColumnName = "kME",corFnc = "bicor")

ensembl=na.omit(colnames(datExpr.UCI))
geneInfo=as.data.frame(cbind(moduleColors, KMEs))
geneInfo2=as.data.frame(cbind(moduleLabels, KMEs))

# merged gene symbol column
colnames(geneInfo)[1]= "Initially.Assigned.Module.Color"
colnames(geneInfo2)[1]= "Initially.Assigned.Module.Color"
write.csv(geneInfo,file=paste('geneInfoSigned_UCI_all_by_exp.time.csv',sep=''))
write.csv(geneInfo2,file=paste('geneInfoSigned_UCI_all_by_exp.time2.csv',sep=''))

PCvalues.UCI=MEs
MEs.snRNAseq=moduleEigengenes(multiExpr[[1]]$data, colors = moduleColors, nPC=1)$eigengenes
MEs.snRNAseq=orderMEs(MEs.snRNAseq)
PCvalues.snRNAseq=MEs.snRNAseq
PCvalues.snRNAseq=PCvalues.snRNAseq[,match(colnames(PCvalues.UCI),colnames(PCvalues.snRNAseq))]



##Final saving
save(list=ls(),file="ClusterPT_all_by_exp.time_3.rda")



#================== Compute module scores in each metacell for all scWGCNA modules =========== 

# load scWGCNA modules #########################################################
cur_celltype <- "allcells" 
n_genes <- 25
geneInfo$GeneSymbol <- rownames(geneInfo)
modules <- unique(geneInfo$Initially.Assigned.Module.Color)
modules <- modules[order(modules)]
module_list <- lapply(modules, function(mod){
  cur <- subset(geneInfo, Initially.Assigned.Module.Color == mod)
  cur[,c('GeneSymbol', paste0('kME', mod))] %>%
    top_n(n_genes) %>% .$GeneSymbol
})
names(module_list) <-modules

# compute module scores:
metacell_seurat <- AddModuleScore(
  metacell_seurat,
  features=module_list,
  pool = rownames(metacell_seurat), k=F, nbin=24,
  name=paste0(cur_celltype, '_module')
)










