library(SingleCellExperiment)
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(doMC)
library(R2HTML)
library(rbokeh)
set.seed(123)
options(width=200)

setwd('~/SCENIC/')
dir.create("int")

### ==================================================================================
### ========================== Extract information from sce object ===================
### ==================================================================================
load('/slingshot/sce')

#extract data matrix
exprMat <- sce@assays@data@listData$counts
saveRDS(exprMat, 'int/exprMat.rds')

#extract cell annotation info
cellInfo <- colData(sce)$GMM
GMMclusters <- cellInfo
cellInfo <- data.frame(GMMclusters)
cellInfo$GMMclusters <- as.factor(cellInfo$GMMclusters)
saveRDS(cellInfo, file="int/cellInfo.Rds")

#extract embeddings
dr_coords <- sce$rd2
saveRDS(rd2, file="int/dr_coords.Rds")


### Load data
exprMat <- readRDS('int/exprMat.rds')
exprMat <- as.matrix(exprMat)


### Initialize SCENIC settings
org <- "mgi" # or hgnc, or dmel
dbDir <- "~/~/" # RcisTarget databases location
myDatasetTitle <- "IRI" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=12) 

# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
scenicOptions@inputDatasetInfo$dr_coords <- "int/dr_coords.Rds"

# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 



### ====================================================================================
### ========================== Co-expression network ===================================
### ====================================================================================
# 1. (Adjust minimum values according to your dataset)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
# 2. Before proceeding to the network inference, check whether any known relevant genes are filtered-out
interestingGenes <- c("Cd74", "Havcr1", "Il1b")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

# 3. We can now **filter the expression matrix** to contain only these `r length(genesKept)` genes. 
#This matrix is now ready for the co-expression analysis.
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
saveRDS(exprMat_filtered, 'int/exprMat_filtered.rds')

# 4. Correlation
#GENIE3/GRNBoost can detect both positive and negative associations. In order to distinguish potential activation from repression, we will split the targets into 
#positive- and negative-correlated targets (i.e. Spearman correlation between the TF and the potential target).
#*(This step can be run either before/after or simultaneously to GENIE3/GRNBoost)*
runCorrelation(exprMat_filtered, scenicOptions)

# 5. GENIE3 
runGenie3(exprMat_filtered, scenicOptions)

# 6. Run the remaining steps using the *wrapper* functions: 
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 6
scenicOptions@settings$seed <- 123

### Build and score the GRN
exprMat_filtered <- as.matrix(exprMat_filtered)
exprMat_log <- log2(exprMat_filtered+1)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
runSCENIC_4_aucell_binarize(scenicOptions)



###=====================================================
###============== Exploring output =====================
###=====================================================
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
print(tsneFileName(scenicOptions))
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")


###============== Regulators for clusters or known cell types =====================
#The regulatory analysis from SCENIC can be combined with other analyses, for example clustering, or focus on regulators for specific cell types. 
#There are multiple options to do these analyses (your imagination is the limit!). Here are some quick examples to start:

#-------- by GMMcluster ---------
#Average Regulon Activity by cluster (Clusters could also be used instead of “cell types”, e.g. with Seurat: cellInfo <- data.frame(seuratCluster=Idents(seuratObject)))
cellInfo <- readRDS("int/cellInfo.Rds")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$GMMclusters),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
write.table(regulonActivity_byCellType_Scaled, 'output/Step4_RegulonActivity_byCellType.tsv', sep='\t', quote=F)
Cairo::CairoPDF("output/Step4_RegulonActivity_PHeatmap_byCellType.pdf", height=30, width=10)
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
write.table(topRegulators, 'output/Step4_RegulonActivity_topRegulators.tsv', sep='\t', quote=F)



