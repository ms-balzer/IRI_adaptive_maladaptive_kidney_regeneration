rm(list=ls())
library(slingshot)
library(Seurat)
library(ggbeeswarm)
library(ggthemes)
library(SingleCellExperiment)
library(RColorBrewer)
set.seed(123)

setwd('~/Slingshot/')

#  ================ 1 LOAD DATA ================
#load Seurat object and convert to SingleCellExpxeriment object
PTtrajectory_A <- readRDS("Seurat/PTtrajectory_A.rds")
barcodes <- colnames(PTtrajectory_A)
table(Idents(PTtrajectory_A))
sce <- as.SingleCellExperiment(PTtrajectory_A)

# Add metadata from Seurat object to sce object
sce@colData$orig.ident <- PTtrajectory_A@ meta.data$orig.ident
sce@colData$orig.ident <- factor(sce@colData$orig.ident, levels = c("norm1", "norm2", "norm3", "norm4", "norm5", "Scl",
                                                                    "IRF2", "IRF6", "IRF3", "IRF11", "IRF8", "IRF10",
                                                                    "IRF1", "IRF5", "IRF4", "IRF7", "IRF9", "IRF12"))
sce@colData$exp.cond <- PTtrajectory_A@ meta.data$exp.cond
sce@colData$time.pod <- PTtrajectory_A@ meta.data$time.pod
sce@colData$IRI.Score <- PTtrajectory_A$IRI.Score
IRIscore_df <- data.frame(PTtrajectory_A$IRI.Score)

# create new meta.data column exp.time and sort levels
sce@colData$exp.time <- factor(paste(sce@colData$exp.cond,sce@colData$time.pod, sep="_"))
sce@colData$exp.time <- factor(sce@colData$exp.time, levels = c("Control_Control", 
                                                                "IRI_short_1", "IRI_short_3", "IRI_short_14", 
                                                                "IRI_long_1", "IRI_long_3", "IRI_long_14"))



#  ================ 2 UPSTREAM ANALYSIS ================
#  ==== 2.1 gene filtering ====
# filter genes down to potential cell-type markers
# at least M (15) reads in at least N (15) cells
geneFilter <- apply(assays(sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sce <- sce[geneFilter, ]


#  ==== 2.2 normalization ====
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)


# ==== 2.3 dimensionality reduction ====
# 2.3a) OPTION 1 (PCA)
pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]
save(rd1, file = 'rd1')

# 2.3b) OPTION 2 (diffusion maps with destiny)
library(destiny, quietly = TRUE)
dm <- DiffusionMap(t(log1p(assays(sce)$norm)))
rd2 <- cbind(DC1 = dm$DC1, DC2 = dm$DC2)
save(rd2, file = 'rd2')

# add both dimensionality reductions to the sce object
reducedDims(sce) <- SimpleList(PCA = rd1, DiffMap = rd2)


# ==== 2.4 Clustering Cells ====
# 2.4a) Add prior clustering information from Seurat to the sce object
sce@colData$Seurat_clusters <- as.character(Idents(PTtrajectory_A))
colData(sce)$Seurat_CLUSTERS <- Idents(PTtrajectory_A)
cl0 <- Idents(PTtrajectory_A)
save(cl0, file = 'cl0')

# 2.4b) OPTION 1 Gaussian mixture modeling
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd2)$classification #found unique(cl1) = 5 clusters
save(cl1, file = 'cl1')
colData(sce)$GMM <- cl1

# 2.4c) OPTION 2 k means (CAVE: need to arbitrarily choose n of centers!!!)
cl2 <- kmeans(rd2, centers = 6)$cluster
save(cl2, file = 'cl2')



#  ================ 3 USING SLINGSHOT ================
# ==== 3.1 RUN 1 ====
# Run Slingshot using the wrapper function "slingshot" (combining "getLineages" and "getCurves")
# with dimensionality reduction 'DiffMap' and cluster labels identified by 'GMM'.
sce <- slingshot(sce, reducedDim = 'DiffMap', clusterLabels = 'GMM')
summary(sce$slingPseudotime_1)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

pdf('slingshot_DiffMap_GMM.pdf')
plot(reducedDims(sce)$DiffMap, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')
legend("right", legend=unique(sce$GMM), 
       title="GMM Cluster", 
       col = unique(brewer.pal(9,'Set1')[sce$GMM]), 
       pch = 20, cex = 1, y.inters=0.7, x.inters=0.5)
dev.off()

# Plot Slingshot pseudotime vs cell stage. 
p1 <- ggplot(as.data.frame(sce@colData), aes(x = slingPseudotime_1, y = as.character(GMM), 
                                             colour = as.character(GMM))) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Cluster") +
  ggtitle("Lineage 1")
p2 <- ggplot(as.data.frame(sce@colData), aes(x = slingPseudotime_2, y = as.character(GMM), 
                                             colour = as.character(GMM))) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Cluster") +
  ggtitle("Lineage 2")
library(cowplot)
pdf('slingshotpseudotime_bycluster_DiffMap_GMM.pdf', width=10, height=5)
plot_grid(p1, p2)
dev.off()



# Save current progress.
save(sce, file = 'sce')



#  ================ 4 DOWNSTREAM SLINGSHOT FUNCTIONALITY ================
# ==== 4.1 Identifying global lineage structure ====
# 4.1a) no end clusters
#rd2=DiffMap, cl1=GMM, start.clus = 
lin1 <- getLineages(rd2, cl1, start.clus = '4')
lin1
save(lin1, file = 'lin1')
#load('lin1')

pdf('getLineages_DiffMap.rd2_GMM.cl1.pdf')
plot(rd2, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(lin1, lwd = 2, col = 'black')
legend("right", legend=unique(sce$GMM), 
       title="GMM Cluster", 
       col = unique(brewer.pal(9,'Set1')[sce$GMM]), 
       pch = 20, cex = 1, y.inters=0.7, x.inters=0.5)
dev.off()

# 4.1b) with known end clusters
lin2 <- getLineages(rd2, cl1, start.clus= '4', end.clus = c('5','1'))
lin2
save(lin2, file = 'lin2')

pdf('getLineages_DiffMap.rd2_GMM.cl1_knownendclusters.pdf')
plot(rd2, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(lin2, lwd = 2, col = 'black')
legend("right", legend=unique(sce$GMM), 
       title="GMM Cluster", 
       col = unique(brewer.pal(9,'Set1')[sce$GMM]), 
       pch = 20, cex = 1, y.inters=0.7, x.inters=0.5)
dev.off()

# ==== 4.2 Constructing smooth curves and ordering cells ====
crv1 <- getCurves(lin1)
crv1
save(crv1, file = 'crv1')



# ==== 4.3 Extract Pseudotime from crv ====
pt1 <- slingPseudotime(crv1, na=TRUE)
save(pt1, file = 'pt1')



#  ================ 5 IDENTIFY TEMPORALLY EXPRESSED GENES ================
require(gam)
t <- sce$slingPseudotime_1

# for time, only look at the 100 most variable genes
Y <- log1p(assays(sce)$norm)
var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
Y <- Y[var100,]

# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  suppressWarnings({
    tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})

#We can then pick out the top genes based on p-value and visualize their expression over developmental time with a heatmap.
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
heatdata <- assays(sce)$norm[topgenes, order(t, na.last = NA)]
heatclus <- sce$GMM[order(t, na.last = NA)]

pdf('DEG_pseudotime_heatmap.pdf', width=15, height=15)
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus],
        margins = c(5, 5))
dev.off()



