library(dplyr)
library(Seurat)
library(data.table)
set.seed(123)

setwd('~/Seurat/')
object <- readRDS("/~/data/LIGER/object_res0.25.rds")
#store LIGER_clusters
object@meta.data$LIGER_clusters <- Idents(object)
#reorder factor levels for LIGER_clusters
object@meta.data$LIGER_clusters <- factor(object@meta.data$LIGER_clusters, 
                                          levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

table_Idents <- table(object@meta.data$LIGER_clusters)
write.table(x = table_Idents, row.names = TRUE, file = "/~/LIGER_clusters_res0.25.csv")



#================= Add experimental setup to metadata =================
#add exp.cond
object@meta.data$exp.cond <- plyr::mapvalues(
  x = object@meta.data$orig.ident,
  from = c("norm1", "norm2", "norm3", "norm4", "norm5", "Scl", 
           "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", "IRF7", "IRF8", "IRF9", "IRF10", "IRF11", "IRF12"),
  to = c("Control", "Control", "Control", "Control", "Control", "Control", 
         "IRI_long", "IRI_short", "IRI_short", "IRI_long", "IRI_long", "IRI_short", "IRI_long", "IRI_short", "IRI_long", "IRI_short", "IRI_short", "IRI_long"))

#add time.pod
object@meta.data$time.pod <- plyr::mapvalues(
  x = object@meta.data$orig.ident,
  from = c("norm1", "norm2", "norm3", "norm4", "norm5", "Scl", 
           "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", "IRF7", "IRF8", "IRF9", "IRF10", "IRF11", "IRF12"),
  to = c("Control", "Control", "Control", "Control", "Control", "Control", 
         "1d", "1d", "3d", "3d", "1d", "1d", "3d", "14d", "14d", "14d", "3d", "14d"))

#add exp.time
object@meta.data$exp.time <- factor(paste(object@meta.data$exp.cond, object@meta.data$time.pod, sep="_"))
object@meta.data$exp.time <- plyr::mapvalues(
  x = object@meta.data$exp.time,
  from = c("Control_Control", "IRI_long_1", "IRI_long_14", "IRI_long_3", "IRI_short_1", "IRI_short_14", 
           "IRI_short_3"),
  to = c("Control", "IRI_long_1d", "IRI_long_14d", "IRI_long_3d", "IRI_short_1d", "IRI_short_14d", 
         "IRI_short_3d"))
object@meta.data$exp.time <- factor(object@meta.data$exp.time, levels = c("Control",
                                                                          "IRI_short_1d", "IRI_short_3d", "IRI_short_14d",
                                                                          "IRI_long_1d", "IRI_long_3d", "IRI_long_14d"))



#================= Scaling the data =================
all.genes <- rownames(object)
object <- ScaleData(object, features = all.genes)



#================== Recode cluster identities =======================
#add LIGER_clusters2
object$LIGER_clusters2 <- plyr::mapvalues(
  x = object$LIGER_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
  to = c("0_Healthy S3", "1_Healthy S3", "2_Healthy S2", "3_Healthy S1", "4_Healthy S1", "5_Healthy S1", "6_Healthy/LOH/DCT", "7_Maladaptive", "8_Injured", "9_Injured", "10_Injured")
)



#================= Investigate Cell cycle genes =================
library(tidyverse)
s <- readRDS("/~/mus_s_genes.RDS") 
g2m <- readRDS("/~/mus_g2m_genes.RDS") 
object <- CellCycleScoring(object, 
                           s.features = s, 
                           g2m.features = g2m, 
                           set.ident = F)



#================= Finding differentially expressed features (cluster biomarkers) =================
object.markers <- FindAllMarkers(object, test.use="MAST", only.pos = F, min.pct = 0.05, logfc.threshold = 0.2)
object.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) %>% print(n = 50) -> tibble_top50_bycluster_MAST
object.markers %>% group_by(cluster) -> tibble_FindAllMarkers_MAST
fwrite(x = tibble_top50_bycluster_MAST, row.names = TRUE, file = "/~/tibble_top50_bycluster_res0.25_MAST.csv")
fwrite(x = tibble_FindAllMarkers_MAST, row.names = TRUE, file = "/~/tibble_FindAllMarkers_res0.25_MAST.csv")




#================= Average feature expression across cell cluster identity, invert log and write out =================
cluster.averages <- AverageExpression(object, return.seurat = FALSE, 
                                      slot = "data", 
                                      verbose = TRUE)



#================= IRI Scoring =================
library(tidyverse)
PT_top100DEG_byexpcond <- read.csv("/~/PT_top100DEG_byexpcond.csv")
Control <- as.vector(PT_top100DEG_byexpcond$top100PT_Control)
IRI_long <- as.vector(PT_top100DEG_byexpcond$top100PT_IRI_long)
intersect(Control, IRI_long) #0 intersecting genes, as expected

#investigate whether expression of genes of these two gene lists correlate or anticorrelate
#load cluster averages for SusztakIRIdataset
avg_Susztak <- cluster.averages
colnames(avg_Susztak) <- c("X","PT_0","PT_1","PT_2","PT_3","PT_4","PT_5","PT_6","PT_7","PT_8","PT_9","PT_10")
rownames(avg_Susztak) <- avg_Susztak$X

#correlation
df <- avg_Susztak[c(Control, IRI_long),] #subset df to only those rows where genes are entailed in both lists 
df2 <- df[,2:12] #delete column X
df3 <- t(df2) #need to transpose so that correlation between genes is correlated, not between PT subclusters!
library("Hmisc")
cor_PCC <- rcorr(as.matrix(df3))
#exclude some genes which do not fulfill criteria for correlation or anticorrelation
Control_corr <- setdiff(Control, c("Rps27", "Rps29", "Rps28", "Mat2a", "Son"))
IRI_long_corr <- setdiff(IRI_long, c("Spp1", "Lcn2", "S100a6", "Plin2", "Tmsb10"))

#repeat correlation
df_corr <- avg_Susztak[c(Control_corr, IRI_long_corr),] #subset df to only those rows where genes are entailed in both lists 
df2_corr <- df_corr[,2:12] #delete column X
df3_corr <- t(df2_corr) #need to transpose so that correlation between genes is correlated, not between PT subclusters!
cor_PCC_corr <- rcorr(as.matrix(df3_corr))
Pvalues <- cor_PCC_corr$P

#scoring
object <- CellCycleScoring(object, 
                           s.features = Control_corr, 
                           g2m.features = IRI_long_corr, 
                           set.ident = F)

#rename phases to IRI score labeling
object@meta.data$IRI_label <- plyr::mapvalues(
  x = object@meta.data$Phase,
  from = c("G1", "G2M", "S"),
  to = c("neither", "IRI", "Control")
)

#reorder factor levels for IRI_label
object$IRI_label <- factor(object$IRI_label, levels = c("Control", "IRI", "neither"))
table(object$IRI_label)

#create new slot and copy IRI.Score
object@meta.data$IRI.Score <- object@meta.data$G2M.Score













