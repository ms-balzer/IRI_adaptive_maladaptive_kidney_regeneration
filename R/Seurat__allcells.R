library(dplyr)
library(Seurat)
library(data.table)
set.seed(123)

object <- readRDS("/~/data/LIGER/object_res0.25.rds")
#store LIGER_clusters
object@meta.data$LIGER_clusters <- Idents(object)
#reorder factor levels for LIGER_clusters
object@meta.data$LIGER_clusters <- factor(object@meta.data$LIGER_clusters, 
                                                        levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15",
                                                                   "16", "17", "18", "19", "20", "21", "22", "23"))

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



#========= Recode cluster identities ===========
object@meta.data$LIGER_clusters3 <- plyr::mapvalues(
  x = object@meta.data$LIGER_clusters,
  from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15",
           "16", "17", "18", "19", "20", "21", "22", "23"),
  to = c("PT S3", "LOH", "PT S1", "DCT", "Mono_1", "Granul_2", "T cell_2", "GEC", "Granul_1", "PC", "T cell_1", "Granul_1", "Macro", "Macro", "Granul_1", "IC",
         "Macro", "B cell", "Endo", "B cell", "Myofib", "Mono_2", "Endo", "Podo")
)
Idents(object) <- object@meta.data$LIGER_clusters3
Idents(object) <- factor(Idents(object),
                         levels = c("GEC", "Endo", "Myofib", "Podo", "PT S1", "PT S3", "LOH", "DCT", 
                                    "PC", "IC", "Granul_1", "Granul_2", "Mono_1", "Mono_2", "Macro", "T cell_1", "T cell_2", "B cell"))



#================= Average feature expression across cell cluster identity, invert log and write out =================
cluster.averages <- AverageExpression(object, return.seurat = FALSE, 
                                      slot = "data", 
                                      verbose = TRUE) # result is a list called RNA












