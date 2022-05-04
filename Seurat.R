library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(rliger)
library(harmony)
library(dplyr)
sObj_FAPS_DMD_merge_081021 <- readRDS("/Volumes/Madre/Articulo Single cell RNA seq/sObj_FAPS_DMD_merge_081021.rds")
DimPlot(sObj_FAPS_DMD_merge_081021,label = T) + NoLegend()
#podemos cambiar la resolución para tener menos grupos
sObj_FAPS_DMD_merge_081021 <- FindClusters(sObj_FAPS_DMD_merge_081021, resolution = 0.4)
DimPlot(sObj_FAPS_DMD_merge_081021,label = T) + NoLegend()
new.cluster.ids <- c("Myofibroblasts", "Proliferative FAPs", "Regulatory FAPs", "Proliferative FAPs", "Proliferative FAPs", "Inflammatory FAPs", "Myofibroblasts")
dmdnewclusters <- sObj_FAPS_DMD_merge_081021
names(new.cluster.ids) <- levels(dmdnewclusters)
dmdnewclusters <- RenameIdents(dmdnewclusters, new.cluster.ids)
DimPlot(dmdnewclusters)
DimPlot(dmdnewclusters, split.by = "group")
#construir tablas que presenten la distribución de las células
count_table <- table(dmdnewclusters@meta.data$seurat_clusters, dmdnewclusters@meta.data$orig.ident)
count_table
#marcadores por cluster
cluster1.markers <- FindMarkers(dmdnewclusters, ident.1 = "Myofibroblasts", min.pct = 0.25)
head(cluster1.markers)
dmd.markers <- FindAllMarkers(dmdnewclusters, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
#ordenar el CSV de mayor a menor
dmd.markers[order(dmd.markers$cluster, -rank(dmd.markers$avg_log2FC), decreasing = FALSE), ]
dmd.markersordered <- dmd.markers[order(dmd.markers$cluster, -rank(dmd.markers$avg_log2FC), decreasing = FALSE), ]
write.csv(dmd.markersordered,"//Volumes/Madre/Articulo Single cell RNA seq/Graficos con menos clusters/dmd_markers.csv", row.names = TRUE)
#markers y heatmap
dmd.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(dmdnewclusters, features = top10$gene) + NoLegend()
#Añadir estadistica
library(ggpubr)
VlnPlot(dmdnewclusters, features = "MYL9", pt.size = 0.1, group.by = "group") + stat_compare_means()
#estudiar expresión de genes diferenciada entres las distintas poblaciones
Idents(dmdnewclusters) <- "group"
avg.t.cells<- as.data.frame(log1p(AverageExpression(dmdnewclusters, verbose = FALSE)$RNA))
avg.t.cells$gene <- rownames(avg.t.cells)
p1 <- ggplot(avg.t.cells, aes(DMD, paedCtl)) + geom_point()
p1
#si queremos hacer un volcanoplot
genes.increased <- FindMarkers(a, ident.1 = "DMD", ident.2 = "paedCtl", verbose = FALSE)
head(genes.increased, n=15)
a.genesordered <- genes.increased[order(-rank(genes.increased$avg_log2FC), decreasing = FALSE), ]
write.csv(a.genesordered,"//Volumes/Madre/Articulo Single cell RNA seq/Graficos con menos clusters/agenesordered.csv", row.names = TRUE)
library(EnhancedVolcano)
library(magrittr)
EnhancedVolcano(a.genesordered, lab = rownames(a.genesordered), x = 'avg_log2FC', y= "p_val_adj")
#elegir un subset determinado
a <- subset(dmd2, idents = "Myofibroblasts")




