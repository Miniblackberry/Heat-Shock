
install.packages("Seurat")

library(Seurat)
library(dplyr)

SR_hs00.data <- Read10X(data.dir = "data/SR_hs00/")
SR_hs60.data <- Read10X(data.dir = "data/SR_hs60/")

SR_hs00.seurat <- CreateSeuratObject(counts = SR_hs00.data, project = "SR_hs")
SR_hs00.seurat$stim <- "hs00"

SR_hs60.seurat <- CreateSeuratObject(counts = SR_hs60.data, project = "SR_hs")
SR_hs60.seurat$stim <- "hs60"

heat_shock <- merge(SR_hs00.seurat, y = SR_hs60.seurat, add.cell.ids = c("hs00", "hs60"), project = "heat_shock")

heat_shock[["percent.mt"]] <- PercentageFeatureSet(heat_shock, pattern = "^MT-")

VlnPlot(heat_shock, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

heat_shock <- subset(heat_shock, subset = nFeature_RNA > 0 & nFeature_RNA < 12000 & percent.mt < 10)

heat_shock <- NormalizeData(heat_shock, normalization.method = "LogNormalize", scale.factor = 10000)

heat_shock <- FindVariableFeatures(heat_shock, selection.method = "vst", nfeatures = 12000)

top10 <- head(VariableFeatures(heat_shock), 10)

plot1 <- VariableFeaturePlot(heat_shock)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

all.genes <- rownames(heat_shock)
heat_shock <- ScaleData(heat_shock, features = all.genes)

heat_shock <- RunPCA(heat_shock, features = VariableFeatures(object = heat_shock))
print(heat_shock[["pca"]], dims = 1:10, nfeatures = 10)


DimPlot(heat_shock, dims = 1:2, reduction = "pca", group.by = "stim", cols = c("blue","red"))

ElbowPlot(heat_shock)


heat_shock<- FindNeighbors(heat_shock, dims = 1:6)


heat_shock<- FindClusters(heat_shock, resolution = 0.5)

heat_shock <- RunUMAP(heat_shock, dims = 1:6)

DimPlot(heat_shock, reduction = "umap")
 DimPlot(heat_shock, reduction = "umap", group.by = "stim", cols = c("blue","red"))

heat_shock <- JoinLayers(heat_shock)

cluster1.0.markers <- FindMarkers(heat_shock, ident.1 = 1)
head(cluster1.0.markers, n = 10)

cluster3.0.markers <- FindMarkers(heat_shock, ident.1 = 3)
head(cluster3.0.markers, n = 10)

FeaturePlot(heat_shock, features = c("ACTB", "GAPDH"))

FeaturePlot(heat_shock, features = c("NSD3", "MAN1A2", "AKAP9", "TNPO3", "RTF1", "PLEKHB2", "SARNP", "TNRC6B", "ELF2", "CTBP2"))

FeaturePlot(heat_shock, features = c("HSPA1A", "HSPA1B"))

FeaturePlot(heat_shock, features = c("FOSL1", "ID1"))

heat_shock.markers <- FindAllMarkers(heat_shock, only.pos = TRUE)
heat_shock.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

heat_shock.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
DoHeatmap(heat_shock, features = top10$gene)
