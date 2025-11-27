
# Umgebung:Pakete laden
#######################

install.packages("dplyr")
install.packages("SeuratObject")
install.packages("Seurat")
install.packages("patchwork")

library(data.table)
library(Matrix)
library(Seurat)
library(ggplot2)
library(dplyr)

# 1) create seurat object

seu <- CreateSeuratObject(
  counts = sparse_mat,
  project = "my_sample"
)

# 2) Calculate percentage of mitochondrial genes
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

# 3) calculate nFeature_RNA, nCount_RNA
meta <- seu@meta.data
meta$barcode <- rownames(meta)

# 4) Einfacher Violinplot: cells (nFeature_RNA) vs counts (nCount_RNA)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)

# 5) Filter and plot again(was für Filtern machen Sinn?)

seu <- subset(
  seu,
  subset =
    nFeature_RNA > 2000 &
    nFeature_RNA < 6000 &
    percent.mt < 70
)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 6) Normalize data
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

# 7) Hoch-variable Gene identifizieren und visualisieren
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
VariableFeaturePlot(seu)

# 8) Daten skalieren
all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)

# 9) Dimensionsreduktion
seu <- RunPCA(seu, features = VariableFeatures(object = seu))

# 10) Dimensionalität des Datensatzes bestimmen (wie viele Dimensionen brauchen wir?)
ElbowPlot(seu)

# Extract the standard deviations from your Seurat object
stdev <- seu@reductions$pca@stdev

# Calculate variance explained by each PC
var_explained <- (stdev^2) / sum(stdev^2)

# Create data frame
df <- data.frame(
  PC = 1:length(var_explained),
  Variance = var_explained
)

# Create plot
ggplot(df, aes(x = PC, y = Variance)) +
  geom_point() +
  geom_line() +
  labs(x = 'PC', 
       y = 'Explained Variance') 


# 11) Clustering the Cells
seu <- FindNeighbors(seu, dims = 1:5)
seu <- FindClusters(seu, resolution = 0.5)

#Nonlinear Dimensional Reduction with UMAP and t-SNE
seu <- RunUMAP(seu, dims = 1:5)
DimPlot(seu, reduction = "umap")

seu <- RunTSNE(seu, dims = 1:8)
DimPlot(seu, reduction = "tsne")

# 12) Find Markers
seu.markers <- FindAllMarkers(seu, only.pos = TRUE)

top3_per_cluster <- seu.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05 )  %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%  # sort within each cluster
  slice_head(n = 3)                        # Top 10 per Cluster

# Distribution via clusters
VlnPlot(seu, features = top3_per_cluster$gene, ncol = 3)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



