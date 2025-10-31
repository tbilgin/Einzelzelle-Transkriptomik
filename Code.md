# Single-Cell RNA-Seq Analyse mit Seurat

Ein Schritt-für-Schritt Tutorial zur Analyse von Einzelzell-Transkriptomdaten mit dem PBMC3k Datensatz.

---

## 1. Arbeitsumgebung einrichten

### Pakete installieren

```r
install.packages("dplyr")
```

```r
install.packages("SeuratObject")
```

```r
install.packages("Seurat")
```

```r
install.packages("patchwork")
```

### Bibliotheken laden

```r
library(dplyr)
library(SeuratObject)
library(Seurat)
library(patchwork)
```

### Datensatz laden

Laden Sie den PBMC-Datensatz von: `https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz`

Entpacken Sie die Datei und stellen Sie sicher, dass der Pfad zum Datenverzeichnis korrekt ist.

```r
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
```

### Seurat-Objekt initialisieren

```r
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
```

---

## 2. Vorverarbeitungs-Workflow

### Qualitätskontrolle und Filterung

Mitochondriale RNAs und niedrige Counts weisen auf Probleme bei der Probenverarbeitung und Zellgesundheit hin.

```r
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```

### QC-Metriken visualisieren

```r
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

### Zellen filtern

Nach Betrachtung der Datenverteilung filtern wir Zellen, die weniger als 200 RNA-Typen oder mehr als 2500 exprimieren, sowie Zellen mit mehr als 5% mitochondrialen Beiträgen.

```r
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

### Daten normalisieren

Nach der Qualitätskontrolle werden die Daten auf log-Skala normalisiert.

```r
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```

### Hoch-variable Gene identifizieren

Identifizieren Sie die am stärksten variablen Features (d.h. die Top 2000 variablen Gene).

```r
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
```

### Die 10 variabelsten Gene identifizieren

```r
top10 <- head(VariableFeatures(pbmc), 10)
```

### Variable Features plotten

```r
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```

### Daten skalieren

Vor der Dimensionsreduktion skalieren wir die Daten, damit stark exprimierte Gene nicht dominieren.

```r
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

---

## 3. Dimensionsreduktion

Vereinfachung des Datensatzes unter Beibehaltung der wesentlichen Komponenten.

### PCA durchführen

```r
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
```

### PCA-Ergebnisse untersuchen und visualisieren

```r
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
```

### Dimensionalität des Datensatzes bestimmen

Welche Dimensionen ermöglichen es dem Datensatz, seine Komplexität zu behalten und gleichzeitig zu vereinfachen?

```r
ElbowPlot(pbmc)
```

**Interpretation:** Es gibt wenig bis keine Variation um die Dimensionen 9-10 herum bei der Hauptkomponentenanalyse.

---

## 4. Clustering der Zellen

### Nachbarn finden

```r
pbmc <- FindNeighbors(pbmc, dims = 1:10)
```

### Cluster identifizieren

```r
pbmc <- FindClusters(pbmc, resolution = 0.5)
```

### Nichtlineare Dimensionsreduktion mit UMAP

```r
pbmc <- RunUMAP(pbmc, dims = 1:10)
```

### Einzelne Cluster visualisieren

```r
DimPlot(pbmc, reduction = "umap")
```

---

## 5. Differenziell exprimierte Features als Cluster-Biomarker finden

### Marker für alle Cluster finden

```r
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
```

### Marker mit hoher log2-Faltänderung filtern

```r
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
```

### Spezifität der Marker messen

```r
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
```

### Differentielle Genexpression zwischen Clustern visualisieren

```r
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
```

### Spezifische Features über Cluster hinweg visualisieren

```r
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
```

---

## 6. Zelltyp-Identitäten den Clustern zuweisen

### Heatmap für Top-10-Marker generieren

```r
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
```

```r
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
```

**Hinweis:** Überprüfen Sie in der Literatur, welche Marker welchem Zelltyp entsprechen. Es gibt auch zusätzliche Programme, die dabei helfen.

### Cluster umbenennen

```r
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
```

### Finale UMAP-Visualisierung mit Zelltyp-Labels

```r
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

---

## 7. Plot speichern

### ggplot2 laden

```r
library(ggplot2)
```

### Finalen Plot erstellen und speichern

```r
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + 
  xlab("UMAP 1") + 
  ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size = 10)))

ggsave(filename = "pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)
```

---

## Zusammenfassung der Analyse-Schritte

1. ✅ **Daten laden und QC**: Qualitätskontrolle und Filterung problematischer Zellen
2. ✅ **Normalisierung**: Log-Normalisierung für vergleichbare Expressionswerte
3. ✅ **Feature-Selektion**: Identifikation hoch-variabler Gene
4. ✅ **Dimensionsreduktion**: PCA zur Erfassung der Hauptvariationsquellen
5. ✅ **Clustering**: Identifikation von Zellpopulationen
6. ✅ **Marker-Identifikation**: Finden charakteristischer Gene für jeden Cluster
7. ✅ **Annotation**: Zuordnung biologischer Zelltypen zu Clustern

---

## Weiterführende Ressourcen

- [Seurat Dokumentation](https://satijalab.org/seurat/)
- [10x Genomics Datensätze](https://www.10xgenomics.com/resources/datasets)
- [PBMC3k Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
