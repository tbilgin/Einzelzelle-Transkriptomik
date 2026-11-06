#---------------------------------------
# GO-Enrichment pro Cluster (BP)
#---------------------------------------

library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

# only significant, positively regulated markers
sig_markers <- seu.markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0)

# Symbol -> EntrezID (Mus musculus)
gene_mapping <- bitr(
  unique(sig_markers$gene),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

sig_markers <- sig_markers %>%
  inner_join(gene_mapping, by = c("gene" = "SYMBOL"))

# Cluster as a character, so there is no trouble with factors
sig_markers$cluster <- as.character(sig_markers$cluster)
clusters <- sort(unique(sig_markers$cluster))

# List of enrichResult objects + combined data frame
go_list <- list()
go_results_list <- list()

for (cl in clusters) {
  message("Running GO for cluster ", cl, " ...")
  
  genes_cl <- sig_markers %>%
    filter(cluster == cl) %>%
    pull(ENTREZID) %>%
    unique()
  
  ego <- enrichGO(
    gene          = genes_cl,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  # can be NULL or empty if nothing significant
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    go_list[[paste0("Cluster_", cl)]] <- ego
    
    df <- as.data.frame(ego)
    df$Cluster <- cl
    go_results_list[[cl]] <- df
  }
}

# All Cluster Tables together
go_results <- dplyr::bind_rows(go_results_list)


top_go_by_cluster <- go_results %>%
  dplyr::group_by(Cluster) %>%
  dplyr::arrange(p.adjust, .by_group = TRUE) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup() %>%
  dplyr::select(Cluster, Description, pvalue, p.adjust, geneID)

top_go_by_cluster


# Die Cluster Bennenen
# change the names in kidney locations
new.cluster.ids <- c(
  "a",
  "b",
  "c",
  "d"
)
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
DimPlot(seu, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
