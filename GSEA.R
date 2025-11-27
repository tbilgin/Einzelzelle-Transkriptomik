
install.packages("enrichR")
library(Seurat)
library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)

pathways_hallmark <- msigdbr::msigdbr(species = "Mus musculus",
                                      category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)


# Run GSEA for all clusters
for (cluster_id in clusters) {
  cat("\n========== CLUSTER", cluster_id, "==========\n")
  
  de_results <- FindMarkers(seu, ident.1 = cluster_id, logfc.threshold = 0)
  
  ranked_genes <- setNames(de_results$avg_log2FC, rownames(de_results))
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  

  
  gsea_results <- fgsea(pathways = pathways_hallmark,
                        stats = ranked_genes,
                        eps = 0,
                        minSize = 5,
                        maxSize = 500)
  
  gsea_results$cluster <- cluster_id
  all_gsea_results[[as.character(cluster_id)]] <- gsea_results
  
  print(head(gsea_results[order(padj)], 10))
}

# Combine results
all_results <- bind_rows(all_gsea_results)


# Visualization: Top pathways per cluster 
# =======================================

# Get top 5 pathways per cluster (by adjusted p-value)
top_per_cluster <- all_results %>%
  group_by(cluster) %>%
  slice_min(padj, n = 5) %>%
  ungroup()

ggplot(top_per_cluster, aes(x = factor(cluster), y = reorder(pathway, NES))) +
  geom_tile(aes(fill = NES)) +
  geom_text(aes(label = round(NES, 2)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  labs(title = "Top Enriched GO Pathways per Cluster",
       x = "Cluster",
       y = "Pathway",
       fill = "NES") +
  theme(axis.text.y = element_text(size = 8))


for (cluster_id in clusters) {
  
  cluster_data <- all_results %>%
    filter(cluster == cluster_id, padj < 0.05, NES > 0) %>%
    arrange(desc(NES)) %>%
    head(10)  # Top 10 enriched pathways
  
  if (nrow(cluster_data) == 0) {
    cat("\nNo significant enriched pathways for cluster", cluster_id, "\n")
    next
  }
  
  p <- ggplot(cluster_data, aes(x = reorder(pathway, NES), y = NES)) +
    geom_col(fill = "red", alpha = 0.7) +
    coord_flip() +
    theme_minimal() +
    labs(title = paste("Cluster", cluster_id, "- Top Enriched Pathways"),
         x = "Pathway",
         y = "NES")
  
  print(p)
}





