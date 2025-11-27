# 1) Setting up your workspace
#Installing the packages necessary for the Single Cell Analysis

#install.packages("dplyr")
#install.packages("SeuratObject")
#install.packages("Seurat")
#install.packages("patchwork")

library(data.table)
library(Matrix)
library(Seurat)
library(ggplot2)
library(dplyr)


# 1) load labeled matrix
labeld <- "Labeled_Matrix.csv"
dt <- fread(labeld)

# 2) name columns
setnames(dt, names(dt)[1:3], c("gene", "cellbarcode", "value"))

# 3) Clean-up
dt[, gene        := trimws(as.character(gene))]
dt[, cellbarcode := trimws(as.character(cellbarcode))]
dt[, value       := as.numeric(value)]

dt <- dt[
  !is.na(gene) & gene != "" &
    !is.na(cellbarcode) & cellbarcode != "" &
    !is.na(value)
]

# 4) sum duplicates
dt_aggr <- dt[, .(value = sum(value, na.rm = TRUE)), by = .(gene, cellbarcode)]


# 5) Sparse-Matrix
genes        <- unique(dt_aggr$gene)
cellbarcodes <- unique(dt_aggr$cellbarcode)  

i <- match(dt_aggr$gene,        genes)
j <- match(dt_aggr$cellbarcode, cellbarcodes)

sparse_mat <- sparseMatrix(
  i = i, j = j, x = dt_aggr$value,
  dims = c(length(genes), length(cellbarcodes)),
  dimnames = list(genes, cellbarcodes)
)

# 6) Mapping Ensembl IDs to Gene Symbols
mapping <- read.csv("gene_mapping_biomart_mouse.csv", stringsAsFactors = FALSE)

symbol_map <- setNames(mapping$symbol, mapping$ensembl_gene_id)
genes_mat <- rownames(sparse_mat)

new_names <- ifelse(
  genes_mat %in% names(symbol_map),
  symbol_map[genes_mat],
  genes_mat
)

new_names <- make.unique(new_names)
rownames(sparse_mat) <- new_names



