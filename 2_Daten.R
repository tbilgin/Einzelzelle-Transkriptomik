# install.packages(c("data.table","Matrix","Seurat","cowplot"))
# install.packages("rmarkdown")
# install.packages("psych")


library(data.table)
library(Matrix)
library(Seurat)
library(ggplot2)
library(dplyr)
library(psych)

# 1) Laden vom Datensatz
labeld <- "Labeled_Matrix.csv"
dt <- fread(labeld)

# 2) Namen der Spalten
setnames(dt, names(dt)[1:3], c("gene", "cellbarcode", "value"))


# 3) Was machen wir Null-Werte?
# Gibt es drin Nulll-Werte?
sum(dt$value == 0)
# Gibt es Werte für jedes Gen in jeder Zelle?
# check ein Gen:
sum(dt$gene == "ENSMUSG00000025917")
# check 2: nimm ein anderes Gen



# 4) erstelle ein Sparse-Matrix
genes        <- unique(dt$gene)
cellbarcodes <- unique(dt$cellbarcode)  

i <- match(dt$gene,        genes)
j <- match(dt$cellbarcode, cellbarcodes)

sparse_mat <- sparseMatrix(
  i = i, j = j, x = dt$value,
  dims = c(length(genes), length(cellbarcodes)),
  dimnames = list(genes, cellbarcodes)
)

# 5) Mapping Ensembl IDs to Gene Symbols
# Ensembl IDs 
mapping <- read.csv("gene_mapping_biomart_mouse.csv", stringsAsFactors = FALSE)

symbol_map <- setNames(mapping$symbol, mapping$ensembl_gene_id)
head(symbol_map)

# Haben alle Ensembl IDs ein GenSymbol?
sum(genes %in% names(symbol_map))
# Haben alle Ensembl IDs ein einzigartiges GenSymbol?
sum(duplicated(symbol_map))
mapping[duplicated(symbol_map),]
# drei Ensembl GenIDs haben alle Ndor1 als Symbol 
mapping[mapping$symbol == "Ndor1",]
# check ein anderes Gen:

# 6) Bennenen mit neuen Namen

# Wechsle die Namen wann möglich
#ifelse(test, 
#       yes, 
#       no)

new_names <- ifelse(
  genes %in% names(symbol_map),
  symbol_map[genes],
  genes
)
# duplizierte Zeilen sind versioniert.
### finde die duplizierte Zeilen:
which(duplicated(new_names))
new_names[duplicated(new_names)]
# versioniere die duplizierte Zeilen mit ".1" usw: 
new_names_uniq <- make.unique(new_names)
# wie sehen sie jetzt aus?
new_names_uniq[duplicated(new_names)]
# sparse_mat hat die aktualisierte Gennamen
rownames(sparse_mat) <- new_names_uniq



