
# Umgebung:Pakete laden
#######################

install.packages("dplyr")
install.packages("SeuratObject")
install.packages("Seurat")
install.packages("patchwork")
install.packages("Matrix")
install.packages("data.table")
install.packages("rmarkdown")
install.packages("psych")

install.packages("remotes")
remotes::install_version("ggplot2", version = "3.5.2", repos = "https://cran.r-project.org")
packageVersion("ggplot2")
# muss version ggplot version 3.5.2 zeigen

# Wichtig: 
# jetzt R restarten, click: Save the changes