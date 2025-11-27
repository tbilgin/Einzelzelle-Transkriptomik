# Bioinformatik & Omics Modul 2025
## 7 & 21.11.25 single Cell Transcriptomics und Pathway Enrichment

### Lernziele
- die Unterschiede zwischen Einzelzell-Sequnenzieren und Bulk-Sequenzieren beschreiben: Anwendung, Methode, Analyse 
- Technologien zur Erfassung einzelner Zellen listen und die entsprechende Ansätze zueinander vergleichen
- Schritte in Einzelzelle-Analyse motivieren, diese in R ausführen und Ausgaben interpretiere
- Einzelzell-Transkriptomik-Matrix erstellen
- Differentielle Expressionsanalyse mit PCA ausführen
- Clusters mit UMAP und t-SNE bilden um die Katogorien zu identifizieren
- Gene und Pathway Enrichment Analyse durchführen
- Zwecke von räumlichen Einzelzell-Sequenzieren beschreiben

## Ablauf

7.11: 
Theorie von Einzelzell-Sequnenzieren (Technologien, Workflow)
Theorie von Hauptkomponentenanalyse und Dimensioinsreduktion in Transkriptomik
Praktikum in R mit einem 10X Tutorial: Single_Cell_RNA_Analysis

21.11:
Theorie von Clusterbildung mit UMAP und t-SNE
Theorie von räumlichem Einzelzell-Sequenzieren
Praktikum in R mit smart-Seq Datensatz:
    1. Datenvorbereitung: Benutze Kode prepare_data.R um den Hausdatensatz (nicht veröffentlich) vozubereiten, gene_mapping_biomart_mouse.csv hat die Ensembl Gene und entsprechende Gene Symbol IDs
    2. Datenbearbeitung um die Katogorien zu identifizieren: Nieren.R
    3. Pathway Enrichment mit SingleCell Transcriptomics: GSEA.R




