counts <- read.csv("bewerkt_countmatrix.csv", row.names = 1)

treatment <- c("control", "control", "control", "control", "reuma", "reuma", "reuma", "reuma")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- c("control1", "control2", "control3", "control4", "reuma1", "reuma2", "reuma3", "reuma4")


colnames(count_matrix)
rownames(treatment_table)



if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")


BiocManager::install("DESeq2")
library(DESeq2)

BiocManager::install("KEGGREST")
library(KEGGREST)

# Maak DESeqDataSet aan
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = treatment_table,
                              design = ~ treatment)

# Voer analyse uit
dds <- DESeq(dds)
resultaten <- results(dds)

# Resultaten opslaan in een bestand
#Bij het opslaan van je tabel kan je opnieuw je pad instellen met `setwd()` of het gehele pad waar je de tabel wilt opslaan opgeven in de code.

write.table(resultaten, file = 'Resultatencasusreuma.csv', row.names = TRUE, col.names = TRUE)

sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)

#Welke genen springen eruit? 
#Nu sorteren we het resultaat om te kijken naar de opvallendste genen:
hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]

#Bekijk nu welke genen het belangrijkst zijn volgens de analyse.
head(laagste_p_waarde)

if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}

install.packages("EnhancedVolcano")

library(EnhancedVolcano)

EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')

# Alternatieve plot zonder p-waarde cutoff (alle genen zichtbaar)
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0)

#opslaan figuur 
dev.copy(png, 'Volcanoplotcasusreuma.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

if (!requireNamespace("pathview", quietly = TRUE)) {
  BiocManager::install("pathview")
}

install.packages("pathview")

library(pathview)

resultaten[1] <- NULL
resultaten[2:5] <- NULL


gene_vector <- resultaten$log2FoldChange
names(gene_vector) <- rownames(resultaten)


pathview(
  gene.data = gene_vector,
  pathway.id = "hsa05323",  # KEGG ID voor Biofilm formation – E. coli
  species = "hsa",          # 'eco' = E. coli in KEGG
  gene.idtype = "SYMBOL"  # Kleurbereik voor log2FC van -5 tot +5
)


# 4. Genereer het PNG-bestand voor de TCA-cyclus
pathview(
  gene.data = gene_vector,
  pathway.id = "eco02026",  # KEGG ID voor Biofilm formation – E. coli
  species = "eco",          # 'eco' = E. coli in KEGG
  gene.idtype = "KEGG",     # Geef aan dat het KEGG-ID's zijn
  limit = list(gene = 5), out.suffix = "TCA_cycle"    # Kleurbereik voor log2FC van -5 tot +5
)

#Met de functie keggLink() uit het KEGGREST-pakket kun je in R koppelingen
#leggen tussen genen en pathways in de KEGG-database.

keggLink("pathway", "eco:b0221")


# Zoek pathways waarin een humaan gen voorkomt (bijv. hsa:7157 = TP53)
keggLink("pathway", "hsa:7157")


keggLink("eco", "path:eco00010")

# KEGG pathway ID voor reuma is hsa05323
keggLink("hsa", "path:hsa05323")

