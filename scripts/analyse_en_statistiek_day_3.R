setwd("C:/Users/dejon/OneDrive - NHL Stenden/Transcriptomics J2P4 CASUS reuma/Data_RA_raw")

counts <- read.csv("bewerkt_countmatrix.csv", row.names = 1)

treatment <- c("control", "control", "control", "control", "reuma", "reuma", "reuma", "reuma")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- c("control1", "control2", "control3", "control4", "reuma1", "reuma2", "reuma3", "reuma4")

count_matrix <- read.csv("bewerkt_countmatrix.csv")

colnames(count_matrix)
rownames(treatment_table)


# Laad benodigde packages
if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("DEseq2")
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

# Laad benodigde packages
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

# Filter significante genen
res_df <- as.data.frame(resultaten)
res_df$gene <- rownames(res_df)
sig_genes <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)

# Zet gen-symbolen om naar Entrez-ID's
entrez_ids <- bitr(sig_genes$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# KEGG verrijkingsanalyse
kegg_res <- enrichKEGG(gene = entrez_ids$ENTREZID,
                       organism = 'hsa',
                       pvalueCutoff = 0.05)

# Resultaten bekijken
head(kegg_res)

# Visualisatie
dotplot(kegg_res, showCategory = 10, title = "KEGG Pathway Enrichment") +
  scale_fill_gradientn(colors = c("mediumorchid4", "#CD1076","#ae017e", "#c51b8a","#dd3497","#FF1493","#fa9fb5", "#fbb4d9"), name = "p.adjust" )+
  labs(y = "Pathways")


# Resultaten omzetten naar data frame
kegg_df <- as.data.frame(kegg_res)

# Opslaan als CSV-bestand
write.csv(kegg_df, file = "kegg_results.csv", row.names = FALSE)

# Optioneel: eerste paar rijen bekijken
head(kegg_df)


pathview(
  gene.data = gene_vector,
  pathway.id = "hsa04064",  # KEGG ID voor 
  species = "hsa",          
  gene.idtype = "SYMBOL"  
)


keggLink("hsa", "path:hsa04064")





