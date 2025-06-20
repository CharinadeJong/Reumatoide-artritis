# KEGG-analyse op basis van DESeq2-resultaten
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)

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






