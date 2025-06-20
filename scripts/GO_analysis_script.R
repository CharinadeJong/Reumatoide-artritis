
# 📦 Benodigde packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# 📥 Data inladen
setwd("C:/Users/dejon/OneDrive - NHL Stenden/Transcriptomics J2P4 CASUS reuma/Data_RA_raw")
sigs <- read.csv("Resultatencasusreuma.csv", sep = " ", header = TRUE, row.names = 1)

# 🔁 Genen filteren
genes_high <- rownames(sigs[sigs$log2FoldChange > 0.5, ])
genes_low <- rownames(sigs[sigs$log2FoldChange < -0.5, ])

# 🧬 GO-analyse voor verhoogde genen
GOresults_high <- enrichGO(
  gene = genes_high,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP"
)

# 🧬 GO-analyse voor verlaagde genen
GOresults_low <- enrichGO(
  gene = genes_low,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP"
)

# Zet de resultaten om naar data frames
GO_df_high <- as.data.frame(GOresults_high)
GO_df_low <- as.data.frame(GOresults_low)

# Sla de resultaten op als CSV-bestanden
write.csv2(GO_df_high, file = "GO_enrichment_results_high.csv", row.names = FALSE)
write.csv2(GO_df_low, file = "GO_enrichment_results_low.csv", row.names = FALSE)

# 🔝 Top 10 GO-termen per groep
top10_high <- GO_df_high %>%
  arrange(p.adjust) %>%
  slice_head(n = 5)

top10_low <- GO_df_low %>%
  arrange(p.adjust) %>%
  slice_head(n = 5)

# 📊 Visualisatie
dotplot_high <- dotplot(GOresults_high, showCategory = 5) +
  theme(axis.text.y = element_text(size = 8)) +  
  scale_fill_gradientn(colors = c("mediumorchid4", "#CD1076","#ae017e", "#c51b8a","#dd3497","#FF1493","#fa9fb5", "#fbb4d9"), name = "p.adjust" ) +
  labs( title = "top 5 GO-termen verhoogd🔬", y = "GO-termen") 


dotplot_low <- dotplot(GOresults_low, showCategory = 5) + 
  theme(axis.text.y = element_text(size = 8)) + 
  scale_fill_gradientn(colors = c("mediumorchid4", "#CD1076","#ae017e", "#c51b8a","#dd3497","#FF1493","#fa9fb5", "#fbb4d9"), name = "p.adjust" ) +
  labs( title = "top 5 GO-termen verlaagd🔬", y = "GO-termen") 

dotplotcombined <- dotplot_high/dotplot_low
dotplotcombined

barplot_low <- barplot(GOresults_low, showCategory = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5), axis.text.y = element_text(size = 10)) + 
  scale_fill_gradientn(colors = c("mediumorchid4", "#CD1076","#ae017e", "#c51b8a","#dd3497","#FF1493","#fa9fb5", "#fbb4d9"), name = "p.adjust" ) + 
  labs( title = "top 5 GO-termen verlaagd🔬", y = "GO-termen") 

barplot_high <- barplot(GOresults_high, showCategory = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5), axis.text.y = element_text(size = 10)) + 
  scale_fill_gradientn(colors = c("#68228B","#ae017e", "#c51b8a","#dd3497","#FF1493","#fa9fb5", "#fbb4d9"), name = "p.adjust" ) + 
  labs( title = "top 5 GO-termen verhoogd🔬", y = "GO-termen") 

combinedplot <- barplot_high/barplot_low
combinedplot
