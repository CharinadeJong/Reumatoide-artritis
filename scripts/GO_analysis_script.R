
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
ggplot(top10_high, aes(x = reorder(Description, p.adjust), y = -log10(p.adjust))) +
  geom_col(fill = "steelblue") +
  labs(title = "Top 10 GO-termen (verhoogde expressie)", x = "GO-term", y = "-log10(p.adjust)")

ggplot(top10_low, aes(x = reorder(Description, p.adjust), y = -log10(p.adjust))) +
  geom_col(fill = "firebrick") +
  labs(title = "Top 10 GO-termen (verlaagde expressie)", x = "GO-term", y = "-log10(p.adjust)")
