
# Laad benodigde packages
library(readr)
library(dplyr)
library(Rsamtools)
library(Rsubread)

BiocManager::install("Rsamtools")
install.packages("Rsamtools")

setwd("C:/Users/dejon/OneDrive - NHL Stenden/Transcriptomics J2P4 CASUS reuma/Data_RA_raw")

# Je definieert een vector met namen van BAM-bestanden. Elke BAM bevat reads van een RNA-seq-experiment (bijv. behandeld vs. controle).

allsamples <- c ("control1.BAM", "control2.BAM", "control3.BAM", "control4.BAM", "reuma1.BAM", "reuma2.BAM", "reuma3.BAM", "reuma4.BAM")

count_matrix <- featureCounts(
  files = allsamples,
  annot.ext = "Homo_sapiens.GRCh38.114.chr.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)

#Count matrix omzetten naar volledige count matrix (afkomstig van school) 
count_matrix <- read.table("count_matrix.txt") 


# Bekijk eerst de structuur van het object
str(count_matrix)


# Stap 2: Stel duidelijke kolomnamen in
colnames(count_matrix) <- c ( "control1", "control2", "control3", "control4", "reuma1", "reuma2", "reuma3", "reuma4")


# Stap 5: Sla de bewerkte matrix op als CSV-bestand
write.csv(count_matrix, "bewerkt_countmatrix.csv")


# Stap 6: Bekijk de eerste paar rijen om te controleren
head(count_matrix)

