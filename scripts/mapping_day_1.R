setwd("C:/Users/dejon/OneDrive - NHL Stenden/Transcriptomics J2P4 CASUS reuma/Data_RA_raw")
getwd()
# Vervang de bestandsnaam hieronder met je eigen zip-bestand

unzip("C:/Users/dejon/OneDrive - NHL Stenden/Transcriptomics J2P4", exdir = "ethanol_data") #Hiermee worden de bestanden uitgepakt in een subm"C:/Users/dejon/OneDrive - NHL Stenden/Transcriptomics J2P4/transfer_3140477_files_2d2312d2.zip"ap 'ethanol_data'

BiocManager::install('Rsubread')
library(Rsubread)

unzip("C:/Users/dejon/OneDrive - NHL Stenden/Transcriptomics J2P4/ncbi_dataset unzip.zip/ncbi_dataset/data/GCF_000005845.2", exdir = "refgene")



buildindex(
  basename = 'ref_Human',
  reference = "Homo_sapiens.GRCh38.dna.toplevel.fa",
  memory = 4000,
  indexSplit = TRUE)

# Ethanol monsters
align.control1 <- align(index = "ref_human", readfile1 = "SRR4785819_1_subset40k.fastq", readfile2 = "SRR4785819_2_subset40k.fastq", output_file = "control1.BAM")
align.control2 <- align(index = "ref_human", readfile1 = "SRR4785820_1_subset40k.fastq", readfile2 = "SRR4785820_2_subset40k.fastq", output_file = "control2.BAM")
align.control3 <- align(index = "ref_human", readfile1 = "SRR4785828_1_subset40k.fastq", readfile2 = "SRR4785828_2_subset40k.fastq", output_file = "control3.BAM")
align.control4 <- align(index = "ref_human", readfile1 = "SRR4785831_1_subset40k.fastq", readfile2 = "SRR4785831_2_subset40k.fastq", output_file = "control4.BAM")
align.reuma1 <- align(index = "ref_human", readfile1 = "SRR4785979_1_subset40k.fastq" , readfile2 = "SRR4785979_2_subset40k.fastq" , output_file = "reuma1.BAM")
align.reuma2 <- align(index = "ref_human", readfile1 = "SRR4785980_1_subset40k.fastq" , readfile2 = "SRR4785980_2_subset40k.fastq" , output_file = "reuma2.BAM")
align.reuma3 <- align(index = "ref_human", readfile1 = "SRR4785986_1_subset40k.fastq" , readfile2 = "SRR4785986_2_subset40k.fastq" , output_file = "reuma3.BAM")
align.reuma4 <- align(index = "ref_human", readfile1 = "SRR4785988_1_subset40k.fastq" , readfile2 = "SRR4785988_2_subset40k.fastq" , output_file = "reuma4.BAM")

# Laad Rsamtools voor sorteren en indexeren
library(Rsamtools)

# Bestandsnamen van de monsters
samples <- c('control1', 'control2', 'control3', 'control4', 'reuma1', 'reuma2', 'reuma3', 'reuma4')

# Voor elk monster: sorteer en indexeer de BAM-file
# Sorteer BAM-bestanden
lapply(samples, function(s) {sortBam(file = paste0(s, '.BAM'), destination = paste0(s, '.sorted'))
})

#INDEXEER 
lapply(samples, function(s) {
  indexBam(paste0(s, ".sorted.bam"))
})


fna_file <- "C:/Users/dejon/OneDrive - NHL Stenden/Transcriptomics J2P4 CASUS reuma/Data_RA_raw/Homo_sapiens.GRCh38.dna.toplevel.fa"

# Maak een .fai indexbestand aan

indexFa(fna_file)




