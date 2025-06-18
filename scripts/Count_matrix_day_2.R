library(readr)
library(dplyr)
library(Rsamtools)
library(Rsubread)



#GFF3 inlezen en voorbereiden
#Let op! Deze stap hoef je niet te herhalen in de casus.

# Inlezen en filteren van GFF3-bestand
gff <- read_tsv("", comment = "#", col_names = FALSE)
# Kolomnamen toevoegen
colnames(gff) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
# Alleen genregels selecteren
gff_gene <- gff %>% filter(type == "gene")
# 'type' aanpassen naar 'exon' zodat featureCounts het accepteert
gff_gene$type <- "exon"

# Extraheer de chromosoomnaam uit BAM-header
bam_chr <- names(scanBamHeader("C:/Users/dejon/OneDrive - NHL Stenden/Transcriptomics J2P4/ncbi_dataset/ncbi_dataset/data/GCF_000005845.2/eth1.sorted.bam")[[1]]$targets)[1]
gff_gene$seqid <- bam_chr

write_delim(gff_gene, "ecoli_ready_for_featureCounts.gtf", delim = "\t", col_names = FALSE)

#featureCounts uitvoeren
#Let op! Deze stap moet je wel uitvoeren in de casus.
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

#Count matrix omzetten naar count matrix gekregen van Dewi
count_matrix <- read.table("count_matrix.txt") 

head(count_matrix$annotation)

head(count_matrix$counts)

# Bekijk eerst de structuur van het object
str(count_matrix)

# Haal alleen de matrix met tellingen eruit, niet voor casus

counts <- count_matrix$counts


head(counts)


# Stap 2: Stel duidelijke kolomnamen in
colnames(count_matrix) <- c ( "control1", "control2", "control3", "control4", "reuma1", "reuma2", "reuma3", "reuma4")

#De eerste kolom bevat de gen-ID’s. Je kunt deze beter als rij-namen gebruiken (in plaats van als gewone kolom).
#rownames(counts) <- counts[, 1]
#colnames(counts) <- c("GeneID", "reuma1.BAM", "reuma2.BAM", "reuma3.BAM", "reuma4.BAM", "control1.BAM", "control2.BAM", "control3.BAM", "control4.BAM")

# 3. Verwijder de eerste kolom (GeneID), want die zit nu in de rij-namen
#counts <- counts[, -1]

# Stap 5: Sla de bewerkte matrix op als CSV-bestand
write.csv(count_matrix, "bewerkt_countmatrix.csv")


# Stap 6: Bekijk de eerste paar rijen om te controleren
head(count_matrix)

