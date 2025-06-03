setwd("C:/Users/indys/OneDrive - NHL Stenden/J2P4 innovatievediagnostiek//Transcriptomics/Casus/Transcriptomics-Indy-J2P4/")
getwd()
unzip("Ruwe Data/Data_RA_raw.zip", exdir = "RA_data")
library(BiocManager)
library(Rsubread)
# bestand te groot dus onderstaande code meegeven in verslag op een of andere manier. 
buildindex(
  basename = 'ref_human',
  reference = 'Human genome/ncbi_dataset (1)/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 4000,
  indexSplit = TRUE)
#controle groepen te vinden op bb
align.ctrl19 <- align(index = "C:/Users/indys/OneDrive - NHL Stenden/J2P4 innovatievediagnostiek/Transcriptomics/Casus/Human genome/ref_human", 
                                      readfile1 = "Ruwe Data/Data_RA_raw/SRR4785819_1_subset40k.fastq", 
                                        readfile2 = "Ruwe Data/Data_RA_raw/SRR4785819_2_subset40k.fastq",
                                        output_file = "Ctrl19.BAM")

align.ctrl20 <- align(index = "C:/Users/indys/OneDrive - NHL Stenden/J2P4 innovatievediagnostiek/Transcriptomics/Casus/Human genome/ref_human", 
                                          readfile1 = "Ruwe Data/Data_RA_raw/SRR4785820_1_subset40k.fastq", 
                                           readfile2 = "Ruwe Data/Data_RA_raw/SRR4785820_2_subset40k.fastq",
                                            output_file = "Ctrl20.BAM")

align.ctrl28 <- align(index = "C:/Users/indys/OneDrive - NHL Stenden/J2P4 innovatievediagnostiek/Transcriptomics/Casus/Human genome/ref_human", 
                                          readfile1 = "Ruwe Data/Data_RA_raw/SRR4785828_1_subset40k.fastq", 
                                          readfile2 = "Ruwe Data/Data_RA_raw/SRR4785828_2_subset40k.fastq",
                                          output_file = "Ctrl28.BAM")

align.ctrl31 <- align(index = "C:/Users/indys/OneDrive - NHL Stenden/J2P4 innovatievediagnostiek/Transcriptomics/Casus/Human genome/ref_human", 
                      readfile1 = "Ruwe Data/Data_RA_raw/SRR4785831_1_subset40k.fastq", 
                      readfile2 = "Ruwe Data/Data_RA_raw/SRR4785831_2_subset40k.fastq",
                      output_file = "Ctrl31.BAM")
#test groepen RA
align.RA79 <- align(index = "C:/Users/indys/OneDrive - NHL Stenden/J2P4 innovatievediagnostiek/Transcriptomics/Casus/Human genome/ref_human", 
                      readfile1 = "Ruwe Data/Data_RA_raw/SRR4785979_1_subset40k.fastq", 
                      readfile2 = "Ruwe Data/Data_RA_raw/SRR4785979_2_subset40k.fastq",
                      output_file = "RA79.BAM")

align.RA80 <- align(index = "C:/Users/indys/OneDrive - NHL Stenden/J2P4 innovatievediagnostiek/Transcriptomics/Casus/Human genome/ref_human", 
                  readfile1 = "Ruwe Data/Data_RA_raw/SRR4785980_1_subset40k.fastq", 
                  readfile2 = "Ruwe Data/Data_RA_raw/SRR4785980_2_subset40k.fastq",
                  output_file = "RA80.BAM")

align.RA86 <- align(index = "C:/Users/indys/OneDrive - NHL Stenden/J2P4 innovatievediagnostiek/Transcriptomics/Casus/Human genome/ref_human", 
                  readfile1 = "Ruwe Data/Data_RA_raw/SRR4785986_1_subset40k.fastq", 
                  readfile2 = "Ruwe Data/Data_RA_raw/SRR4785986_2_subset40k.fastq",
                  output_file = "RA86.BAM")

align.RA88 <- align(index = "C:/Users/indys/OneDrive - NHL Stenden/J2P4 innovatievediagnostiek/Transcriptomics/Casus/Human genome/ref_human", 
                  readfile1 = "Ruwe Data/Data_RA_raw/SRR4785988_1_subset40k.fastq", 
                  readfile2 = "Ruwe Data/Data_RA_raw/SRR4785988_2_subset40k.fastq",
                  output_file = "RA88.BAM")

# Laad Rsamtools voor sorteren en indexeren
library(Rsamtools)
samples <- c('ctrl19', 'ctrl20', 'ctrl28', 'ctrl31', 'RA79', 'RA80', 'RA86', 'RA88')

# Voor elk monster: sorteer en indexeer de BAM-file
# Sorteer BAM-bestanden
lapply(samples, function(s) {sortBam(file = paste0(s, '.BAM'), destination = paste0(s, '.sorted'))
})

lapply(samples, function(s) {indexBam(file = paste0(s, '.sorted.bam'))
})


# voor werkcollege 2 ispaired end = TRUE en niet false 
#annot.ext is de gtF van humaan esamble
# Je definieert een vector met namen van BAM-bestanden. Elke BAM bevat reads van een RNA-seq-experiment (bijv. behandeld vs. controle).
library(readr)
library(dplyr)
library(Rsamtools)
library(Rsubread)
  allsamples <- c("BAM bestanden/ctrl19.BAM", "BAM bestanden/ctrl20.BAM", "BAM bestanden/ctrl28.BAM", "BAM bestanden/ctrl31.BAM", "BAM bestanden/RA79.BAM", "BAM bestanden/RA80.BAM", "BAM bestanden/RA86.BAM", "BAM bestanden/RA88.BAM")

count_matrix <- featureCounts(
  files = allsamples,
    annot.ext = "GTF Humaan/Homo_sapiens.GRCh38.114.chr.gtf.gz",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)

head(count_matrix$annotation)
head(count_matrix$counts)
# Bekijk eerst de structuur van het object
str(count_matrix)
# Haal alleen de matrix met tellingen eruit
counts <- count_matrix$counts
colnames(counts) <- c( "ctrl19.BAM", "ctrl20.BAM", "ctrl28.BAM", "ctrl31.BAM", "RA79.BAM", "RA80.BAM", "RA86.BAM", "RA88.BAM")
write.csv(counts, "bewerkt_countmatrix.csv")
#count matrix krijg de goede van dewi dus die toevoegen in mapje!

countstest <- read.table("Ruwe Data/count_matrix.txt")
treatment <- c("control", "control", "control","control", "RA", "RA", "RA", "RA")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- c('SRR4785819', 'SRR4785820', 'SRR4785828', 'SRR4785831', 'SRR4785979', 'SRR4785980', 'SRR4785986', 'SRR4785988')
library(DESeq2)
library(KEGGREST)
# Maak DESeqDataSet aan
dds <- DESeqDataSetFromMatrix(countData = round(countstest),
                              colData = treatment_table,
                              design = ~ treatment)
# Voer analyse uit
dds <- DESeq(dds)
resultaten <- results(dds)
# Resultaten opslaan in een bestand
#Bij het opslaan van je tabel kan je opnieuw je pad instellen met `setwd()` of het gehele pad waar je de tabel wilt opslaan opgeven in de code.

write.csv(resultaten, 'Ruwe Data/Resultaten/ResultatenRA.csv')
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)
#bovenste is hoeveel genen meer in expressie en de onderste hoeveel genen minder expressie hebben
hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]
head(laagste_p_waarde)

if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}
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
dev.copy(png, 'VolcanoplotRA.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()
if (!requireNamespace("goseq", quietly = TRUE)) {
  BiocManager::install("goseq")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
library(goseq)
library(org.Hs.eg.db)
library(DESeq2)




# TEST
# data
resultaten <- read.csv("Ruwe Data/Resultaten/ResultatenRA.csv", row.names = 1)
head(resultaten)

sig <- resultaten %>%
  filter(padj < 0.05)
head(sig)
sig <- rownames(sig)
sig

all <- rownames(resultaten)


gene.vector <- as.integer(all%in%sig)
names(gene.vector) <- all


pwf <- nullp(gene.vector, "hg19", "geneSymbol")
GO.wall = goseq(pwf, "hg19", "geneSymbol")
enriched.GO=GO.wall$category[GO.wall$over_represented_pvalue<.05]
#NOTE: They recommend using a more stringent multiple testing corrected p value here
capture.output(for(go in enriched.GO[1:258]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
, file="SigGo.txt")

GO.results = goseq(pwf, "hg19", "geneSymbol")
GO.results %>%
  top_n(10, wt = -over_represented_pvalue) %>%
  mutate(hitsPerc = numDEInCat*100/numInCat) %>%
  ggplot(aes(x=hitsPerc, 
             y = term,
             colour = over_represented_pvalue,
             size = numDEInCat)) +
  geom_point() +
  expand_limits(x =0)+
  labs(x = "hits(%)", y = "GO term", colour = "p value", size = "count")
  
library(GO.db)

GOTERM [GO.results$category[1]]

#KEGG PATHWAY ANALYSE
if (!requireNamespace("pathview", quietly = TRUE)) {
  BiocManager::install("pathview")
}
library(pathview)
resultaten[1] <- NULL
resultaten[2:5] <- NULL

pathview(
  gene.data = gene_vector,
  pathway.id = "03260",  # KEGG ID voor virion - human immunodeficiency virus 
  species = "hsa",          # 'hgn' = human genome 
  gene.idtype = "ENTREZ",     # Geef aan dat het KEGG-ID's zijn
  limit = list(gene = 5)    # Kleurbereik voor log2FC van -5 tot +5
)

