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

