setwd("C:/Users/indys/OneDrive - NHL Stenden/J2P4/Transcriptomics/Casus")
getwd()
unzip("Ruwe Data/Data_RA_raw.zip", exdir = "RA_data")
library(BiocManager)
library(Rsubread)
buildindex(
  basename = 'ref_human',
  reference = '',
  memory = 4000,
  indexSplit = TRUE)