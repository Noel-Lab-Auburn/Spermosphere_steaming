#run the below two lines before running the script
#source /opt/asn/etc/asn-bash-profiles-special/modules.sh
#module load R/4.1.0

library(dada2)
library(Biostrings)

taxonomy.file.path <- "~/noel_shared/db_fungi/sh_general_release_dynamic_s_all_25.07.2023_mockseqadded.fasta"
otus.file.path <- "clustered/otus.fasta"

# Fasta
FASTA.otus <- readDNAStringSet(otus.file.path, format="fasta", seek.first.rec=TRUE, use.names=TRUE)

taxa <- assignTaxonomy(FASTA.otus, taxonomy.file.path, multithread=TRUE, tryRC = TRUE)

write.csv(taxa, file = "NBCtaxa.csv")