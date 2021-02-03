args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]

all_blast <- read.delim(paste0(filename, ".txt"), header = F)
all_blast_colnames <- c("molecule_reads", "refseq_ids", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
colnames(all_blast)[1:12] <- all_blast_colnames
save(all_blast, file = paste0(filename, ".RData"))
