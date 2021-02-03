
rm(list = ls())

all_dir <- list.dirs(path = "data/split")
all_dir <- all_dir[-1] ##since first one its own directory

all_blast <- c()
missing_file <- c()
for(d in 1:length(all_dir)) {
  
  all_blst_files <- list.files(path = all_dir[d], pattern = "*.blst", full.names = T)
  
  for(i in 1:length(all_blst_files)) {
    
    if(file.info(all_blst_files[i])$size > 0) {
      temp_blst <- read.delim(all_blst_files[i], header = F)
      all_blast <- rbind(all_blast, temp_blst)
    } else 
      missing_file <- c(missing_file, all_blst_files[i])
    
  }
  
  print(d)
}

all_blast_colnames <- c("molecule_reads", "refseq_ids", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
colnames(all_blast)[1:12] <- all_blast_colnames
save(all_blast, file = "data/all_blast.RData")

#test_data <- read.delim("data/split/all.txt", header = F)
