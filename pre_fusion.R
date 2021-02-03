
rm(list = ls())

library(parallel)

#start-------------------merge refseq and ucsc refseq to gene symbols---------------------
# refseq2gene <- read.delim(file = "data/all_refseq_gene.txt", header = T)
# ucsc_refseq_genes <- read.delim(file = "data/all_refseq_gene_ucsc.txt", header = T)
# 
# refseq2gene <- refseq2gene[-which(duplicated(refseq2gene[ ,"name"])), ]
# refseq2gene <- as.matrix(refseq2gene)
# 
# refseq2gene[ ,"name"] <- unlist(lapply(strsplit(refseq2gene[ ,"name"], split = "[.]"), function(x) {return(x[1])}))
# 
# ucsc_refseq_genes <- ucsc_refseq_genes[ucsc_refseq_genes[ ,"name"] %in% setdiff(ucsc_refseq_genes[ ,"name"], refseq2gene[ ,"name"]), ]
# ucsc_refseq_genes <- ucsc_refseq_genes[-which(duplicated(ucsc_refseq_genes[ ,"name"])), ]
# ucsc_refseq_genes <- as.matrix(ucsc_refseq_genes)
# 
# combined_refseq2gene <- rbind(refseq2gene, ucsc_refseq_genes)

# save(combined_refseq2gene, file = "data/combined_refseq2gene.RData")
#end-------------------merge refseq and ucsc refseq to gene symbols---------------------

#start-------Update blast files with gene_symbol, and chr-----------------------------------
Update_all_blast <- function(all_blast, combined_refseq2gene) {
  
  empty_mat <- matrix("", nrow = nrow(all_blast), ncol = 5)
  colnames(empty_mat) <- c("gene_symbol", "chr", "strand", "startPos", "endPos")
  
  all_blast <- as.matrix(all_blast)
  all_blast <- cbind(all_blast, empty_mat)
  all_blast[ ,"refseq_ids"] <- unlist(lapply(strsplit(all_blast[ ,"refseq_ids"], split = "[.]"), function(x) {return(x[1])}))
  
  matched_index <- match(all_blast[ ,"refseq_ids"], combined_refseq2gene[ ,"name"])
  
  all_blast[ ,c("gene_symbol", "chr", "strand", "startPos", "endPos")] <- combined_refseq2gene[matched_index, c("name2", "chrom", "strand", "txStart", "txEnd")]
  
  return(all_blast)
  
}

load("data/combined_refseq2gene.RData")
data_path <-"data/data_control/"
blast_filename <- "P15059_sample5_counting_all_blast"

load(paste0(data_path, blast_filename, ".RData"))
all_blast <- Update_all_blast(all_blast, combined_refseq2gene)

output_dir <- paste0(data_path, blast_filename)
dir.create(output_dir)
#------------Remove RefseqIDs didn't find any Gene Symbol------------------
write_info <- file(paste0(output_dir, "/sample_info.txt"), open = "w")

writeLines("Refseq IDs matched with GeneSymbol: ", write_info, sep = "\n")
writeLines(paste0("Total Refseq IDs matched with GeneSymbol: ", length(all_blast[ ,"refseq_ids"])), write_info, sep = "\n")
writeLines(paste0("Total Unique Refseq IDs: ", length(unique(all_blast[ ,"refseq_ids"]))), write_info, sep = "\n")
writeLines(paste0("Total missed Refseq IDs (not matched with GeneSymbol): ", length(all_blast[is.na(all_blast[ ,"gene_symbol"]), "refseq_ids"])), write_info, sep = "\n")
writeLines(paste0("Total missed Unique Refseq IDs: ", length(unique(all_blast[is.na(all_blast[ ,"gene_symbol"]), "refseq_ids"]))), write_info, sep = "\n")

all_blast <- all_blast[!is.na(all_blast[ ,"gene_symbol"]), ]

close(write_info)
#---------------Filter quality blast results only----------------------

pdf(file = paste0(output_dir, "/bitscore_hist.pdf"))
hist(as.numeric(all_blast[ ,"bitscore"]), xlab = "bitscore", main = "")
#lines(density(as.numeric(all_blast[ ,"bitscore"])), lwd = 2, col = "chocolate3")
abline(v = mean(as.numeric(all_blast[ ,"bitscore"])), col = "royalblue", lwd = 2)
abline(v = median(as.numeric(all_blast[ ,"bitscore"])), col = "red", lwd = 2)
legend(x = "topleft", c("Mean", "Median"), col = c("royalblue", "red"), lwd = c(2, 2, 2))
dev.off()
param_bitscore <- 200
if(param_bitscore >= median(as.numeric(all_blast[ ,"bitscore"]))) param_bitscore <- median(as.numeric(all_blast[ ,"bitscore"]))
# test_mat <- all_blast[as.numeric(all_blast[ ,"mismatch"]) + as.numeric(all_blast[ ,"gapopen"]) <= 1 & as.numeric(all_blast[ ,"evalue"]) <= 1e-50 & as.numeric(all_blast[ ,"bitscore"]) >= 200, ]
all_blast <- all_blast[as.numeric(all_blast[ ,"pident"]) >= 100 & as.numeric(all_blast[ ,"evalue"]) <= 1e-50 & as.numeric(all_blast[ ,"bitscore"]) >= param_bitscore, ]


#----------------add molecules column-----------
all_blast <- cbind(substr(all_blast[ ,"molecule_reads"], 1, 12), all_blast)
colnames(all_blast)[1] <- "molecules"
save(all_blast, file = paste0(output_dir, "/all_blast.RData"))

rm(combined_refseq2gene)
gc()

#--------------------------
print("list of molecules per gene calculation ......")
all_genes_symbols <- unique(all_blast[ ,"gene_symbol"])
all_genes_mols_list <- vector(mode = "list", length = length(all_genes_symbols))
names(all_genes_mols_list) <- all_genes_symbols
print(paste0("Total genes: ", length(all_genes_symbols)))
# all_genes_mols_list <- lapply(all_genes_symbols, function(x) {return(all_genes_mols_list[x] <- unique(all_blast[all_blast[ ,"gene_symbol"] == x, "molecules"]))})
all_genes_mols_list <- mclapply(all_genes_symbols, function(x) {return(all_genes_mols_list[x] <- unique(all_blast[all_blast[ ,"gene_symbol"] == x, "molecules"]))}, mc.cores = 10)

print("gene gene molecules sharing calculation......")
gene_gene_mol_sharing_mat <- matrix(NA, nrow = length(all_genes_symbols), ncol = length(all_genes_symbols))
colnames(gene_gene_mol_sharing_mat) <- all_genes_symbols
rownames(gene_gene_mol_sharing_mat) <- all_genes_symbols

total_genes <- length(all_genes_mols_list)
total_work <- ((total_genes*total_genes)-total_genes)/2
per_done <- unlist(lapply(seq(10, 100, by=10), function(x) {return(round((total_work*x)/100))}))
names(per_done) <- seq(10, 100, by=10)
it_counter <- 0

for(i in 1:(length(all_genes_mols_list)-1)) {
  
  for(j in (i+1):length(all_genes_mols_list)) {
    
    gene_gene_mol_sharing_mat[i, j] <- length(intersect(all_genes_mols_list[[i]], all_genes_mols_list[[j]]))
    it_counter <- it_counter + 1
    if(it_counter %in% per_done) print(paste0(names(per_done[per_done==it_counter]), "% is done ..."))
  }
}

gene_gene_mol_sharing_mat[gene_gene_mol_sharing_mat==0] <- NA
temp_pos <- which(!is.na(gene_gene_mol_sharing_mat), arr.ind = T)
gene_gene_mol_sharing_mat <- cbind(rownames(gene_gene_mol_sharing_mat)[temp_pos[ ,"row"]], colnames(gene_gene_mol_sharing_mat)[temp_pos[ ,"col"]], gene_gene_mol_sharing_mat[temp_pos])
colnames(gene_gene_mol_sharing_mat) <- c("gene1", "gene2", "#shared_mole")
gene_gene_mol_sharing_mat <- gene_gene_mol_sharing_mat[order(as.numeric(gene_gene_mol_sharing_mat[ ,"#shared_mole"]), decreasing = T), ]

pdf(file = paste0(output_dir, "/gene_gene_mol_sharing_hist.pdf"))
hist(as.numeric(gene_gene_mol_sharing_mat[ ,"#shared_mole"]), main = "", xlab = "#gene-gene shared molecule")
dev.off()

save(gene_gene_mol_sharing_mat, file = paste0(output_dir, "/gene_gene_mol_sharing_mat.RData"))
print("Success End......")

