
rm(list = ls())

Edge_naming <- function(conn_mat_p) {
  
  conn_mat_p <- cbind(conn_mat_p, paste(conn_mat_p[ ,"gene1"], conn_mat_p[ ,"gene2"], sep = "-vs-"))
  colnames(conn_mat_p)[4] <- "edges"
  rownames(conn_mat_p) <- paste(conn_mat_p[ ,"gene2"], conn_mat_p[ ,"gene1"], sep = "-vs-")
  
  return(conn_mat_p)
}

Common_edges <- function(l_edges_p, l_rnames_p) {
  
  t_common_edges <- Reduce(intersect, l_edges_p)
  t_common_edges <- strsplit(t_common_edges, split = "-vs-")
  t_common_edges <- unlist(lapply(t_common_edges, function(x) {return(paste(x[2], x[1], sep = "-vs-"))}))
  
  common_edges <- Reduce(intersect, l_rnames_p) 
  common_edges <- union(common_edges, t_common_edges)
  
  return(common_edges)
}

#---------------Prepare each tumor and control data---------------
gg_conn_t1 <- get(load("data/data_tumor/P14754_CGCTCATT_sample1_counting_all_blast/gene_gene_mol_sharing_mat.RData"))
gg_conn_t2 <- get(load("data/data_tumor/P14754_CGCTCATT_sample2_counting_all_blast/gene_gene_mol_sharing_mat.RData"))
gg_conn_t3 <- get(load("data/data_tumor/P14754_CGCTCATT_sample3_counting_all_blast/gene_gene_mol_sharing_mat.RData"))
gg_conn_t4 <- get(load("data/data_tumor/P14754_CGCTCATT_sample4_counting_all_blast/gene_gene_mol_sharing_mat.RData"))
gg_conn_t5 <- get(load("data/data_tumor/P14754_CGCTCATT_sample5_counting_all_blast/gene_gene_mol_sharing_mat.RData"))
gg_conn_t6 <- get(load("data/data_tumor/P14754_CGCTCATT_sample6_counting_all_blast/gene_gene_mol_sharing_mat.RData"))


#gg_conn_c5 <- get(load("data/data_control/P15059_sample5_counting_all_blast/gene_gene_mol_sharing_mat.RData"))
gg_conn_c6 <- get(load("data/data_control/P15059_sample6_counting_all_blast/gene_gene_mol_sharing_mat.RData"))
gg_conn_c8 <- get(load("data/data_control/P15059_sample8_counting_all_blast/gene_gene_mol_sharing_mat.RData"))


gg_conn_t1 <- Edge_naming(gg_conn_t1)
gg_conn_t2 <- Edge_naming(gg_conn_t2)
gg_conn_t3 <- Edge_naming(gg_conn_t3)
gg_conn_t4 <- Edge_naming(gg_conn_t4)
gg_conn_t5 <- Edge_naming(gg_conn_t5)
gg_conn_t6 <- Edge_naming(gg_conn_t6)


l_edges <- list(gg_conn_t1[ ,"edges"], gg_conn_t2[ ,"edges"])
l_rnames <- list(rownames(gg_conn_t1), rownames(gg_conn_t2))
common_edges <- Common_edges(l_edges, l_rnames)
merged_tumor1_conn_mat <- cbind(gg_conn_t1[common_edges, ], gg_conn_t2[common_edges, "#shared_mole"])

l_edges <- list(gg_conn_t3[ ,"edges"], gg_conn_t4[ ,"edges"])
l_rnames <- list(rownames(gg_conn_t3), rownames(gg_conn_t4))
common_edges <- Common_edges(l_edges, l_rnames)
merged_tumor2_conn_mat <- cbind(gg_conn_t3[common_edges, ], gg_conn_t4[common_edges, "#shared_mole"])

l_edges <- list(gg_conn_t5[ ,"edges"], gg_conn_t6[ ,"edges"])
l_rnames <- list(rownames(gg_conn_t5), rownames(gg_conn_t6))
common_edges <- Common_edges(l_edges, l_rnames)
merged_tumor3_conn_mat <- cbind(gg_conn_t5[common_edges, ], gg_conn_t6[common_edges, "#shared_mole"])


#gg_conn_c5 <- Edge_naming(gg_conn_c5)
gg_conn_c6 <- Edge_naming(gg_conn_c6)
gg_conn_c8 <- Edge_naming(gg_conn_c8)

l_edges <- list(gg_conn_c6[ ,"edges"], gg_conn_c8[ ,"edges"])
l_rnames <- list(rownames(gg_conn_c6), rownames(gg_conn_c8))
common_edges <- Common_edges(l_edges, l_rnames)
merged_control_conn_mat <- cbind(gg_conn_c6[common_edges, ], gg_conn_c8[common_edges, "#shared_mole"])

#----------------calculate tumor conn only based on control set but not other tumor sets---------------------------
only_t1_conn_mat <- merged_tumor1_conn_mat[setdiff(rownames(merged_tumor1_conn_mat), rownames(merged_control_conn_mat)), ]
only_t1_conn_mat <- cbind(only_t1_conn_mat, "")
only_t1_conn_mat <- cbind(only_t1_conn_mat, "")
colnames(only_t1_conn_mat)[c(6, 7)] <- c("mean", "median")
temp_only_t_conn_mat <- apply(only_t1_conn_mat[ ,c(3,5)], c(1,2), as.numeric) 
only_t1_conn_mat[ ,"mean"] <- round(apply(temp_only_t_conn_mat, 1, mean))
only_t1_conn_mat[ ,"median"] <- round(apply(temp_only_t_conn_mat, 1, median))

only_t2_conn_mat <- merged_tumor2_conn_mat[setdiff(rownames(merged_tumor2_conn_mat), rownames(merged_control_conn_mat)), ]
only_t2_conn_mat <- cbind(only_t2_conn_mat, "")
only_t2_conn_mat <- cbind(only_t2_conn_mat, "")
colnames(only_t2_conn_mat)[c(6, 7)] <- c("mean", "median")
temp_only_t_conn_mat <- apply(only_t2_conn_mat[ ,c(3,5)], c(1,2), as.numeric) 
only_t2_conn_mat[ ,"mean"] <- round(apply(temp_only_t_conn_mat, 1, mean))
only_t2_conn_mat[ ,"median"] <- round(apply(temp_only_t_conn_mat, 1, median))

only_t3_conn_mat <- merged_tumor3_conn_mat[setdiff(rownames(merged_tumor3_conn_mat), rownames(merged_control_conn_mat)), ]
only_t3_conn_mat <- cbind(only_t3_conn_mat, "")
only_t3_conn_mat <- cbind(only_t3_conn_mat, "")
colnames(only_t3_conn_mat)[c(6, 7)] <- c("mean", "median")
temp_only_t_conn_mat <- apply(only_t3_conn_mat[ ,c(3,5)], c(1,2), as.numeric) 
only_t3_conn_mat[ ,"mean"] <- round(apply(temp_only_t_conn_mat, 1, mean))
only_t3_conn_mat[ ,"median"] <- round(apply(temp_only_t_conn_mat, 1, median))

rm(temp_only_t_conn_mat)

# save(only_t1_conn_mat, file = "data/t_analysis/only_t1_conn_mat.RData")
# save(only_t2_conn_mat, file = "data/t_analysis/only_t2_conn_mat.RData")
# save(only_t3_conn_mat, file = "data/t_analysis/only_t3_conn_mat.RData")


#--------------------Individual tumor analysis----------
rm(list = ls())

#--------Process_Genes_with_Multiple_Chr_Pos-------------
Process_Genes_with_Multiple_Chr_Location <- function(gene_symbols_with_diff_Chr_region, temp_mat_1) {
  
  for(i in 1:length(gene_symbols_with_diff_Chr_region)) {
    x <- gene_symbols_with_diff_Chr_region[i]
    t_start_pos <- min(temp_mat_1[temp_mat_1[ ,"Gene.name"] == x, "Gene.start..bp."])
    t_end_pos <- max(temp_mat_1[temp_mat_1[ ,"Gene.name"] == x, "Gene.end..bp."])
    temp_mat_1[min(which(temp_mat_1[ ,"Gene.name"] == x)), c("Gene.start..bp.", "Gene.end..bp.")] <- c(t_start_pos, t_end_pos)
    temp_mat_1 <- temp_mat_1[-setdiff(which(temp_mat_1[ ,"Gene.name"] == x), min(which(temp_mat_1[ ,"Gene.name"] == x))), ]
  }
  
  return(temp_mat_1)
}
#------end-----------

only_t_conn_mat <- get(load("data/t_analysis/only_t3_conn_mat.RData"))
only_t_conn_mat <- only_t_conn_mat[ ,c(1,2,4,3,5,6,7)]
colnames(only_t_conn_mat) <- c("Gene1", "Gene2", "Edge", "#Shared_mole_sample1", "#Shared_mole_sample2", "Mean", "Median")

#--------------measure gene-gene distance---------------------

all_f_genes <- union(only_t_conn_mat[ ,"Gene1"], only_t_conn_mat[ ,"Gene2"])

mart_gene_chr_pos <- read.delim2(file = "data/mart_export_hg38_unique.txt", header = T, sep = "\t")
mart_gene_chr_pos <- mart_gene_chr_pos[ ,c("Gene.name", "Chromosome.scaffold.name", "Gene.start..bp.", "Gene.end..bp.", "Gene.type")]
mart_gene_chr_pos <-  mart_gene_chr_pos[mart_gene_chr_pos[ ,"Gene.name"] %in% all_f_genes, ]

gene_symbols_with_diff_Chr_region <- as.character(mart_gene_chr_pos[ ,"Gene.name"][duplicated(mart_gene_chr_pos[ ,"Gene.name"])])
mart_gene_chr_pos <- Process_Genes_with_Multiple_Chr_Location(gene_symbols_with_diff_Chr_region, mart_gene_chr_pos)

mart_gene_chr_pos <- cbind(mart_gene_chr_pos, c(abs(mart_gene_chr_pos[ ,"Gene.end..bp."] - mart_gene_chr_pos[ ,"Gene.start..bp."])))
colnames(mart_gene_chr_pos) <- c("GeneSymbol","Chr","startPos", "endPos", "GeneType", "GeneLength")
mart_gene_chr_pos <- as.matrix(mart_gene_chr_pos)

missing_gene_symbols <- setdiff(all_f_genes, as.character(mart_gene_chr_pos[ ,"GeneSymbol"]))

ncbi_gene_info <- read.delim2(file = "data/Homo_sapiens.gene_info", header = T, sep = "\t")
ncbi_gene_info <- ncbi_gene_info[ncbi_gene_info[ ,"Symbol"] %in% missing_gene_symbols, c("Symbol", "chromosome", "type_of_gene")]
ncbi_gene_info <- cbind(ncbi_gene_info[ ,c(1:2)], "", "", ncbi_gene_info[ ,3], "")
colnames(ncbi_gene_info) <- c("GeneSymbol","Chr","startPos", "endPos", "GeneType", "GeneLength")
ncbi_gene_info <- as.matrix(ncbi_gene_info)

genes_info_mat <- rbind(mart_gene_chr_pos, ncbi_gene_info)
genes_info_mat[ ,"Chr"] <- paste0("chr", genes_info_mat[ ,"Chr"])

missing_gene_symbols <- setdiff(all_f_genes, as.character(genes_info_mat[ ,"GeneSymbol"]))

load("data/combined_refseq2gene.RData")
combined_refseq2gene <- combined_refseq2gene[combined_refseq2gene[ ,"name2"] %in% missing_gene_symbols, ]
combined_refseq2gene <- combined_refseq2gene[ ,c("chrom", "name2")]
combined_refseq2gene <- unique.matrix(combined_refseq2gene)
combined_refseq2gene <- combined_refseq2gene[ ,c(2,1)]
combined_refseq2gene <- cbind(combined_refseq2gene, "", "", "", "")
colnames(combined_refseq2gene) <- c("GeneSymbol","Chr","startPos", "endPos", "GeneType", "GeneLength")

genes_info_mat <- rbind(genes_info_mat, combined_refseq2gene)
rownames(genes_info_mat) <- genes_info_mat[ ,"GeneSymbol"]

#--------------Set gene-gene ConnType-----------------------------
only_t_conn_mat <- cbind(only_t_conn_mat, "")
colnames(only_t_conn_mat)[ncol(only_t_conn_mat)] <- "ConnType"
for(i in 1:nrow(only_t_conn_mat)) {
  
  t_res <- unique(genes_info_mat[genes_info_mat[ ,"GeneSymbol"] %in% only_t_conn_mat[i, c("Gene1", "Gene2")], "Chr"])
  if(length(t_res) == 1) only_t_conn_mat[i, "ConnType"] <- "cis" else only_t_conn_mat[i, "ConnType"] <- "trans"
}

#---------Set gene-gene Dist---------------------
only_t_conn_mat <- cbind(only_t_conn_mat, "")
colnames(only_t_conn_mat)[ncol(only_t_conn_mat)] <- "Dist"

for(i in 1:nrow(only_t_conn_mat)) {
  
  if(only_t_conn_mat[i, "ConnType"] == "cis") {
    g1_pos <- as.numeric(genes_info_mat[only_t_conn_mat[i, "Gene1"], c("startPos", "endPos")])
    g2_pos <- as.numeric(genes_info_mat[only_t_conn_mat[i, "Gene2"], c("startPos", "endPos")])
    if(is.na(g1_pos) | is.na(g2_pos)) {
      only_t_conn_mat[i, "Dist"] <- "missing"
    } else {
      if(length(intersect(seq(g1_pos[1], g1_pos[2]), seq(g2_pos[1], g2_pos[2]))) > 0) {
        only_t_conn_mat[i, "Dist"] <- "overlapped"
      } else {
        if(g1_pos[2] < g2_pos[1]) only_t_conn_mat[i, "Dist"] <- g2_pos[1] - g1_pos[2]
        else only_t_conn_mat[i, "Dist"] <- g1_pos[1] - g2_pos[2]
      }
    }
    
  }
  
}

#----------Filter based on gene Type-----------
rm_genes <- genes_info_mat[genes_info_mat[ ,"GeneType"] %in% c("rRNA"), "GeneSymbol"]
only_t_conn_mat <- only_t_conn_mat[!only_t_conn_mat[ ,"Gene1"] %in% rm_genes & !only_t_conn_mat[ ,"Gene2"] %in% rm_genes, ]

#----------Filte based on gene length------------
rm_genes <- setdiff(genes_info_mat[as.numeric(genes_info_mat[ ,"GeneLength"]) >= 1e+06, "GeneSymbol"], rm_genes)
rm_genes <- rm_genes[!is.na(rm_genes)] #---for some genes we dont have any length/size, theose produce NA
only_t_conn_mat <- only_t_conn_mat[!only_t_conn_mat[ ,"Gene1"] %in% rm_genes & !only_t_conn_mat[ ,"Gene2"] %in% rm_genes, ]

# save(only_t_conn_mat, file = "data/t_analysis/t3_res/only_t3_conn_mat.RData")
# save(genes_info_mat, file = "data/t_analysis/t3_res/t3_genes_info_mat.RData")

#end--------------measure gene-gene distance ()---------------------


#----------------Old Visualization of analysis------------------------

rm(list = ls())
only_t_conn_mat <- get(load("data/t_analysis/t3_res/only_t3_conn_mat.RData"))
genes_info_mat <- get(load("data/t_analysis/t3_res/t3_genes_info_mat.RData"))
fig_path <- "data/t_analysis/t3_res/"

pdf(paste0(fig_path, "1. Genes_by_length.pdf"))
par(mfrow=c(1,1))
hist(log(as.numeric(genes_info_mat[ ,"GeneLength"])+1), main = "Gene Size/Length", xlab = "Gene length/size (in log scale)", ylab = "Frequency")
dev.off()

temp_conn_mat <- only_t_conn_mat[only_t_conn_mat[ ,"ConnType"] == "cis" & !only_t_conn_mat[ ,"Dist"] %in% c("overlapped", "missing"), ]
pdf(paste0(fig_path, "2. cis_connection_by_dist.pdf"))
par(mfrow=c(1,3))
plot(log(as.numeric(temp_conn_mat[ ,"Mean"])+1), log(as.numeric(temp_conn_mat[ ,"Dist"])+1), main = "Distance vs #shared molecules \nonly for cis-connection", xlab = "log(#shared Molecules)", ylab = "log(Distance)", type = "p")
hist(log(as.numeric(temp_conn_mat[ ,"Dist"])+1), main = "Distribution by distance \nonly for cis-connection", xlab = "log(Distance)", ylab = "Frequency", labels = T)
hist(log(as.numeric(temp_conn_mat[log(as.numeric(temp_conn_mat[ ,"Mean"])+1) >=3 ,"Dist"])+1), main = "Distribution by distance \nonly for top cis-connection", xlab = "log(Distance)", ylab = "Frequency", labels = T)
dev.off()

pdf(paste0(fig_path, "3. #Shared_mole_in_samples.pdf"))
par(mfrow=c(1,3), pty="s")
plot(log(as.numeric(only_t_conn_mat[ ,"#Shared_mole_sample1"])), log(as.numeric(only_t_conn_mat[ ,"#Shared_mole_sample2"])), main = "tumor_sample1 \nvs \ntumor_sample2",xlab = "sample1 log(#shared molecules)", ylab = "sample2 log(#shared molecules)")
plot(log(as.numeric(only_t_conn_mat[ ,"#Shared_mole_sample1"])), log(as.numeric(only_t_conn_mat[ ,"Mean"])), main = "tumor_sample1 \nvs \n tumor_mean",xlab = "sample1 log(#shared molecules)", ylab = "mean log(#shared molecules)")
plot(log(as.numeric(only_t_conn_mat[ ,"#Shared_mole_sample2"])), log(as.numeric(only_t_conn_mat[ ,"Mean"])), main = "tumor_sample2 \nvs \n tumor_mean",xlab = "sample2 log(#shared molecules)", ylab = "mean log(#shared molecules)")
dev.off()

pdf(paste0(fig_path, "4. Count_Histogram.pdf"))
par(mfrow=c(1,2))
hist(log(as.numeric(only_t_conn_mat[only_t_conn_mat[ ,"ConnType"] == "cis", "Mean"])), main = "within same chromosome", xlab = "log(#shared molecules)", labels = T)
hist(log(as.numeric(only_t_conn_mat[only_t_conn_mat[ ,"ConnType"] == "trans", "Mean"])), main = "in different chromosome", xlab = "log(#shared molecules)", labels = T)
dev.off()

# round(exp(c(0:10)))
gene_gene_net <- only_t_conn_mat[ ,c("Gene1", "Gene2", "Mean", "ConnType", "Dist")]
colnames(gene_gene_net)[3] <- "Conn"
gene_gene_net[ ,"Conn"] <- round(log(as.numeric(gene_gene_net[ ,"Conn"])))
write.table(gene_gene_net, file = paste0(fig_path, "gene_gene_net.txt"), quote = F, sep = "\t", row.names = FALSE)

gene_deg_conn <- sort(table(gene_gene_net[ ,"Gene1"]), decreasing = T)
t_gene_deg_conn <- sort(table(gene_gene_net[ ,"Gene2"]), decreasing = T)
all_f_genes <- (union(only_t_conn_mat[ ,"Gene1"], only_t_conn_mat[ ,"Gene2"]))
degree_all_f_genes <- vector(mode = "integer", length = length(all_f_genes))
names(degree_all_f_genes) <- all_f_genes
degree_all_f_genes[names(gene_deg_conn)] <- degree_all_f_genes[names(gene_deg_conn)] + gene_deg_conn
degree_all_f_genes[names(t_gene_deg_conn)] <- degree_all_f_genes[names(t_gene_deg_conn)] + t_gene_deg_conn
degree_all_f_genes <- sort(degree_all_f_genes, decreasing = T)

pdf(paste0(fig_path, "5. Degree_of_Genes_Conn.pdf"))
par(mfrow=c(1,2))
plot(log(degree_all_f_genes), xlab = "Gene #", ylab = "log(Degree)")
hist(log(degree_all_f_genes), breaks = 50, labels = T, main = "", xlab = "log(Degree)")
dev.off()

#---------with threshold log(2(7), 3(20), 4(55))-------------------
rm(list = ls())
only_t_conn_mat <- get(load("data/t_analysis/t1_res/only_t1_conn_mat.RData"))
fig_path <- "data/t_analysis/t1_res/"

only_t_conn_mat <- only_t_conn_mat[only_t_conn_mat[ ,"Dist"] != "overlapped", ]

gene_gene_net <- only_t_conn_mat[ ,c("Gene1", "Gene2", "Mean", "ConnType", "Dist")]
colnames(gene_gene_net)[3] <- "Conn"
gene_gene_net[ ,"Conn"] <- round(log(as.numeric(gene_gene_net[ ,"Conn"])))
write.table(gene_gene_net, file = paste0(fig_path, "gene_gene_net.txt"), quote = F, sep = "\t", row.names = FALSE)

gene_gene_net <- only_t_conn_mat[as.numeric(only_t_conn_mat[ ,"Mean"]) >= 7, c("Gene1", "Gene2", "Mean", "ConnType", "Dist")]
colnames(gene_gene_net)[3] <- "Conn"
gene_gene_net[ ,"Conn"] <- round(log(as.numeric(gene_gene_net[ ,"Conn"])))
write.table(gene_gene_net, file = paste0(fig_path, "gene_gene_net_log2.txt"), quote = F, sep = "\t", row.names = FALSE)

gene_gene_net <- only_t_conn_mat[as.numeric(only_t_conn_mat[ ,"Mean"]) >= 20, c("Gene1", "Gene2", "Mean", "ConnType", "Dist")]
colnames(gene_gene_net)[3] <- "Conn"
gene_gene_net[ ,"Conn"] <- round(log(as.numeric(gene_gene_net[ ,"Conn"])))
write.table(gene_gene_net, file = paste0(fig_path, "gene_gene_net_log3.txt"), quote = F, sep = "\t", row.names = FALSE)

gene_gene_net <- only_t_conn_mat[as.numeric(only_t_conn_mat[ ,"Mean"]) >= 55, c("Gene1", "Gene2", "Mean", "ConnType", "Dist")]
colnames(gene_gene_net)[3] <- "Conn"
gene_gene_net[ ,"Conn"] <- round(log(as.numeric(gene_gene_net[ ,"Conn"])))
write.table(gene_gene_net, file = paste0(fig_path, "gene_gene_net_log4.txt"), quote = F, sep = "\t", row.names = FALSE)


# print(dim(only_t_conn_mat)) # total gene pairs
# print(table(only_t_conn_mat[ ,"ConnType"])) #total gene pairs within and between chromosomes
