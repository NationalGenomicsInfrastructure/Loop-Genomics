
#======================== Source code for Figures in the technical report of Loop Fusion Gene Analysis ========================

#************************ Figure 1: Distribution of gene lengths *********************************
#------------------
rm(list = ls())
library(ggplot2)
library(gridExtra)

genes_info_mat <- get(load("data/t_analysis/t1_res/t1_genes_info_mat.RData"))

genes_info_mat <- as.data.frame(genes_info_mat, stringsAsFactors = FALSE)
genes_info_mat[ ,c("startPos", "endPos", "GeneLength")] <- 
  apply(genes_info_mat[ ,c("startPos", "endPos", "GeneLength")], c(1,2), as.integer)

t_breaks <- seq(from = 0.0, to = 20.0, by = 2)
t_labels <- paste0(t_breaks, "\n(", round(exp(t_breaks)), ")")


g1 <- ggplot(genes_info_mat, aes(x=log(GeneLength))) + #, fill=ConnType
  geom_histogram(binwidth = 0.5, color = 'white') +
  xlab("Gene Length") +
  ylab("Count") +
  theme_bw() +
  theme(aspect.ratio = 1, legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  theme(legend.position = c(0.99,0.25), legend.justification = c("right", "bottom")) +
  theme(axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 14)) +
  labs(subtitle = "Tumor1") +
  scale_x_continuous(breaks=t_breaks, labels=t_labels)

g1


genes_info_mat <- get(load("data/t_analysis/t2_res/t2_genes_info_mat.RData"))

genes_info_mat <- as.data.frame(genes_info_mat, stringsAsFactors = FALSE)
genes_info_mat[ ,c("startPos", "endPos", "GeneLength")] <- 
  apply(genes_info_mat[ ,c("startPos", "endPos", "GeneLength")], c(1,2), as.integer)

g2 <- ggplot(genes_info_mat, aes(x=log(GeneLength))) + #, fill=ConnType
  geom_histogram(binwidth = 0.5, color = 'white') +
  #stat_bin(aes(label=..count..), geom="text", binwidth = 0.5, position=position_stack(vjust = 0.5), angle = c(90), size=5, color="orange") +
  xlab("Gene Length") +
  ylab("Count") +
  theme_bw() +
  theme(aspect.ratio = 1, legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  theme(legend.position = c(0.99,0.25), legend.justification = c("right", "bottom")) +
  theme(axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 14)) +
  labs(subtitle = "Tumor2") +
  scale_x_continuous(breaks=t_breaks, labels=t_labels)

g2

genes_info_mat <- get(load("data/t_analysis/t3_res/t3_genes_info_mat.RData"))

genes_info_mat <- as.data.frame(genes_info_mat, stringsAsFactors = FALSE)
genes_info_mat[ ,c("startPos", "endPos", "GeneLength")] <- 
  apply(genes_info_mat[ ,c("startPos", "endPos", "GeneLength")], c(1,2), as.integer)


g3 <- ggplot(genes_info_mat, aes(x=log(GeneLength))) + #, fill=ConnType
  geom_histogram(binwidth = 0.5, color = 'white') +
  #stat_bin(aes(label=..count..), geom="text", binwidth = 0.5, position=position_stack(vjust = 0.5), angle = c(90), size=5, color="orange") +
  xlab("Gene Length") +
  ylab("Count") +
  theme_bw() +
  theme(aspect.ratio = 1, legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  theme(legend.position = c(0.99,0.25), legend.justification = c("right", "bottom")) +
  theme(axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 14)) +
  labs(subtitle = "Tumor3") +
  scale_x_continuous(breaks=t_breaks, labels=t_labels)

#pdf(file = "test.pdf")
g3

grid.arrange(g1, g2, g3, nrow=1)

#------------------
#************************ Figure 2: Read counts (#shared moledules) between replicates *********************************
#------------------

rm(list = ls())
library(ggplot2)
library(gridExtra)

only_t_conn_mat <- get(load("data/t_analysis/t3_res/only_t3_conn_mat.RData"))

only_t_conn_mat <- as.data.frame(only_t_conn_mat, stringsAsFactors = FALSE)
only_t_conn_mat[ ,c("#Shared_mole_sample1", "#Shared_mole_sample2", "Mean", "Median")] <- 
  apply(only_t_conn_mat[ ,c("#Shared_mole_sample1", "#Shared_mole_sample2", "Mean", "Median")], c(1,2), as.integer)


t_breaks <- seq(from = 0.0, to = 10.0, by = 1.5)
t_labels <- paste0(t_breaks, "\n(", round(exp(t_breaks)), ")")

r2_val <- round(cor(log(only_t_conn_mat$`#Shared_mole_sample1`), log(only_t_conn_mat$`#Shared_mole_sample2`))^2, 2)
g1 <- ggplot(only_t_conn_mat, aes(x=log(`#Shared_mole_sample1`), y=log(`#Shared_mole_sample2`))) + 
  geom_point(size=0.5, aes(col=ConnType)) + 
  xlab("Sample 1") +
  ylab("Sample 2") +
  theme_bw() +
  theme(aspect.ratio = 1, legend.text = element_text(size = 8), legend.title = element_text(size = 9)) +
  theme(legend.position = c(0.99,0.25), legend.justification = c("right", "bottom")) +
  annotate("text", x = 2, y = 7, label = paste0("italic(R)^2 == ", r2_val), parse=T) +
  labs(subtitle = "")

g1 <- g1 + scale_x_continuous(breaks=t_breaks, labels=t_labels) + scale_y_continuous(breaks=t_breaks, labels=t_labels)

g1

r2_val <- round(cor(log(only_t_conn_mat$`#Shared_mole_sample1`), log(only_t_conn_mat$Mean))^2, 2)
g2 <- ggplot(only_t_conn_mat, aes(x=log(`#Shared_mole_sample1`), y=log(Mean))) + 
  geom_point(size=0.5, aes(col=ConnType)) + 
  xlab("Sample 1") +
  ylab("Mean") +
  theme_bw() +
  theme(aspect.ratio = 1, legend.text = element_text(size = 8), legend.title = element_text(size = 9)) +
  theme(legend.position = c(0.99,0.25), legend.justification = c("right", "bottom")) +
  annotate("text", x = 2, y = 7, label = paste0("italic(R)^2 == ", r2_val), parse=T) +
  labs(subtitle = "")

g2 <- g2 + scale_x_continuous(breaks=t_breaks, labels=t_labels) + scale_y_continuous(breaks=t_breaks, labels=t_labels)

g2

r2_val <- round(cor(log(only_t_conn_mat$`#Shared_mole_sample2`), log(only_t_conn_mat$Mean))^2, 2)
g3 <- ggplot(only_t_conn_mat, aes(x=log(`#Shared_mole_sample2`), y=log(Mean))) + 
  geom_point(size=0.5, aes(col=ConnType)) + 
  xlab("Sample 2") +
  ylab("Mean") +
  theme_bw() +
  theme(aspect.ratio = 1, legend.text = element_text(size = 8), legend.title = element_text(size = 9)) +
  theme(legend.position = c(0.99,0.25), legend.justification = c("right", "bottom")) +
  annotate("text", x = 2, y = 7, label = paste0("italic(R)^2 == ", r2_val), parse=T) +
  labs(subtitle = "")

g3 <- g3 + scale_x_continuous(breaks=t_breaks, labels=t_labels) + scale_y_continuous(breaks=t_breaks, labels=t_labels)

grid.arrange(g1, g2, g3, nrow=1)

# pdf(file = "test.pdf")
# grid.arrange(g1, g2, g3, nrow=1)
# dev.off()


#------------------
#************************ Figure 3: Read counts (#shared moledules) between Tumors *********************************
#------------------

rm(list = ls())
library(ggplot2)
library(gridExtra)

Common_edges <- function(l_edges_p, l_rnames_p) {
  
  t_common_edges <- Reduce(intersect, l_edges_p)
  t_common_edges <- strsplit(t_common_edges, split = "-vs-")
  t_common_edges <- unlist(lapply(t_common_edges, function(x) {return(paste(x[2], x[1], sep = "-vs-"))}))
  
  common_edges <- Reduce(intersect, l_rnames_p) 
  common_edges <- union(common_edges, t_common_edges)
  
  return(common_edges)
}

only_t_conn_mat_1 <- get(load("data/t_analysis/t1_res/only_t1_conn_mat.RData"))
only_t_conn_mat_2 <- get(load("data/t_analysis/t2_res/only_t2_conn_mat.RData"))
only_t_conn_mat_3 <- get(load("data/t_analysis/t3_res/only_t3_conn_mat.RData"))

l_edges <- list(only_t_conn_mat_1[ ,"Edge"], only_t_conn_mat_2[ ,"Edge"])
l_rnames <- list(rownames(only_t_conn_mat_1), rownames(only_t_conn_mat_2))
common_edges <- Common_edges(l_edges, l_rnames)

only_t_conn_mat <- cbind(only_t_conn_mat_1[common_edges ,c("Mean")], only_t_conn_mat_2[common_edges ,c("Mean", "ConnType")])
colnames(only_t_conn_mat) <- c("Mean1", "Mean2", "ConnType")

only_t_conn_mat <- as.data.frame(only_t_conn_mat, stringsAsFactors = FALSE)
only_t_conn_mat[ ,c("Mean1", "Mean2")] <- apply(only_t_conn_mat[ ,c("Mean1", "Mean2")], c(1,2), as.integer)


t_breaks <- seq(from = 0.0, to = 10.0, by = 1.5)
t_labels <- paste0(t_breaks, "\n(", round(exp(t_breaks)), ")")

r2_val <- round(cor(log(only_t_conn_mat$Mean1), log(only_t_conn_mat$Mean2))^2, 2)
g1 <- ggplot(only_t_conn_mat, aes(x=log(Mean1), y=log(Mean2))) + 
  geom_point(size=0.5, aes(col=ConnType)) + 
  xlab("Tumor1") +
  ylab("Tumor2") +
  theme_bw() +
  theme(aspect.ratio = 1, legend.text = element_text(size = 8), legend.title = element_text(size = 9)) +
  theme(legend.position = c(0.99,0.25), legend.justification = c("right", "bottom")) +
  annotate("text", x = 2, y = 7, label = paste0("italic(R)^2 == ", r2_val), parse=T) +
  labs(subtitle = "")

g1 <- g1 + scale_x_continuous(breaks=t_breaks, labels=t_labels) + scale_y_continuous(breaks=t_breaks, labels=t_labels)

g1

l_edges <- list(only_t_conn_mat_1[ ,"Edge"], only_t_conn_mat_3[ ,"Edge"])
l_rnames <- list(rownames(only_t_conn_mat_1), rownames(only_t_conn_mat_3))
common_edges <- Common_edges(l_edges, l_rnames)

only_t_conn_mat <- cbind(only_t_conn_mat_1[common_edges ,c("Mean")], only_t_conn_mat_3[common_edges ,c("Mean", "ConnType")])
colnames(only_t_conn_mat) <- c("Mean1", "Mean2", "ConnType")

only_t_conn_mat <- as.data.frame(only_t_conn_mat, stringsAsFactors = FALSE)
only_t_conn_mat[ ,c("Mean1", "Mean2")] <- apply(only_t_conn_mat[ ,c("Mean1", "Mean2")], c(1,2), as.integer)

r2_val <- round(cor(log(only_t_conn_mat$Mean1), log(only_t_conn_mat$Mean2))^2, 2)
g2 <- ggplot(only_t_conn_mat, aes(x=log(Mean1), y=log(Mean2))) + 
  geom_point(size=0.5, aes(col=ConnType)) + 
  xlab("Tumor1") +
  ylab("Tumor3") +
  theme_bw() +
  theme(aspect.ratio = 1, legend.text = element_text(size = 8), legend.title = element_text(size = 9)) +
  theme(legend.position = c(0.99,0.25), legend.justification = c("right", "bottom")) +
  annotate("text", x = 2, y = 7, label = paste0("italic(R)^2 == ", r2_val), parse=T) +
  labs(subtitle = "")

g2 <- g2 + scale_x_continuous(breaks=t_breaks, labels=t_labels) + scale_y_continuous(breaks=t_breaks, labels=t_labels)

g2

l_edges <- list(only_t_conn_mat_2[ ,"Edge"], only_t_conn_mat_3[ ,"Edge"])
l_rnames <- list(rownames(only_t_conn_mat_2), rownames(only_t_conn_mat_3))
common_edges <- Common_edges(l_edges, l_rnames)

only_t_conn_mat <- cbind(only_t_conn_mat_2[common_edges ,c("Mean")], only_t_conn_mat_3[common_edges ,c("Mean", "ConnType")])
colnames(only_t_conn_mat) <- c("Mean1", "Mean2", "ConnType")

only_t_conn_mat <- as.data.frame(only_t_conn_mat, stringsAsFactors = FALSE)
only_t_conn_mat[ ,c("Mean1", "Mean2")] <- apply(only_t_conn_mat[ ,c("Mean1", "Mean2")], c(1,2), as.integer)

r2_val <- round(cor(log(only_t_conn_mat$Mean1), log(only_t_conn_mat$Mean2))^2, 2)
g3 <- ggplot(only_t_conn_mat, aes(x=log(Mean1), y=log(Mean2))) + 
  geom_point(size=0.5, aes(col=ConnType)) + 
  xlab("Tumor2") +
  ylab("Tumor3") +
  theme_bw() +
  theme(aspect.ratio = 1, legend.text = element_text(size = 8), legend.title = element_text(size = 9)) +
  theme(legend.position = c(0.99,0.25), legend.justification = c("right", "bottom")) +
  annotate("text", x = 2, y = 7, label = paste0("italic(R)^2 == ", r2_val), parse=T) +
  labs(subtitle = "")

g3 <- g3 + scale_x_continuous(breaks=t_breaks, labels=t_labels) + scale_y_continuous(breaks=t_breaks, labels=t_labels)

grid.arrange(g1, g2, g3, nrow=1)




#------------------
#************************ Figure 4: Count Histogram *********************************
#------------------

rm(list = ls())
library(ggplot2)
library(gridExtra)

only_t_conn_mat <- get(load("data/t_analysis/t3_res/only_t3_conn_mat.RData"))
gene_gene_net <- only_t_conn_mat[ ,c("Gene1", "Gene2", "Mean", "ConnType", "Dist")] # for next figure

only_t_conn_mat <- as.data.frame(only_t_conn_mat, stringsAsFactors = FALSE)
only_t_conn_mat[ ,c("#Shared_mole_sample1", "#Shared_mole_sample2", "Mean", "Median")] <- 
  apply(only_t_conn_mat[ ,c("#Shared_mole_sample1", "#Shared_mole_sample2", "Mean", "Median")], c(1,2), as.integer)

t_breaks <- seq(from = 0.0, to = 10.0, by = 0.5)
t_labels <- paste0(t_breaks, "\n(", round(exp(t_breaks)), ")")


g1 <- ggplot(only_t_conn_mat, aes(x=log(Mean), fill=ConnType)) + #, fill=ConnType
  geom_histogram(binwidth = 0.5, color = 'white') +
  stat_bin(aes(label=ifelse(..count.. > 0, ..count.., ""), y=ifelse(..count..<=50, 100, ..count..)), geom="text", binwidth = 0.5, position=position_stack(vjust = 0.5), angle = c(45), size=6) +
  xlab("#Shared Molecules") +
  ylab("Count") +
  theme_bw() +
  theme(aspect.ratio = 1, legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  theme(legend.position = c(0.99,0.25), legend.justification = c("right", "bottom")) +
  theme(axis.text = element_text(size = 10, face = "bold"), axis.title = element_text(size = 12)) +
  labs(subtitle = "") +
  scale_x_continuous(breaks=t_breaks, labels=t_labels)

#pdf(file = "test.pdf")
g1
#dev.off()

colnames(gene_gene_net)[3] <- "Conn"
gene_gene_net[ ,"Conn"] <- round(log(as.numeric(gene_gene_net[ ,"Conn"])))
gene_deg_conn <- sort(table(gene_gene_net[ ,"Gene1"]), decreasing = T)
t_gene_deg_conn <- sort(table(gene_gene_net[ ,"Gene2"]), decreasing = T)
all_f_genes <- (union(only_t_conn_mat[ ,"Gene1"], only_t_conn_mat[ ,"Gene2"]))
degree_all_f_genes <- vector(mode = "integer", length = length(all_f_genes))
names(degree_all_f_genes) <- all_f_genes
degree_all_f_genes[names(gene_deg_conn)] <- degree_all_f_genes[names(gene_deg_conn)] + gene_deg_conn
degree_all_f_genes[names(t_gene_deg_conn)] <- degree_all_f_genes[names(t_gene_deg_conn)] + t_gene_deg_conn
degree_all_f_genes <- sort(degree_all_f_genes, decreasing = T)

degree_all_f_genes <- as.data.frame(degree_all_f_genes)

t_breaks <- seq(from = 0.0, to = 10.0, by = 0.5)
t_labels <- paste0(t_breaks, "\n(", round(exp(t_breaks)), ")")

g2 <- ggplot(degree_all_f_genes, aes(x=log(degree_all_f_genes))) + #, fill=ConnType
  geom_histogram(binwidth = 0.5, color = "white", fill="lightblue") +
  stat_bin(aes(label=..count..), geom="text", binwidth = 0.5, position=position_stack(vjust = 0.5), angle = c(45), size=6, color="black") +
  xlab("Degree (#fusions) ") +
  ylab("Counts (Gene)") +
  theme_bw() +
  theme(aspect.ratio = 1, legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  theme(legend.position = c(0.99,0.25), legend.justification = c("right", "bottom")) +
  theme(axis.text = element_text(size = 10, face = "bold"), axis.title = element_text(size = 12)) +
  labs(subtitle = "") +
  scale_x_continuous(breaks=t_breaks, labels=t_labels)

grid.arrange(g1, g2, nrow=1, ncol=2)


#------------------
#************************ Figure 5: The proximity of cis-connections and distribution of cis-connections by distance *********************************
#------------------

rm(list = ls())
library(ggplot2)
library(gridExtra)

only_t_conn_mat <- get(load("data/t_analysis/t1_res/only_t1_conn_mat.RData"))

temp_conn_mat <- only_t_conn_mat[only_t_conn_mat[ ,"ConnType"] == "cis" & !only_t_conn_mat[ ,"Dist"] %in% c("overlapped", "missing"), ]

temp_conn_mat <- as.data.frame(temp_conn_mat, stringsAsFactors = FALSE)
temp_conn_mat[ ,c("#Shared_mole_sample1", "#Shared_mole_sample2", "Mean", "Median", "Dist")] <- 
  apply(temp_conn_mat[ ,c("#Shared_mole_sample1", "#Shared_mole_sample2", "Mean", "Median", "Dist")], c(1,2), as.integer)

t_breaks <- seq(from = 0.0, to = 20.0, by = 0.5)
t_labels <- paste0(t_breaks, "\n(", round(exp(t_breaks)), ")")

g1 <- ggplot(temp_conn_mat, aes(x=log(Mean+1), y=log(Dist+1))) + 
  geom_point(size=0.5) + 
  xlab("#Shared Molecules") +
  ylab("Distance") +
  theme_bw() +
  theme(aspect.ratio = 1, legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  theme(axis.text = element_text(size = 10, face = "bold"), axis.title = element_text(size = 14)) +
  labs(subtitle = "#Shared Molecules vs Distance") +
  scale_x_continuous(breaks=t_breaks, labels=t_labels) +
  scale_y_continuous(breaks=t_breaks, labels=t_labels)

t_breaks <- seq(from = 0.0, to = 20.0, by = 2)
t_labels <- paste0(t_breaks, "\n(", round(exp(t_breaks)), ")")

g2 <- ggplot(temp_conn_mat, aes(x=log(Dist))) + 
  geom_histogram(binwidth = 2, color="white") + 
  stat_bin(aes(label=..count..), geom="text", binwidth = 2, position=position_stack(vjust = 0.5), angle = c(45), size=6, color="green") +
  xlab("Distance") +
  ylab("Counts") +
  theme_bw() +
  theme(aspect.ratio = 1, legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  theme(legend.position = c(0.99,0.25), legend.justification = c("right", "bottom")) +
  theme(axis.text = element_text(size = 10, face = "bold"), axis.title = element_text(size = 12)) +
  labs(subtitle = "Distribution by distance (only for cis-connections)") +
  scale_x_continuous(breaks=t_breaks, labels=t_labels) 

grid.arrange(g1, g2, nrow=1, ncol=2)

#------------------
#************************ Figure : Generate files for cytoscpae figure with threshold log(2(7), 3(20), 4(55)) *********************************
#------------------
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

#------------------
#************************ Table : cis vs overlapped vs non-overlapped *********************************
#------------------
rm(list = ls())

only_t_conn_mat <- get(load("data/t_analysis/t3_res/only_t3_conn_mat.RData"))

temp_conn_mat <- only_t_conn_mat[only_t_conn_mat[ ,"ConnType"] == "cis", ]
print(paste0("Total cis-connections in Tumor ", nrow(temp_conn_mat)))
print(paste0("Total cis-connections where genes are non-overlapped in Tumor ", nrow(temp_conn_mat[!temp_conn_mat[ ,"Dist"] %in% c("overlapped", "missing"), ])))
print(paste0("Total cis-connections where genes are overlapped in Tumor ", nrow(temp_conn_mat[temp_conn_mat[ ,"Dist"] %in% c("overlapped"), ])))
print(paste0("Total cis-connections where genes are positions were missing in Tumor ", nrow(temp_conn_mat[temp_conn_mat[ ,"Dist"] %in% c("missing"), ])))

        

