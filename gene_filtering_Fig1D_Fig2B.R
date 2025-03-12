# Generate data tables for heatmaps in Fig1D and Fig2B
# Load libraries
library(ggplot2)
library(reshape2)
library(reshape2)

# Load data from Mageck, gene level data
L11vsNons1 <- read.csv("./L11vsNons1_test.gene_summary.txt", sep = "\t",)
L51vsNons1 <- read.csv("./L51vsNons1_test.gene_summary.txt", sep = "\t",)
L12vsNons2 <- read.csv("./L12vsNons2_test.gene_summary.txt", sep = "\t",)
L52vsNons2 <- read.csv("./L52vsNons2_test.gene_summary.txt", sep = "\t",)

H11vsNons1 <- read.csv("./H11vsNons1_test.gene_summary.txt", sep = "\t",)
H12vsNons2 <- read.csv("./H12vsNons2_test.gene_summary.txt", sep = "\t",)
H51vsNons1 <- read.csv("./H51vsNons1_test.gene_summary.txt", sep = "\t",)
H52vsNons2 <- read.csv("./H52vsNons2_test.gene_summary.txt", sep = "\t",)

# Produce list based on fdr <0.05 or twice <0.25
L11vsNons1 <- subset(L11vsNons1, pos.fdr < 0.25)[,c(1,9:14)]
colnames(L11vsNons1) <- c("id", "L11_score", "L11_p.value", "L11_fdr", "L11_rank", "L11_goodsrna", "L11_lfc")
L51vsNons1 <- subset(L51vsNons1, pos.fdr < 0.25)[,c(1,9:14)]
colnames(L51vsNons1) <- c("id", "L51_score", "L51_p.value", "L51_fdr", "L51_rank", "L51_goodsrna", "L51_lfc")
L12vsNons2 <- subset(L12vsNons2, pos.fdr < 0.25)[,c(1,9:14)]
colnames(L12vsNons2) <- c("id", "L12_score", "L12_p.value", "L12_fdr", "L12_rank", "L12_goodsrna", "L12_lfc")
L52vsNons2 <- subset(L52vsNons2, pos.fdr < 0.25)[,c(1,9:14)]
colnames(L52vsNons2) <- c("id", "L52_score", "L52_p.value", "L52_fdr", "L52_rank", "L52_goodsrna", "L52_lfc")
Low <- merge(L11vsNons1, L51vsNons1, by = "id", all = TRUE)
Low_2 <- merge(L12vsNons2, L52vsNons2, by = "id", all = TRUE)
Low <- merge(Low, Low_2, by = "id", all = TRUE)
# write.csv2(Low, "/Users/fernkorn/Desktop/20211215_low_spry4_enriched_OR-gate_fdr<0.05.csv")
Low <- subset(Low, rowSums(is.na(Low[,])) < 15)
# write.csv2(Low, "/Users/fernkorn/Desktop/20211215_low_spry4_enriched_AND2-gate_fdr<0.05.csv")

L11vsNons1 <- subset(L11vsNons1, pos.fdr < 0.05)[,c(1,9:14)]
colnames(L11vsNons1) <- c("id", "L11_score", "L11_p.value", "L11_fdr", "L11_rank", "L11_goodsrna", "L11_lfc")
L51vsNons1 <- subset(L51vsNons1, pos.fdr < 0.05)[,c(1,9:14)]
colnames(L51vsNons1) <- c("id", "L51_score", "L51_p.value", "L51_fdr", "L51_rank", "L51_goodsrna", "L51_lfc")
L12vsNons2 <- subset(L12vsNons2, pos.fdr < 0.05)[,c(1,9:14)]
colnames(L12vsNons2) <- c("id", "L12_score", "L12_p.value", "L12_fdr", "L12_rank", "L12_goodsrna", "L12_lfc")
L52vsNons2 <- subset(L52vsNons2, pos.fdr < 0.05)[,c(1,9:14)]
colnames(L52vsNons2) <- c("id", "L52_score", "L52_p.value", "L52_fdr", "L52_rank", "L52_goodsrna", "L52_lfc")
Low_3 <- merge(L11vsNons1, L51vsNons1, by = "id", all = TRUE)
Low_4 <- merge(L12vsNons2, L52vsNons2, by = "id", all = TRUE)
Low_2 <- merge(Low_3, Low_4, by = "id", all = TRUE)
Low_total <- merge(Low, Low_2,  by = "id", all = TRUE)


#high
H11vsNons1 <- subset(H11vsNons1, pos.fdr < 0.25)[,c(1,9:14)]
colnames(H11vsNons1) <- c("id", "H11_score", "H11_p.value", "H11_fdr", "H11_rank", "H11_goodsrna", "H11_lfc")
H12vsNons2 <- subset(H12vsNons2, pos.fdr < 0.25)[,c(1,9:14)]
colnames(H12vsNons2) <- c("id", "H51_score", "H51_p.value", "H51_fdr", "H51_rank", "H51_goodsrna", "H51_lfc")
H51vsNons1 <- subset(H51vsNons1, pos.fdr < 0.25)[,c(1,9:14)]
colnames(H51vsNons1) <- c("id", "H12_score", "H12_p.value", "H12_fdr", "H12_rank", "H12_goodsrna", "H12_lfc")
H52vsNons2 <- subset(H52vsNons2, pos.fdr < 0.25)[,c(1,9:14)]
colnames(H52vsNons2) <- c("id", "H52_score", "H52_p.value", "H52_fdr", "H52_rank", "H52_goodsrna", "H52_lfc")
High <- merge(H11vsNons1, H12vsNons2, by = "id", all = TRUE)
High_2 <- merge(H51vsNons1, H52vsNons2, by = "id", all = TRUE)
High <- merge(High, High_2, by = "id", all = TRUE)
High <- subset(High, rowSums(is.na(High[,])) < 15)
H11vsNons1 <- subset(H11vsNons1, pos.fdr < 0.05)[,c(1,9:14)]
colnames(H11vsNons1) <- c("id", "H11_score", "H11_p.value", "H11_fdr", "H11_rank", "H11_goodsrna", "H11_lfc")
H12vsNons2 <- subset(H12vsNons2, pos.fdr < 0.05)[,c(1,9:14)]
colnames(H12vsNons2) <- c("id", "H51_score", "H51_p.value", "H51_fdr", "H51_rank", "H51_goodsrna", "H51_lfc")
H51vsNons1 <- subset(H51vsNons1, pos.fdr < 0.05)[,c(1,9:14)]
colnames(H51vsNons1) <- c("id", "H12_score", "H12_p.value", "H12_fdr", "H12_rank", "H12_goodsrna", "H12_lfc")
H52vsNons2 <- subset(H52vsNons2, pos.fdr < 0.05)[,c(1,9:14)]
colnames(H52vsNons2) <- c("id", "H52_score", "H52_p.value", "H52_fdr", "H52_rank", "H52_goodsrna", "H52_lfc")
High_3 <- merge(H11vsNons1, H12vsNons2, by = "id", all = TRUE)
High_4 <- merge(H51vsNons1, H52vsNons2, by = "id", all = TRUE)
High_2 <- merge(High_3, High_4, by = "id", all = TRUE)
High_total <- merge(High, High_2,  by = "id", all = TRUE)


# Extract complete information for list of genes
rownames(L11vsNons1) <- L11vsNons1$id
rownames(L51vsNons1) <- L51vsNons1$id
rownames(L12vsNons2) <- L12vsNons2$id
rownames(L52vsNons2) <- L52vsNons2$id
Low_all_data <- data.frame(L11vsNons1[Low_total$id,c(1,9:14)],L51vsNons1[Low_total$id,c(9:14)],L12vsNons2[Low_total$id,c(9:14)],L52vsNons2[Low_total$id,c(9:14)])
colnames(Low_all_data) <- c("id", "L11_score", "L11_p.value", "L11_fdr", "L11_rank", "L11_goodsrna", "L11_lfc", "L51_score", "L51_p.value", "L51_fdr", "L51_rank", "L51_goodsrna", "L51_lfc", 
                            "L12_score", "L12_p.value", "L12_fdr", "L12_rank", "L12_goodsrna", "L12_lfc", "L52_score", "L52_p.value", "L52_fdr", "L52_rank", "L52_goodsrna", "L52_lfc")
write.csv2(Low_all_data, "./20220113_low_all_data_spry4_enriched_AND2-gate_fdr<0.25_OR_fdr<0.05.csv", row.names = FALSE)

rownames(H11vsNons1) <- H11vsNons1$id
rownames(H12vsNons2) <- H12vsNons2$id
rownames(H51vsNons1) <- H51vsNons1$id
rownames(H52vsNons2) <- H52vsNons2$id
High_all_data <- data.frame(H11vsNons1[High_total$id,c(1,9:14)],H12vsNons2[High_total$id,c(9:14)],H51vsNons1[High_total$id,c(9:14)],H52vsNons2[High_total$id,c(9:14)])
colnames(High_all_data) <- c("id", "H11_score", "H11_p.value", "H11_fdr", "H11_rank", "H11_goodsrna", "H11_lfc", "H51_score", "H51_p.vaHue", "H51_fdr", "H51_rank", "H51_goodsrna", "H51_lfc", 
                            "H12_score", "H12_p.value", "H12_fdr", "H12_rank", "H12_goodsrna", "H12_lfc", "H52_score", "H52_p.vaHue", "H52_fdr", "H52_rank", "H52_goodsrna", "H52_lfc")
write.csv2(High_all_data, "./20220113_high_all_data_spry4_enriched_AND2-gate_fdr<0.25_OR_fdr<0.05.csv", row.names = FALSE)


