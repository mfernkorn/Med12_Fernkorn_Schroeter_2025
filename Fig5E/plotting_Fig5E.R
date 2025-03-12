#PLot Figure 5E
#load libraries
library(DESeq2)
library(ggplot2)
library(reshape2)
library(ggh4x)
library(pheatmap)
library(stringr)
library(dplyr)

# First part is same as for Figure 5D
# For Differential gene expression analysis DESeq will be used, which takes unnormalized raw counts of reads as an input
# Load Raw counts
cts <- read.csv(file = "Annotated Probe Report for All Probes_Raw_Counts.txt",sep="\t")
# Make Probe names unique by including Positional information
rownames(cts) <- paste(cts$Probe, cts$Chromosome, cts$Start, cts$End, sep = "_") 
cts <- as.matrix(cts[,13:24])

# Define column with all metadata information (FGF state and Med12 genotype)
FGF_state = rep(c(rep("_2i", 3),rep("N2", 3)),2)
Med12_state = c(rep("wt", 6), rep("Med12_KO", 6))
coldata <- data.frame(State = paste(FGF_state, Med12_state, sep = "_"),
                      row.names = colnames(cts))

# Combine Counts and Metadata to DESeq object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ State)

# Filter lowly expressed genes (Raw count >= 100 in at least one sample)
keep <- rowSums(counts(dds) >= 100) >= 1
dds <- dds[keep,]

# Run Deseq2 and 
dds <- DESeq(dds)

# extract results as a list of all combinations+ applying log foldchange shrinkage
results <- list()
for (i in 1:6){
  cur_name <- paste(combn(unique(coldata$State)[c(2,4,1,3)],2)[,i][1], 
                    combn(unique(coldata$State)[c(2,4,1,3)],2)[,i][2], sep = "_vs_")
  results[cur_name] <- lfcShrink(dds, contrast = c("State", combn(unique(coldata$State)[c(2,4,1,3)],2)[,i]), type="ashr")
}

# Produce volcano plots
Comparison <- rbind(data.frame("Gene" = rownames(results$N2_wt_vs_N2_Med12_KO), 
                               "Comparison" = "wt vs. Med12 KO in N2B27",
                               "FC" = results$N2_wt_vs_N2_Med12_KO$log2FoldChange,
                               "p" = results$N2_wt_vs_N2_Med12_KO$padj),
                    data.frame("Gene" = rownames(results$N2_wt_vs_N2_Med12_KO),
                               "Comparison" = "wt vs. Med12 KO in 2i",
                               "FC" = results$`_2i_wt_vs__2i_Med12_KO`$log2FoldChange,
                               "p" = results$`_2i_wt_vs__2i_Med12_KO`$padj),
                    data.frame("Gene" = rownames(results$N2_wt_vs__2i_wt),
                               "Comparison" = "2i vs N2B27 in wt",
                               "FC" = results$N2_wt_vs__2i_wt$log2FoldChange,
                               "p" = results$N2_wt_vs__2i_wt$padj),
                    data.frame("Gene" = rownames(results$N2_Med12_KO_vs__2i_Med12_KO),
                               "Comparison" = "2i vs N2B27 in Med12 KO",
                               "FC" = results$N2_Med12_KO_vs__2i_Med12_KO$log2FoldChange,
                               "p" = results$N2_Med12_KO_vs__2i_Med12_KO$padj))
Comparison$Gene <- sub("_.*", "", Comparison$Gene) 
Comparison$Comparison <- factor(Comparison$Comparison, 
                                levels = c("wt vs. Med12 KO in 2i", "wt vs. Med12 KO in N2B27",
                                           "2i vs N2B27 in wt","2i vs N2B27 in Med12 KO"),
                                ordered = TRUE)

table(subset(Comparison, FC>1 & p<0.01)$Comparison)
table(subset(Comparison, FC<(-1) & p<0.01)$Comparison)

Comparison$signif <- ifelse((Comparison$FC<(-1) | Comparison$FC>1) & Comparison$p<0.01, 
                            "significant", "not significant")
y_limit <- 1e+100
Comparison_edge <- subset(Comparison, p<(1/y_limit))
Comparison_edge$p <-1/ y_limit


DEG_2i_vs_N2B27_Med12KO <- subset(Comparison, (FC<(-1) | FC>1) & p<0.01 & Comparison == "2i vs N2B27 in Med12 KO")
DEG_2i_vs_N2B27_wt <- subset(Comparison, (FC<(-1) | FC>1) & p<0.01 & Comparison == "2i vs N2B27 in wt")

DEGs <- merge(DEG_2i_vs_N2B27_Med12KO, DEG_2i_vs_N2B27_wt, by = "Gene")

pluripotency_markers <- c("Klf4", "Nanog", "Esrrb", "Tbx3", "Tfcp2l1","Prdm14", "Zfp42")
Top10down <- DEGs[order(DEGs$FC.y)[1:100], ]$Gene
Top10down_mutant <- DEGs[order(DEGs$FC.x)[1:100], ]$Gene
Top10up <- DEGs[order(DEGs$FC.y)[(length(DEGs$FC.y)-99):length(DEGs$FC.y)], ]$Gene
Top10up_mutant <- DEGs[-order(DEGs$FC.x)[(length(DEGs$FC.y)-99):length(DEGs$FC.y)], ]$Gene


# Load own raw data
Med12comp <- read.csv(file = "./RNA-seq_Comparison_with_Lackner_et_al/Data/Annotated Probe Report for All Probes_Raw_Counts.txt",sep="\t")
Med12comp <- Med12comp[!duplicated(Med12comp$Probe),] # Remove duplicate Genes
Med12comp_raw <- Med12comp[,c(1,13:24)]
rownames(Med12comp_raw) <- Med12comp_raw$Probe
Med12comp_raw$Probe <- NULL

# Define column with all metadata information (FGF state and Med12 genotype)
FGF_state = rep(c(rep("twoi", 3),rep("N2", 3)),2)
Med12_state = c(rep("wt", 6), rep("Med12_KO", 6))
coldata <- data.frame(State = paste(FGF_state, Med12_state, sep = "_"),
                      row.names = colnames(Med12comp_raw))

# Combine Counts and Metadata to DESeq object
dds <- DESeqDataSetFromMatrix(countData = Med12comp_raw,
                              colData = coldata,
                              design = ~ State)

# Filter lowly expressed genes (Raw count >= 100 in at least one sample)
keep <- rowSums(counts(dds) >= 100) >= 1
dds <- dds[keep,]

# Run Deseq2 and extract results as a list of all combinations+ applying log foldchange shrinkage
dds <- DESeq(dds)
results <- list()
for (i in 1:6){
  cur_name <- paste(combn(unique(coldata$State)[c(2,4,1,3)],2)[,i][1], 
                    combn(unique(coldata$State)[c(2,4,1,3)],2)[,i][2], sep = "_vs_")
  results[cur_name] <- lfcShrink(dds, contrast = c("State", combn(unique(coldata$State)[c(2,4,1,3)],2)[,i]), type="ashr")
}
Comparison_own <- rbind(data.frame("Gene" = rownames(results$N2_wt_vs_N2_Med12_KO), 
                                   "FC_wt" = results$N2_wt_vs_twoi_wt$log2FoldChange,
                                   "FC_Med12KO" = results$N2_Med12_KO_vs_twoi_Med12_KO$log2FoldChange))

Comparison_own_N2vs2i <- rbind(data.frame("Gene" = rownames(results$N2_wt_vs_twoi_wt), 
                                          "FC_wt" = results$N2_wt_vs_twoi_wt$log2FoldChange,
                                          "FC_Med12KO" = results$N2_Med12_KO_vs_twoi_Med12_KO$log2FoldChange))

Comparison_own_wtvsKO <- rbind(data.frame("Gene" = rownames(results$twoi_wt_vs_twoi_Med12_KO), 
                                          "FC_2i" = results$twoi_wt_vs_twoi_Med12_KO$log2FoldChange,
                                          "FC_N2B27" = results$N2_wt_vs_N2_Med12_KO$log2FoldChange))

Med12_TPM <- read.csv(file = "./RNA-seq_Comparison_with_Lackner_et_al/Data/Annotated Probe Report for All Probes_Log2TPMs.txt",sep="\t")
Med12_TPM <- Med12_TPM[!duplicated(Med12_TPM$Probe),] # Remove duplicate Genes
Med12_TPM <- Med12_TPM[,c(1,13:24)]


# Downregulated wt
Pluri_TPM <- subset(Med12_TPM, Probe %in% Top10down)[,2:13]
rownames(Pluri_TPM) <- subset(Med12_TPM, Probe %in% Top10down)[,1]
# Pluri_TPM <- 2**Pluri_TPM

Pluri_TPM <- data.frame(t(Pluri_TPM))
Pluri_TPM$Rep <- as.numeric(str_sub(rownames(Pluri_TPM), -1,-1))
Pluri_TPM$Condition <- str_sub(rownames(Pluri_TPM), -10,-3)

Pluri_TPM <- melt(Pluri_TPM, id.vars = c("Condition", "Rep"))
colnames(Pluri_TPM) <- c("Condition", "Rep", "Gene", "Expression")

Pluri_TPM$Media <- sub(".*_","",Pluri_TPM$Condition)
Pluri_TPM$Genotype <- sub("_.*","",Pluri_TPM$Condition)
mean_data <- Pluri_TPM %>%
  group_by(Gene, Media, Genotype) %>%
  summarize(mean_value = mean(Expression),
            sd_value = sd(Expression))
mean_data <- data.frame(mean_data)

Pluri_TPM_wide <- reshape(mean_data, direction = "wide", timevar = "Genotype", idvar = c("Gene","Media"))

Pluri_TPM_wide_wide <- reshape(Pluri_TPM_wide, direction = "wide", timevar = c("Media"), idvar = c("Gene"))
Pluri_TPM_wide_wide$slope <- (Pluri_TPM_wide_wide$mean_value.M12KO.2i-Pluri_TPM_wide_wide$mean_value.M12KO.N2)/(Pluri_TPM_wide_wide$mean_value.wt.2i-Pluri_TPM_wide_wide$mean_value.wt.N2)
Slope <- data.frame(Condition = "Most downregulated in wt", Slope = Pluri_TPM_wide_wide$slope)





# Downregulated mutant
Pluri_TPM <- subset(Med12_TPM, Probe %in% Top10down_mutant)[,2:13]
rownames(Pluri_TPM) <- subset(Med12_TPM, Probe %in% Top10down_mutant)[,1]
# Pluri_TPM <- 2**Pluri_TPM

Pluri_TPM <- data.frame(t(Pluri_TPM))
Pluri_TPM$Rep <- as.numeric(str_sub(rownames(Pluri_TPM), -1,-1))
Pluri_TPM$Condition <- str_sub(rownames(Pluri_TPM), -10,-3)

Pluri_TPM <- melt(Pluri_TPM, id.vars = c("Condition", "Rep"))
colnames(Pluri_TPM) <- c("Condition", "Rep", "Gene", "Expression")

Pluri_TPM$Media <- sub(".*_","",Pluri_TPM$Condition)
Pluri_TPM$Genotype <- sub("_.*","",Pluri_TPM$Condition)
mean_data <- Pluri_TPM %>%
  group_by(Gene, Media, Genotype) %>%
  summarize(mean_value = mean(Expression),
            sd_value = sd(Expression))
mean_data <- data.frame(mean_data)

Pluri_TPM_wide <- reshape(mean_data, direction = "wide", timevar = "Genotype", idvar = c("Gene","Media"))

Pluri_TPM_wide_wide <- reshape(Pluri_TPM_wide, direction = "wide", timevar = c("Media"), idvar = c("Gene"))
Pluri_TPM_wide_wide$slope <- (Pluri_TPM_wide_wide$mean_value.M12KO.2i-Pluri_TPM_wide_wide$mean_value.M12KO.N2)/(Pluri_TPM_wide_wide$mean_value.wt.2i-Pluri_TPM_wide_wide$mean_value.wt.N2)
Slope <- rbind(Slope, data.frame(Condition = "Most downregulated in mutant", Slope = Pluri_TPM_wide_wide$slope))



# Pluripotency markers
Pluri_TPM <- subset(Med12_TPM, Probe %in% pluripotency_markers)[,2:13]
rownames(Pluri_TPM) <- subset(Med12_TPM, Probe %in% pluripotency_markers)[,1]
# Pluri_TPM <- 2**Pluri_TPM

Pluri_TPM <- data.frame(t(Pluri_TPM))
Pluri_TPM$Rep <- as.numeric(str_sub(rownames(Pluri_TPM), -1,-1))
Pluri_TPM$Condition <- str_sub(rownames(Pluri_TPM), -10,-3)

Pluri_TPM <- melt(Pluri_TPM, id.vars = c("Condition", "Rep"))
colnames(Pluri_TPM) <- c("Condition", "Rep", "Gene", "Expression")

Pluri_TPM$Media <- sub(".*_","",Pluri_TPM$Condition)
Pluri_TPM$Genotype <- sub("_.*","",Pluri_TPM$Condition)
mean_data <- Pluri_TPM %>%
  group_by(Gene, Media, Genotype) %>%
  summarize(mean_value = mean(Expression),
            sd_value = sd(Expression))
mean_data <- data.frame(mean_data)

Pluri_TPM_wide <- reshape(mean_data, direction = "wide", timevar = "Genotype", idvar = c("Gene","Media"))
ggplot(Pluri_TPM_wide, aes(x=mean_value.wt, y= mean_value.M12KO, color = Gene, group = Gene, shape = Media, fill = Media)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +
  geom_path() + 
  geom_errorbar(aes(ymin = mean_value.M12KO - sd_value.M12KO, ymax = mean_value.M12KO + sd_value.M12KO), width = 0.2, color = "grey") +
  geom_errorbar(aes(xmin = mean_value.wt - sd_value.wt, xmax = mean_value.wt + sd_value.wt), width = 0.2, color = "grey") +
  geom_point(size = 1) +
  # scale_x_continuous(limits = c(-0.1,9.3), expand = c(0,0)) +
  # scale_y_continuous(limits = c(-0.1,9.3), expand = c(0,0)) +
  scale_shape_manual(values = c(21,20)) +
  scale_fill_manual(values = c("white","NA"))+
  xlab("Expression (log2-TPM) in wild type") +
  ylab("Expression (log2-TPM) in Med12 mutant") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        text = element_text(color = 'black'),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        # legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA),
        strip.background = element_blank()) +
  force_panelsizes(rows = unit(6, "cm"),
                   cols = unit(6, "cm"))
# ggsave("/Users/fernkorn/Documents/RNA-seq_data/Own_Data/bulk_RNA-seq/Med12KO_24h_Comparison/Lackner_like_analysis/RNAseq_data/20240109_Expression_of_Top100down_in_mutant_Dotplot.pdf", height = 5, width = 5)



# Upregulated mutant
Pluri_TPM <- subset(Med12_TPM, Probe %in% Top10up_mutant)[,2:13]
rownames(Pluri_TPM) <- subset(Med12_TPM, Probe %in% Top10up_mutant)[,1]
# Pluri_TPM <- 2**Pluri_TPM

Pluri_TPM <- data.frame(t(Pluri_TPM))
Pluri_TPM$Rep <- as.numeric(str_sub(rownames(Pluri_TPM), -1,-1))
Pluri_TPM$Condition <- str_sub(rownames(Pluri_TPM), -10,-3)

Pluri_TPM <- melt(Pluri_TPM, id.vars = c("Condition", "Rep"))
colnames(Pluri_TPM) <- c("Condition", "Rep", "Gene", "Expression")

Pluri_TPM$Media <- sub(".*_","",Pluri_TPM$Condition)
Pluri_TPM$Genotype <- sub("_.*","",Pluri_TPM$Condition)
mean_data <- Pluri_TPM %>%
  group_by(Gene, Media, Genotype) %>%
  summarize(mean_value = mean(Expression),
            sd_value = sd(Expression))
mean_data <- data.frame(mean_data)

Pluri_TPM_wide <- reshape(mean_data, direction = "wide", timevar = "Genotype", idvar = c("Gene","Media"))

Pluri_TPM_wide_wide <- reshape(Pluri_TPM_wide, direction = "wide", timevar = c("Media"), idvar = c("Gene"))
Pluri_TPM_wide_wide$slope <- (Pluri_TPM_wide_wide$mean_value.M12KO.2i-Pluri_TPM_wide_wide$mean_value.M12KO.N2)/(Pluri_TPM_wide_wide$mean_value.wt.2i-Pluri_TPM_wide_wide$mean_value.wt.N2)
Slope <- rbind(Slope, data.frame(Condition = "Most upregulated in mutant", Slope = Pluri_TPM_wide_wide$slope))

# Upregulated wildtype
Pluri_TPM <- subset(Med12_TPM, Probe %in% Top10up)[,2:13]
rownames(Pluri_TPM) <- subset(Med12_TPM, Probe %in% Top10up)[,1]
# Pluri_TPM <- 2**Pluri_TPM

Pluri_TPM <- data.frame(t(Pluri_TPM))
Pluri_TPM$Rep <- as.numeric(str_sub(rownames(Pluri_TPM), -1,-1))
Pluri_TPM$Condition <- str_sub(rownames(Pluri_TPM), -10,-3)

Pluri_TPM <- melt(Pluri_TPM, id.vars = c("Condition", "Rep"))
colnames(Pluri_TPM) <- c("Condition", "Rep", "Gene", "Expression")

Pluri_TPM$Media <- sub(".*_","",Pluri_TPM$Condition)
Pluri_TPM$Genotype <- sub("_.*","",Pluri_TPM$Condition)
mean_data <- Pluri_TPM %>%
  group_by(Gene, Media, Genotype) %>%
  summarize(mean_value = mean(Expression),
            sd_value = sd(Expression))
mean_data <- data.frame(mean_data)

Pluri_TPM_wide <- reshape(mean_data, direction = "wide", timevar = "Genotype", idvar = c("Gene","Media"))

Pluri_TPM_wide_wide <- reshape(Pluri_TPM_wide, direction = "wide", timevar = c("Media"), idvar = c("Gene"))
Pluri_TPM_wide_wide$slope <- (Pluri_TPM_wide_wide$mean_value.M12KO.2i-Pluri_TPM_wide_wide$mean_value.M12KO.N2)/(Pluri_TPM_wide_wide$mean_value.wt.2i-Pluri_TPM_wide_wide$mean_value.wt.N2)
Slope <- rbind(Slope, data.frame(Condition = "Most upregulated in wt", Slope = Pluri_TPM_wide_wide$slope))


Slope$Condition <- factor(Slope$Condition, levels = c("Most downregulated in wt","Most downregulated in mutant","Most upregulated in wt","Most upregulated in mutant"),
                          ordered = TRUE)








# Plot Slope for Top 100 genes (Fig 5E, left panel)
Slope_down <- subset(Slope, Condition %in% c("Most downregulated in wt","Most downregulated in mutant"))

ggplot(Slope_down, aes(x = Condition, y = Slope)) + 
  geom_violin(linewidth = 0.2, fill = "darkgrey") +
  geom_boxplot(position = position_dodge(width = 0.9),width = 0.2,outlier.size = 0.3,linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  ylab("downregulation slope") +
  xlab("condition") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        text = element_text(color = 'black'),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, vjust = 1),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        # legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA),
        strip.background = element_blank()) +
  force_panelsizes(rows = unit(3.5, "cm"),
                   cols = unit(3, "cm"))
ggsave("./20240131_Distribution_Slopes_Top100_Downregualted_wide.pdf", height = 4, width = 3)


# Plot Slope for Top 100 genes (Fig 5E, right panel)
Slope_up <- subset(Slope, Condition %in% c("Most upregulated in wt","Most upregulated in mutant"))

ggplot(Slope_up, aes(x = Condition, y = Slope)) + 
  geom_violin(linewidth = 0.2, fill = "darkgrey") +
  geom_boxplot(position = position_dodge(width = 0.9),width = 0.2,outlier.size = 0.3,linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  ylab("upregulation slope") +
  xlab("condition") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        text = element_text(color = 'black'),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, vjust = 1),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        # legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA),
        strip.background = element_blank()) +
  force_panelsizes(rows = unit(3.5, "cm"),
                   cols = unit(3, "cm"))
ggsave("./20240131_Distribution_Slopes_Top100_Upregualted_wide.pdf", height = 4, width = 3)


