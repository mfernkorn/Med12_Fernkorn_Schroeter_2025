# Generate Figure panels 1B, 1C, 2A, S2, S3
# Load libraries
library(ggplot2)
library(reshape2)
library(ggrepel)
library(ggh4x)

# Load data, output from MAGeCK, sgRNA level data
L11vsNons1 <- read.csv("./L11vsNons1_test.sgRNA_summary.txt", sep = "\t",)
L51vsNons1 <- read.csv("./L51vsNons1_test.sgrna_summary.txt", sep = "\t",)
L12vsNons2 <- read.csv("./L12vsNons2_test.sgrna_summary.txt", sep = "\t",)
L52vsNons2 <- read.csv("./L52vsNons2_test.sgrna_summary.txt", sep = "\t",)

H11vsNons1 <- read.csv("./H11vsNons1_test.sgrna_summary.txt", sep = "\t",)
H12vsNons2 <- read.csv("./H12vsNons2_test.sgrna_summary.txt", sep = "\t",)
H51vsNons1 <- read.csv("./H51vsNons1_test.sgrna_summary.txt", sep = "\t",)
H52vsNons2 <- read.csv("./H52vsNons2_test.sgrna_summary.txt", sep = "\t",)

# Plot LFC distribution (Figure 1B, S2A-C, S3A-D)
plot_LFC_distribution <- function(dataset, filename = quote(dataset)){
  dataset_plot <- dataset[order(-dataset$LFC),]
  dataset_plot$order <- c(1:length(dataset$LFC))
  dataset_plot$group <- "targeting"
  dataset_plot[which(substring(dataset_plot$Gene,1,4) == "ntc_"),17] <- "non-targeting"
  dataset_plot <- dataset_plot[with(dataset_plot, order(dataset_plot$group, decreasing = TRUE)), ]
  ggplot(data=dataset_plot, aes(x=order, y=LFC)) + 
    geom_point(aes(color=group, shape=group), size = 0.7) + 
    scale_color_manual(values = c("#9CBCD3","#144976")) +
    scale_shape_manual(values = c(20,19)) +
    scale_x_continuous(trans='log10') +
    xlab(label = "Number of gRNAs") + ylab(label = "Log2-Foldchange") +
    theme_bw(base_size = 10) + 
    theme(aspect.ratio = 1, text = element_text(color = 'black'), 
          axis.text = element_text(colour = "black"), 
          panel.border = element_rect(color = "black", linewidth = 1),
          # panel.background = element_rect(fill = NA),
          plot.background = element_rect(fill = NA, linewidth = 0),
          legend.background = element_rect(fill = NA, linewidth = 0),
          legend.key = element_rect(fill = NA)) +
    force_panelsizes(rows = unit(5, "cm"),
                     cols = unit(4, "cm"))
  ggsave(paste0(filename,".png"), height = 12, width = 12, dpi = 1000, units = "cm")
  ggsave(paste0(filename,".pdf"), height = 12, width = 12, dpi = 1000, units = "cm")
}

plot_LFC_distribution(L11vsNons1, filename = "./20230707_Fig1B")
plot_LFC_distribution(L51vsNons1, filename = "./20230707_FigS2A")
plot_LFC_distribution(L12vsNons2, filename = "./20230707_Fig2B")
plot_LFC_distribution(L52vsNons2, filename = "./20230707_Fig2C")

plot_LFC_distribution(H11vsNons1, filename = "./20230707_FigS3A")
plot_LFC_distribution(H51vsNons1, filename = "./20230707_FigS3B")
plot_LFC_distribution(H12vsNons2, filename = "./20230707_FigS3C")
plot_LFC_distribution(H52vsNons2, filename = "./20230707_FigS3D")

# Load data, output from MAGeCK, gene level data
L11vsNons1 <- read.csv("./L11vsNons1/L11vsNons1_test.gene_summary.txt", sep = "\t",)
L51vsNons1 <- read.csv("./L51vsNons1/L51vsNons1_test.gene_summary.txt", sep = "\t",)
L12vsNons2 <- read.csv("./L12vsNons2/L12vsNons2_test.gene_summary.txt", sep = "\t",)
L52vsNons2 <- read.csv("./L52vsNons2/L52vsNons2_test.gene_summary.txt", sep = "\t",)

H11vsNons1 <- read.csv("./H11vsNons1/H11vsNons1_test.gene_summary.txt", sep = "\t",)
H12vsNons2 <- read.csv("./H12vsNons2/H12vsNons2_test.gene_summary.txt", sep = "\t",)
H51vsNons1 <- read.csv("./H51vsNons1/H51vsNons1_test.gene_summary.txt", sep = "\t",)
H52vsNons2 <- read.csv("./H52vsNons2/H52vsNons2_test.gene_summary.txt", sep = "\t",)

# RRA plot for genes (Figure 1C, S2D-F, 2A, S3E-G)
plot_RRA <- function(dataset, filename = quote(dataset)){
  dataset_plot <- dataset[order(dataset$pos.score),]
  dataset_plot$order <- c(1:length(dataset$pos.score))
  threshold_0_05 <- (dataset_plot$pos.score[which(dataset_plot$pos.fdr > 0.05)[1]-1] + dataset_plot$pos.score[which(dataset_plot$pos.fdr > 0.05)[2]-1])/2
  threshold_0_2 <- (dataset_plot$pos.score[which(dataset_plot$pos.fdr > 0.2)[1]-1] + dataset_plot$pos.score[which(dataset_plot$pos.fdr > 0.2)[2]-1])/2
  
  ggplot(data=dataset_plot, aes(x=order, y=pos.score)) +
    annotate("rect", xmin=-Inf,xmax=Inf,ymin=threshold_0_05,ymax=0, alpha=0.2, fill="#166420") +
    annotate("rect", xmin=-Inf,xmax=Inf,ymin=threshold_0_2,ymax=threshold_0_05, alpha=0.2, fill="#324147") +
    geom_point(colour="#144976", size = 0.7) + scale_y_continuous(trans='log10') + 
    xlab(label = "Number of Genes") + ylab(label = "RRA-score") + 
    geom_text_repel(data=subset(dataset_plot, pos.fdr<0.05), aes(label=subset(dataset_plot, pos.fdr<0.05)[,1]),
                    max.overlaps = 30, box.padding = 0.8, size = 8*0.36) +
    annotate("text", x=20000,y=threshold_0_05,label = "FDR < 0.05", vjust = 1.4, hjust = "right", size = 8*0.36) +
    annotate("text", x=20000,y=threshold_0_2,label = "FDR < 0.2", vjust = 1.4, hjust = "right", size = 8*0.36) + 
    theme_bw(base_size = 10) + 
    theme(text = element_text(color = 'black'), 
          axis.text = element_text(colour = "black"), 
          panel.border = element_rect(color = "black", linewidth = 1),
          # panel.background = element_rect(fill = NA),
          plot.background = element_rect(fill = NA, linewidth = 0),
          legend.background = element_rect(fill = NA, linewidth = 0),
          legend.key = element_rect(fill = NA)) +
    force_panelsizes(rows = unit(5, "cm"),
                     cols = unit(4, "cm"))
  ggsave(paste0(filename,".png"), height = 12, width = 12, dpi = 1000, units = "cm")
  ggsave(paste0(filename,".pdf"), height = 12, width = 12, dpi = 1000, units = "cm")
}

plot_RRA(L11vsNons1, filename = "./20230707_Fig1C")
plot_RRA(L51vsNons1, filename = "./20230707_FigS2D")
plot_RRA(L12vsNons2, filename = "./20230707_FigS2E")
plot_RRA(L52vsNons2, filename = "./20230707_FigS2F")

plot_RRA(H11vsNons1, filename = "./20230707_Fig2A")
plot_RRA(H12vsNons2, filename = "./20230707_FigS3E")
plot_RRA(H51vsNons1, filename = "./20230707_FigS3F")
plot_RRA(H52vsNons2, filename = "./20230707_FigS3G")

