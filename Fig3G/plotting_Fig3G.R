# Create Figure Panel Fig3G
# Load libraries
library(DESeq2)
library(ggrepel)
library(ggplot2)
library(ggh4x) # for force_panelsize option

# Load TPMs for PCA 
cts <- read.csv(file = "Annotated Probe Report for Value above 50_TPMs.txt",sep="\t")
# Make Probe names unique by including Positional information
rownames(cts) <- paste(cts$Probe, cts$Chromosome, cts$Start, cts$End, sep = "_") 
cts <- as.matrix(cts[,13:28])

# Calculate PCA contributions and plot
PCAs <- prcomp(t(cts))
PCAs$sdev^2 / sum(PCAs$sdev^2)

PCAs <- data.frame(PCAs$x)
PCAs$Replicate <- c(rep(c("1","2","3"), 4),1,1,1,1)
PCAs$Condition <- c(rep("wt", 3), rep("wt + FGF", 3), 
                    rep("Med12-mutant", 3), rep("Med12-mutant + FGF", 3),
                    "Med12-mutant", "Med12-mutant + FGF", "Med12-mutant", "Med12-mutant + FGF")

ggplot(PCAs, aes(x=PC1, y=PC2, color=Condition, shape=Replicate)) +
  geom_point(size=2) +
  xlab("PC1 (33.2%)") +
  ylab("PC2 (12.8%)") +
  scale_color_manual(values = c( "#E98624", "#F7D1A5","#2E6A67", "#91E2DA"))+
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
  force_panelsizes(rows = unit(4, "cm"),
                   cols = unit(4, "cm"))
ggsave("./20231102_PCA_Plot_TPMs_PC1_vs_PC2.pdf", width = 4, height = 4)

ggplot(PCAs, aes(x=PC1, y=PC3, color=Condition, shape=Replicate)) +
  geom_point(size=2) +
  xlab("PC1 (33.2%)") +
  ylab("PC3 (9.3%)") +
  scale_color_manual(values = c( "#E98624", "#F7D1A5","#2E6A67", "#91E2DA"))+
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
  force_panelsizes(rows = unit(4, "cm"),
                   cols = unit(4, "cm"))
ggsave("./20231102_PCA_Plot_TPMs_PC1_vs_PC3.pdf", width = 4, height = 4)
