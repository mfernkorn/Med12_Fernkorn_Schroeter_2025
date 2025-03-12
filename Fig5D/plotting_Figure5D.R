# Plot Figure 5D
#load libraries
library(reshape2)
library(ggplot2)
library(stringr)
library(ggpubr)
library(ggh4x)
library(dplyr)

pluripotency_markers <- c("Klf4", "Nanog", "Esrrb", "Tbx3", "Tfcp2l1","Prdm14", "Zfp42")


Med12_TPM <- read.csv(file = "./Annotated Probe Report for All Probes_Log2TPMs.txt",sep="\t")
Med12_TPM <- Med12_TPM[!duplicated(Med12_TPM$Probe),] # Remove duplicate Genes
Med12_TPM <- Med12_TPM[,c(1,13:24)]

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

# Plot Figure 5D
ggplot(Pluri_TPM_wide, aes(x=mean_value.wt, y= mean_value.M12KO, color = Gene, group = Gene, shape = Media, fill = Media)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +
  geom_path() + 
  geom_errorbar(aes(ymin = mean_value.M12KO - sd_value.M12KO, ymax = mean_value.M12KO + sd_value.M12KO), width = 0.2, color = "grey") +
  geom_errorbar(aes(xmin = mean_value.wt - sd_value.wt, xmax = mean_value.wt + sd_value.wt), width = 0.2, color = "grey") +
  geom_point(size = 1) +
  scale_x_continuous(limits = c(-0.1,9.3), expand = c(0,0)) +
  scale_y_continuous(limits = c(-0.1,9.3), expand = c(0,0)) +
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
# ggsave("./20231110_Expression_of_pluripotency_markers_Dotplot.pdf", height = 5, width = 5)
