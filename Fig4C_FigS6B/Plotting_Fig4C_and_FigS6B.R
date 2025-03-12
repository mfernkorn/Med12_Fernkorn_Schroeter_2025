# Plot Figure panels 4C and S6B
# load libraries
library(ggplot2)
library(ggh4x) # for force_panelsize option
library(dplyr)

# Load Fraction data from all single replicates
Rep0 <- readRDS("./Rep0.Rds")
Rep1 <- readRDS("./Rep1.Rds")
Rep2 <- readRDS("./Rep2.Rds")

# Merge into one single data frame all meta data
df_list <- list(Rep0, Rep1, Rep2)      
all_reps <- Reduce(function(x, y) merge(x, y, all=TRUE, by = c("Description", "Celltype", "Genotype", "Condition")), df_list)
all_reps <- all_reps[rev(order(all_reps[,2])), ]
# Calculate Mean of Fractions over replicates and Standard Error
all_reps$Mean_Fraction <- rowMeans(all_reps[,5:7], na.rm = TRUE)
all_reps$sd_Fraction <- apply(all_reps[,5:7], 1, sd, na.rm = TRUE)

all_reps <- all_reps %>% group_by(Genotype, Condition, Description) %>% 
  mutate(SDpos = cumsum(Mean_Fraction))

all_reps$Celltype <- factor(all_reps$Celltype, levels = c("double negative", "SOX17-; NANOG+","double positive","SOX17+; NANOG-"), ordered = TRUE)
all_reps <- subset(all_reps, !(Condition %in% c("nonreseeded_nodox_4h","nonreseeded_nodox_4.5h")))
all_reps$Condition <- factor(all_reps$Condition)
levels(all_reps$Condition) <- c("no FGF, 4.5h dox", "no FGF, 4h dox","no FGF, 8h dox","10 ng/mL FGF, 4.5h dox","10 ng/mL FGF, 4h dox","10 ng/mL FGF, 8h dox","no FGF, no dox")
all_reps$Condition <- factor(all_reps$Condition, levels = c("no FGF, no dox","no FGF, 4h dox","no FGF, 4.5h dox","no FGF, 8h dox","10 ng/mL FGF, 4h dox","10 ng/mL FGF, 4.5h dox","10 ng/mL FGF, 8h dox"), ordered = TRUE)
all_reps$Genotype <- factor(all_reps$Genotype, levels = c("wt", "Med12 KO"), ordered = TRUE)

# Focus on samples with additional FGF and with 8h induction (for Panel Fig4C)
all_reps_Fig5E <- subset(all_reps, (Condition  == "no FGF, 8h dox"))

ggplot(all_reps_Fig5E, aes(fill=Celltype, y=Mean_Fraction, x=Genotype)) + 
  geom_bar(position="stack", stat="identity", colour="black", linewidth = 0.2) +
  theme_classic() + 
  labs(x = "Genotype", y = "Fraction") +
  scale_fill_manual(values=c("#8BC5EB","#00A76E","#D1BB2B","#D86FAA")) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black"), 
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        # legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA)) +
  geom_errorbar(aes(ymin=SDpos-sd_Fraction, ymax=SDpos), width=.4,
                position=position_dodge(0), linewidth = 0.24) +
  force_panelsizes(rows = unit(3, "cm"),
                   cols = unit(1, "cm"))
#ggsave("./20231130_Compare_wt(4h)_vs_Med12KO(8h).pdf", width = 4, height = 4)

# Plot sll samples (with and without FGF) for 4h and 8h induction (for Panel FigS6B)
all_reps_plot <- subset(all_reps, (Condition %in% c("no FGF, 8h dox","no FGF, 4h dox","10 ng/mL FGF, 4h dox","10 ng/mL FGF, 8h dox")))

ggplot(all_reps_plot, aes(fill=Celltype, y=Mean_Fraction, x=Condition, width=.8)) + 
  geom_bar(position="stack", stat="identity", colour="black", linewidth = 0.2) +
  theme_classic() + 
  labs(x = "Condition", y = "Fraction") +
  scale_fill_manual(values=c("#8BC5EB","#00A76E","#D1BB2B","#D86FAA")) +
  facet_grid(~Genotype, scales = "free",space = "free") + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black", angle = 90), 
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        # legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA)) +
  geom_errorbar(aes(ymin=SDpos-sd_Fraction, ymax=SDpos), width=.4, linewidth=0.2,
                position=position_dodge(0))
#ggsave("./20231115_Compare_8hvs4h_Rep0to2_CS.pdf", width = 5, height = 3)

# Print n's for each condition
all_reps_plot %>% group_by(Description, Celltype) %>%
  summarise(mean = mean())
