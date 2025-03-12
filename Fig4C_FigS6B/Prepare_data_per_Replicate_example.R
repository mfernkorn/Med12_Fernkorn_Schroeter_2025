# Example for exporting fractions from one replicate in an Rds file
# load libraries
library(ggplot2)
library(stringr)
library(ggh4x) # for force_panelsize option
library(reshape2)

# Load data from csvs, produced with Segmentation_and_Measurement.ijm
filenames <- list.files(path = "./Csvs/")
Condition <- c(rep("nonreseeded_0FGF8hdox", 3), rep("nonreseeded_10FGF8hdox", 3),
               rep("nonreseeded_nodox_8h", 3), 
               rep("nonreseeded_0FGF4.5hdox", 3), rep("nonreseeded_0FGF4hdox", 3), rep("nonreseeded_0FGF8hdox", 3),
               rep("nonreseeded_10FGF4.5hdox", 3),rep("nonreseeded_10FGF4hdox", 3),rep("nonreseeded_10FGF8hdox", 3),
               rep("nonreseeded_nodox_8h", 3),rep("nonreseeded_nodox_4.5h", 3),rep("nonreseeded_nodox_4h", 3))
               

# Loop over all files
dat <- data.frame()
for (i in 1:length(filenames)){
  temp_dat <- read.csv(paste0("./Csvs/", filenames[i]))[,2:6]
  temp_dat$Channel <- str_sub(filenames[i], start = -8, end = -5)
  temp_dat$Genotype <- str_sub(filenames[i], start = 13, end = 15)
  temp_dat$Condition <- Condition[i]
  temp_dat$Cell <- 1:length(temp_dat$Condition)
  dat <- rbind(dat, temp_dat)
}

# rename levels of Genotype and Channel to more meaningful names
dat$Genotype <- as.factor(dat$Genotype)
levels(dat$Genotype) <- c('Med12 KO', 'wt')
dat$Channel <- as.factor(dat$Channel)
levels(dat$Channel) <- c('Nuclei', 'Nanog', 'Sox17')

# Filter for minimal area of a cell
dat <- subset(dat, Area >=30)
#Rearrange columns and only use IntDens
dat_reshape <- reshape(dat[,c(2,6,7,8,9)],
                       idvar = c("Genotype", "Condition", "Cell"),
                       timevar = "Channel",
                       direction = "wide")
SOX17threshold <- 25 
NANOGthreshold <- 30

# Plot intensities between channels for one condition for both genotypes
ggplot(dat_reshape, aes(x=Mean.Nanog, y=Mean.Sox17, color = Genotype)) +
  geom_point(size = 0.3) +
  geom_hline(yintercept = SOX17threshold) + geom_vline(xintercept = NANOGthreshold) +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_manual(values=c("#B6D0DD","#164876")) +
  annotation_logticks(sides = "lb", outside = TRUE) +
  facet_wrap(vars(Condition), nrow = 1) + 
  ylab("SOX17 Fluorescence") + 
  xlab("NANOG  Fluorescence") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black"), 
        panel.border = element_rect(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        # legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA)) +
  force_panelsizes(rows = unit(4, "cm"),
                   cols = unit(4, "cm"))

# Get frequencies of cell types
dat_reshape$Description <- factor(paste(dat_reshape$Genotype, dat_reshape$Condition, sep = "__"))
for (i in 1:length(dat_reshape$Genotype)){
  if (dat_reshape$Mean.Sox17[i] > SOX17threshold && dat_reshape$Mean.Nanog[i] > NANOGthreshold){
    dat_reshape$class[i] <- "double positive"
  } else if (dat_reshape$Mean.Sox17[i] < SOX17threshold && dat_reshape$Mean.Nanog[i] > NANOGthreshold){
    dat_reshape$class[i] <- "SOX17-; NANOG+"
  } else if (dat_reshape$Mean.Sox17[i] > SOX17threshold && dat_reshape$Mean.Nanog[i] < NANOGthreshold){
    dat_reshape$class[i] <- "SOX17+; NANOG-"
  } else {
    dat_reshape$class[i] <- "double negative"
  }
}

dat_frequencies <- table(dat_reshape$Description, dat_reshape$class)
dat_frequencies <- dat_frequencies/rowSums(dat_frequencies)
dat_frequencies <- melt(dat_frequencies)
colnames(dat_frequencies) <- c("Description", "Celltype", "Fraction")
dat_frequencies[c('Genotype', 'Condition')] <- str_split_fixed(dat_frequencies$Description, '__', 2)
dat_frequencies$Celltype <- factor(dat_frequencies$Celltype, 
                                   levels = c("double negative","SOX17-; NANOG+",
                                              "double positive","SOX17+; NANOG-"), 
                                   ordered = TRUE)

# Plot frequencies in barplots similar to Raina et al
ggplot(dat_frequencies, aes(fill=Celltype, y=Fraction, x=Condition, width=.8)) + 
  geom_bar(position="stack", stat="identity", colour="black") +
  theme_classic() + 
  labs(x = "Condition", y = "Fraction") +
  scale_fill_manual(values=c("#8BC5EB","#00A76E","#D1BB2B","#D86FAA")) +
  facet_wrap(vars(Genotype), nrow = 1, scales = "free_x") + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black", angle = 90), 
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        # legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA)) +
  force_panelsizes(rows = unit(8, "cm"),
                   cols = unit(10, "cm"))

# Save frequencies of cell types as Rds for comparison of multiple replicates
saveRDS(dat_frequencies, "/Users/fernkorn/Documents/Immunostaining_Analysis/Comparison_8vs4h_PrE_Epi/Rep2.Rds")

