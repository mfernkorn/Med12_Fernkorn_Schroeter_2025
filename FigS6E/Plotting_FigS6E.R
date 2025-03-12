# PLot FigS6E
# load libraries
library(ggplot2)
library(stringr)
library(ggh4x) # for force_panelsize option
library(reshape2)

# Load data from csvs, produced with Segmentation_and_Measurement.ijm
# necessary to type in condition names manually, take care about filename
filenames <- list.files(path = "./Csvs/Rep1/")
Condition <- c(rep("no FGF", 3), rep("10 ng/µL FGF", 3),
               rep("no FGF", 3), rep("10 ng/µL FGF", 3))

# Loop over all files
dat_1 <- data.frame()
for (i in 1:length(filenames)){
  temp_dat <- read.csv(paste0("./Csvs/Rep1/", filenames[i]))[,2:6]
  temp_dat$Channel <- str_sub(filenames[i], start = -8, end = -5)
  temp_dat$Genotype <- str_sub(filenames[i], start = 13, end = 15)
  temp_dat$Condition <- Condition[i]
  temp_dat$Cell <- 1:length(temp_dat$Condition)
  dat_1 <- rbind(dat_1, temp_dat)
}
dat_1$Rep <- "#1"

filenames <- list.files(path = "./Csvs/Rep2/")
# Loop over all files
dat_2 <- data.frame()
for (i in 1:length(filenames)){
  temp_dat <- read.csv(paste0("./Csvs/Rep2/", filenames[i]))[,2:6]
  temp_dat$Channel <- str_sub(filenames[i], start = -8, end = -5)
  temp_dat$Genotype <- str_sub(filenames[i], start = 13, end = 15)
  temp_dat$Condition <- Condition[i]
  temp_dat$Cell <- 1:length(temp_dat$Condition)
  dat_2 <- rbind(dat_2, temp_dat)
}
dat_2$Rep <- "#2"

dat <- rbind(dat_1, dat_2)
# rename levels of Genotype and Channel to more meaningful names
dat$Genotype <- as.factor(dat$Genotype)
levels(dat$Genotype) <- c('Med12 KO', 'wt')
dat$Channel <- as.factor(dat$Channel)
levels(dat$Channel) <- c('Nuclei', 'Nanog', 'Sox17')

# Filter for minimal area of a cell
hist(dat$Area)
dat <- subset(dat, Area >=40)

#Rearrange columns and only use IntDens
dat_reshape <- reshape(dat[,c(2,6,7,8,9,10)],
                       idvar = c("Genotype", "Condition", "Cell", "Rep"),
                       timevar = "Channel",
                       direction = "wide")

dat_reshape <- subset(dat_reshape, Mean.Nanog >=10)
dat_reshape <- subset(dat_reshape, Mean.Sox17 >=10)

SOX17threshold <- 40 
NANOGthreshold <- 35

# Plot intensities between channels for one condition for both genotypes
ggplot(dat_reshape, aes(x=Mean.Nanog, y=Mean.Sox17, color = Rep)) +
  geom_point(size = 0.3) +
  geom_hline(yintercept = SOX17threshold) + geom_vline(xintercept = NANOGthreshold) +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_manual(values=c("#B6D0DD","#164876")) +
  annotation_logticks(sides = "lb", outside = TRUE) +
  facet_grid(Genotype~Condition) + 
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
ggsave("./Plots/20231115_Nanog_Sox17_Intensities_Sorting_Both_Reps.pdf", width = 30, height = 6)

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

dat_frequencies <- table(paste(dat_reshape$Description, dat_reshape$Rep, sep = "__"), dat_reshape$class)
dat_frequencies <- dat_frequencies/rowSums(dat_frequencies)
dat_frequencies <- melt(dat_frequencies)
colnames(dat_frequencies) <- c("Description", "Celltype", "Fraction")

dat_frequencies$Celltype <- factor(dat_frequencies$Celltype, 
                                   levels = c("double negative","SOX17-; NANOG+",
                                              "double positive","SOX17+; NANOG-"), 
                                   ordered = TRUE)
dat_frequencies$Genotype <- factor(dat_frequencies$Genotype, levels = c("wt", "Med12 KO"), ordered = TRUE)
dat_frequencies$Condition <- factor(dat_frequencies$Condition, levels = c("no FGF", "10 ng/µL FGF"), ordered = TRUE)
dat_frequencies <- dat_frequencies[order(dat_frequencies$Celltype,decreasing=T),]
                
ggplot(dat_frequencies, aes(x=Rep, y=Fraction, fill= Celltype, group = Rep)) + 
  geom_bar(position="stack",stat = "identity", width = 1.2) +
  # scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values=c("#8BC5EB","#00A76E","#D1BB2B","#D86FAA")) +
  facet_wrap(vars(Genotype, Condition, Rep),scales = "free_x", nrow = 1, ncol = 8) +
  theme(panel.spacing = unit(0, "lines"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.48),
        panel.border = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        # legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA),
        strip.text = element_text(color = "black"),
        panel.spacing.x = unit(c(0,5,0,10,0,5,0), "pt")) +
  force_panelsizes(rows = unit(3, "cm"),
                   cols = unit(0.3, "cm"))
ggsave("./Plots/20231115_Celltype_fractions_Sorting_Both_Reps.pdf", width = 12, height = 6)


