# Plot FigS4H
# load libraries
library(dplyr)
library(stringr)
library(cowplot)
library(ggplot2)
library(ggh4x)
library(ggpubr)

# Function to read in any replicate
# countpath shell be the path to the folder with all count files
# nanogstring shell be the string in the filename identifying the Nanog staining (others will be called Spry4)
# repname should could contain an identifier for the sample
loadCounts <- function(countpath, nanogstring, repname){
  filenames <- list.files(path = countpath)
  
  # Loop over all files to read in counts
  # there is a separate file for each cell within each position and channel
  # this lists spots, and the number of spots per file corresponds to the mRNA counts per file
  dat <- data.frame()
  for (i in 1:length(filenames)){
    # Predefine meta data columns, based on the file names. Channels might have to be modified per replicate individually
    if (grepl("2i", filenames[i])){
      Media <- "2i" 
    } else {
      Media <- "N2B27" 
    }
    if (grepl("M12|Med12KO", filenames[i])){
      Genotype <- "Med12 mutant" 
    } else {
      Genotype <- "wild type" 
    }
    if (grepl(nanogstring, filenames[i])){
      Channel <- "Nanog" 
    } else {
      Channel <- "Spry4" 
    }
    if (grepl("negative_control", filenames[i])){
      Sampletype <- "negative control" 
    } else {
      Sampletype <- "sample" 
    }
    Position <- str_extract(filenames[i], "(?<=Pos)\\d+")
    Cell <- str_extract(filenames[i], "\\d+(?=\\.csv)")
    
    # Test if file or cell does contain spots. If not add one line per cell with 0 as spot number
    if (file.info(paste0(countpath, filenames[i]))$size == 0) {
      temp_dat <- data.frame(x = NA, y = NA, z = NA, intensity = NA, media = Media, genotype = Genotype, channel = Channel, position = Position, cell = Cell, sampletype = Sampletype, spot = 0)
    } else {
      temp_dat <- read.csv(paste0(countpath, filenames[i]))[,c(1,2,3,6)]
      temp_dat$media <- Media
      temp_dat$genotype <- Genotype
      temp_dat$channel <- Channel
      temp_dat$position <- Position
      temp_dat$cell <- Cell
      temp_dat$sampletype <- Sampletype
      temp_dat$spot <- 1:length(temp_dat[,1])
    }
    # Merge cell to exisiting cells
    dat <- rbind(dat, temp_dat)
  }
  dat$replicate <- repname
  return(dat)
}

# Load data for all replicates
Rep1_countpath = "./Analysis_Rep1/Counts_new/"
spot_counts_rep1 <- loadCounts(countpath = Rep1_countpath, nanogstring = "GFP", repname = "Rep1")
Rep2_countpath = "./Analysis_Rep2/Counts_new/"
spot_counts_rep2 <- loadCounts(countpath = Rep2_countpath, nanogstring = "A647", repname = "Rep2")
Rep3_countpath = "./Analysis_Rep3/Counts/"
spot_counts_rep3 <- loadCounts(countpath = Rep3_countpath, nanogstring = "A647", repname = "Rep3")
Rep4_countpath = "./Analysis_Rep4/Counts/"
spot_counts_rep4 <- loadCounts(countpath = Rep4_countpath, nanogstring = "GFP", repname = "Rep4")

# Combine all replicates to a single dataframe
spot_counts <- rbind(spot_counts_rep1, spot_counts_rep2, spot_counts_rep3, spot_counts_rep4)

# Summarize numbers of spots for each cell by taking the maximum of the spot column per cell (by group)
spot_counts_summarized <- spot_counts %>% group_by(media, genotype, channel, position, cell, sampletype, replicate) %>%
  summarise(spot_count = max(spot))

# Manually delete positions
# In this case, in rep1 in wild type at position 36 and 37 in Spry4 channel in N2B27 is for unknown reason overexposed
spot_counts_summarized<- subset(spot_counts_summarized, !(media == "N2B27" & channel == "Spry4" & genotype == "wild type" & replicate == "Rep1" & (position == 36 | position == 37)))

# Display sample numbers per replicate and condition
table(paste(spot_counts_summarized$media, spot_counts_summarized$genotype, spot_counts_summarized$channel, spot_counts_summarized$sampletype, spot_counts_summarized$replicate))

# Plot number of counts per genotype and media condition
spot_counts_summarized$genotype <- factor(spot_counts_summarized$genotype, levels = c("wild type", "Med12 mutant"), ordered = TRUE)

# Example plot for one replicate (Fig S4H, left panel)
spot_counts_summarized_rep1 <- subset(spot_counts_summarized, replicate == "Rep1" & sampletype == "sample")
ggplot(spot_counts_summarized_rep1, aes(x = genotype, y = spot_count, fill = genotype)) + 
  geom_violin(linewidth = 0.3) +
  scale_fill_manual(values = c("#2E6A67","#E98624")) +
  geom_boxplot(width = 0.2, size = 0.4, outlier.size = 0.3, fill = "white", linewidth = 0.48) +
  facet_grid(cols = vars(media), rows = vars(channel), scales = "free", space = "free_x")  +
  xlab("genotype") +
  ylab("median mRNA counts per cell") +
  theme_bw(base_size = 10) + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black"), 
        panel.border = element_rect(color = "black", linewidth = 0.48),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        # panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA)) +
  force_panelsizes(rows = unit(3, "cm"),
                   cols = unit(3, "cm"))
ggsave("./240814_RNA_Count_distributions_Rep1.pdf", width = 6, height = 6)




# Define function for coefficient of variation
cv <- function(x) {
  cv_value <- sd(x) / mean(x)
  return(cv_value)
}

# Calculate median spot counts individually for each replicate (and condition)
spot_counts_medians <- spot_counts_summarized %>% 
  group_by(media, genotype, channel, sampletype, replicate) %>%
  summarise(median = median(spot_count), CoV = cv(spot_count))

# Plot medians of samples (Figure S4H, right panel)
spot_counts_medians <- subset(spot_counts_medians, sampletype != "negative control")
spot_counts_medians$genotype <- factor(spot_counts_medians$genotype, levels = c("wild type", "Med12 mutant"), ordered = TRUE)
ggplot(spot_counts_medians, aes(x=genotype, y=median, fill = genotype)) +
  geom_bar(position = "dodge", stat = "summary", fun = "mean", width = 0.7) +
  geom_jitter(aes(shape = factor(replicate)), position = position_jitter(0.1), color = "#4C4C4C") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.3, color = "black") +
  ylab("median mRNA counts per cell") +
  facet_grid(cols = vars(media), rows = vars(channel), scales = "free", space = "free_x") +
  scale_fill_manual(values = c("#2E6A67","#E98624")) +
  scale_y_continuous(expand = expansion(mult = c(0, .2))) +
  stat_compare_means(
    method = "t.test",
    paired = TRUE,
    label = "p.signif", hjust = 0.5, label.x = 1.5  # Adjust the position of the labels
  ) +
  theme_bw(base_size = 10) + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black"), 
        panel.border = element_rect(color = "black", linewidth = 0.48),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        # panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA)) +
  force_panelsizes(rows = unit(3, "cm"),
                   cols = unit(2, "cm"))
ggsave("./240814_Median_RNA_Count_barplot.pdf", width = 4, height = 6)

