# PLot Figure S9E
# load libraries
library(ggplot2)
library(stringr)
library(ggh4x) # for force_panelsize option
library(reshape2)
library(ggpointdensity)
library(dplyr)

# Load data from csvs, produced with Segmentation_and_Measurement.ijm
# right now necessary to type in condition names manually, take care about filename
filenames <- list.files(path = "./Csvs/")

# Loop over all files
dat <- data.frame()
for (i in 1:length(filenames)){
  temp_dat <- read.csv(paste0("./Csvs/", filenames[i]))[,2:6]
  temp_dat$Channel <- str_sub(filenames[i], start = -8, end = -5)
  if (grepl("2i", filenames[i])){
    temp_dat$Media <- "2i" 
  } else if(grepl("N2", filenames[i])){
    temp_dat$Media <- "N2B27" 
  } else {
    temp_dat$Media <- "negative control"
  }
  if (grepl("wt", filenames[i])){
    temp_dat$Genotype <- "wt" 
  } else {
    temp_dat$Genotype <- "Med12_mutant" 
  }
  temp_dat$Image <- str_extract(filenames[i], "(?<=.)\\d(?=_Merg)")
  temp_dat$Cell <- 1:length(temp_dat$Channel)
  dat <- rbind(dat, temp_dat)
}

# rename levels of Channel to more meaningful names
dat$Channel <- as.factor(dat$Channel)
levels(dat$Channel) <- c('Nuclei', 'Oct4', 'Spry4', "Nanog")

# Filter for minimal area of a cell
hist(dat$Area, breaks = 200)
dat <- subset(dat, Area >=30)

#Rearrange columns and only use IntDens
dat_reshape <- reshape(dat,
                       idvar = c("Media","Cell","Image","Genotype", "Area"),
                       timevar = "Channel",
                       direction = "wide")

# obtian median values for each condition for text
table <- dat_reshape %>%
  group_by(Media, Genotype) %>%
  summarize(number_of_cells = max(Cell),
            median_Area = median(Area),
            median_mean_Spry4 = median(Mean.Spry4), median_IntDen_Spry4 = median(IntDen.Spry4), 
            median_mean_Nanog = median(Mean.Nanog), median_IntDen_Nanog = median(IntDen.Nanog),
            median_mean_Oct4 = median(Mean.Oct4),median_IntDen_Oct4 = median(IntDen.Oct4))
# write.csv2(table, file ="./Plots/20240603_Summary.csv")

# Plot Oct4 intensities (Figure S9E)
dat_reshape$Genotype <- factor(dat_reshape$Genotype, levels = c("wt", "Med12_mutant"), ordered = TRUE)
ggplot(subset(dat_reshape, Media != "negative control"), aes(x = Genotype, y=IntDen.Oct4, fill = Genotype)) + 
  geom_violin(linewidth = 0.3) +
  scale_fill_manual(values = c("#2E6A67","#E98624")) +
  geom_boxplot(width = 0.3, outlier.size = 0.2, fill = "white", linewidth = 0.2) +  
  scale_y_log10() +
  ylab("Oct4 fluorescence [au]") +
  xlab("genotype") +
  facet_wrap(facets = vars(Media)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black"), 
        panel.border = element_rect(color = "black", linewidth = 0.48),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA)) +
  force_panelsizes(rows = unit(4, "cm"),
                   cols = unit(2, "cm"))
ggsave("./Plots/20240814_IntDen_Oct4_Figure.pdf", width = 4, height = 4)
