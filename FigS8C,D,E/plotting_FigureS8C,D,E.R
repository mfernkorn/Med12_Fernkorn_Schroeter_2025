# Plot Figure S8 C (lower panel), D and E
# Load packages
library(ggplot2)
library(stringr)
library(ggh4x) # for force_panelsize option
library(reshape2)
library(ggpointdensity)
library(viridis)
library(mclust)
library(tidyverse)

# Load data from all replicates
# Define read in function
read_in_Intensities <- function(rel_filepath, Rep_ID){
  filenames <- list.files(path = rel_filepath)
  
  dat <- data.frame()
  for (i in 1:length(filenames)){
    temp_dat <- read.csv(paste0(rel_filepath, filenames[i]))[,2:6]
    temp_dat$Channel <- str_sub(filenames[i], start = -8, end = -5)
    temp_dat$Condition <- str_extract(filenames[i], ".*?(?=_Cer)")
    if (grepl("1h", filenames[i])){
      temp_dat$Labeling_Time <- "1h" 
    } else if (grepl("30min", filenames[i])){
      temp_dat$Labeling_Time <- "30min" 
    } else {
      temp_dat$Labeling_Time <- "1h" 
    }
    if (grepl("2i", filenames[i])){
      temp_dat$Media <- "2i" 
    } else {
      temp_dat$Media <- "N2B27" 
    }
    if (grepl("wtCer", filenames[i])){
      temp_dat$Labeled <- "wt-Cerulean" 
    } else {
      temp_dat$Labeled <- "Med12_mutant-Cerulean" 
    }
    temp_dat$Image <- str_extract(filenames[i], "(?<=.)\\d(?=_Merg)")
    temp_dat$Cell <- 1:length(temp_dat$Condition)
    dat <- rbind(dat, temp_dat)
  }
  
  # rename levels of Channel to more meaningful names
  dat$Channel <- as.factor(dat$Channel)
  levels(dat$Channel) <- c('Cerulean', 'Spry4', 'RNA', "Nuclei")
  dat$Rep <- Rep_ID
  return(dat)
}

# Load data from all replicates 
# (generated from immunostainig images with the ImageJ macro "Segmentation_and_Measurement_v2.0_adjusted_for_4channels")
# Replicate 1
rep1 <- read_in_Intensities("./20240425_Rep1/Analysis/Csvs/", "Rep1")
# Replicate 2
rep2 <- read_in_Intensities("./20240503_Rep2/Analysis/Csvs/", "Rep2")
# Replicate 3
rep3 <- read_in_Intensities("./20240531_Rep3/Analysis/Csvs/", "Rep3")
# Replicate 4
rep4 <- read_in_Intensities("./20240604_Rep4/Analysis/Csvs/", "Rep4")

# Merge Data of replicates
dat <- bind_rows(rep1, rep2, rep3, rep4)

# Filter for minimal area of a cell
hist(dat$Area, breaks = 200)
dat <- subset(dat, Area >=40)

#Rearrange columns (reformat to wide format for each channel)
dat_reshape <- pivot_wider(dat[,c(7:13, 1, 6, 2, 4)], 
                           names_from = Channel,
                           values_from = c(Mean, IntDen))

# Apply GMM for each condition (per Condition, Rep, Media, Labeled, Labeling_Time)
dat_mclust <- dat_reshape %>% 
  group_by(Condition, Rep, Media, Labeled, Labeling_Time) %>%
  summarize(GMM = list(Mclust(log10(Mean_Cerulean), G = 2))) # Fixing G=2 for mixture of two gaussians
# Assign most important GMM parameters to dataframe columns manually
dat_mclust$Sd_1 <- sapply(dat_mclust$GMM, \(x) sqrt(x$parameters$variance$sigmasq))[1,]
dat_mclust$Sd_2 <- sapply(dat_mclust$GMM, \(x) sqrt(x$parameters$variance$sigmasq))[2,]
dat_mclust$Mean_1 <- sapply(dat_mclust$GMM, \(x) x$parameters$mean)[1,]
dat_mclust$Mean_2 <- sapply(dat_mclust$GMM, \(x) x$parameters$mean)[2,]
dat_mclust$Weight_1 <- sapply(dat_mclust$GMM, \(x) x$parameters$pro)[1,]
dat_mclust$Weight_2 <- sapply(dat_mclust$GMM, \(x) x$parameters$pro)[2,]

# Define function as function for probability
posterior_prob <- function(x, mean1, sd1, mean2, sd2, weight1, weight2) {
  density1 <- weight1 * dnorm(x, mean1, sd1)
  density2 <- weight2 * dnorm(x, mean2, sd2)
  total_density <- density1 + density2
  # Handle cases where total_density is zero to avoid division by zero
  if (total_density == 0) {
    return(NA)
  }
  post_prob1 <- density1 / total_density
  return(post_prob1)
}

threshold_low <- c()
threshold_high <- c()
for (i in 1:length(dat_mclust$Condition)){
  # Find the thresholds where posterior probability is 0.05 and 0.95
  # Define data range based on maximum of probability function (use range between maximum and max of data)
  lowerbound <- optimize(function(x) posterior_prob(x, dat_mclust$Mean_1[i], dat_mclust$Sd_1[i], dat_mclust$Mean_2[i], dat_mclust$Sd_2[i], dat_mclust$Weight_1[i], dat_mclust$Weight_2[i]), interval = range(dat_mclust$GMM[[i]]$data), maximum = TRUE)$maximum
  higherbound <- max(dat_mclust$GMM[[i]]$data)
  # Calculate thresholds
  threshold_low[i] <- uniroot(function(x) posterior_prob(x, dat_mclust$Mean_1[i], dat_mclust$Sd_1[i], dat_mclust$Mean_2[i], dat_mclust$Sd_2[i], dat_mclust$Weight_1[i], dat_mclust$Weight_2[i]) - 0.15,
                         interval = range(lowerbound, higherbound), extendInt = "yes")$root

  threshold_high[i] <- uniroot(function(x) posterior_prob(x, dat_mclust$Mean_1[i], dat_mclust$Sd_1[i], dat_mclust$Mean_2[i], dat_mclust$Sd_2[i], dat_mclust$Weight_1[i], dat_mclust$Weight_2[i]) - 0.85,
                          interval = range(lowerbound, higherbound), extendInt = "yes")$root

  print(ggplot(subset(dat_reshape, Condition == dat_mclust$Condition[i] & Rep == dat_mclust$Rep[i] & Media == dat_mclust$Media[i] & Labeled == dat_mclust$Labeled[i] & Labeling_Time == dat_mclust$Labeling_Time[i]), 
         aes(x=log10(Mean_Cerulean))) +
    geom_histogram(aes(y = after_stat(count / sum(count))),bins = 300) +
    geom_function(fun = function(x) dnorm(x, mean = dat_mclust$Mean_1[i], sd = dat_mclust$Sd_1[i])*(dat_mclust$Weight_1[i])*0.005, colour = "red") +
    geom_function(fun = function(x) dnorm(x, mean = dat_mclust$Mean_2[i], sd = dat_mclust$Sd_2[i])*(dat_mclust$Weight_2[i])*0.005, colour = "blue") +
    geom_vline(xintercept = threshold_low[i], colour = "blue", linetype = "dotted") +
    geom_vline(xintercept = threshold_high[i], colour = "red",linetype = "dotted") +
    ggtitle(paste(dat_mclust$Rep[i], dat_mclust$Media[i], dat_mclust$Labeled[i], dat_mclust$Labeling_Time[i])))
  # ggsave(paste0("./Compare_Replicates_Plots/20240605_histogram_gmm_with_thresholds_", paste(dat_mclust$Rep[i], dat_mclust$Media[i], dat_mclust$Labeled[i], dat_mclust$Labeling_Time[i], sep = "_"), ".pdf"), width = 7, height = 7)
  
}


# Assign Genotypes based on Cerulean expression in combination with Condition
dat_reshape <- data.frame(dat_reshape)
rownames(dat_reshape) <- 1:length(rownames(dat_reshape)) #?
dat_reshape$Assigned_Genotype <- "Nan"
for (j in 1:length(dat_mclust$Condition)){
  current_rows <- as.numeric(rownames(subset(dat_reshape, Condition == dat_mclust$Condition[j] & Rep == dat_mclust$Rep[j] & Media == dat_mclust$Media[j] & Labeled == dat_mclust$Labeled[j] & Labeling_Time == dat_mclust$Labeling_Time[j])))
  for (i in current_rows){
    if (dat_reshape$Mean_Cerulean[i] > 10^threshold_low[j] && dat_reshape$Labeled[i] == "Med12_mutant-Cerulean"){
      dat_reshape$Assigned_Genotype[i] <- "-> Med12_mutant"
    } else if (dat_reshape$Mean_Cerulean[i] > 10^threshold_low[j] && dat_reshape$Labeled[i] == "wt-Cerulean"){
      dat_reshape$Assigned_Genotype[i] <- "-> wt"
    } else if (dat_reshape$Mean_Cerulean[i] < 10^threshold_high[j] && dat_reshape$Labeled[i] == "Med12_mutant-Cerulean"){
      dat_reshape$Assigned_Genotype[i] <- "-> wt"
    } else if (dat_reshape$Mean_Cerulean[i] < 10^threshold_high[j] && dat_reshape$Labeled[i] == "wt-Cerulean"){
      dat_reshape$Assigned_Genotype[i] <- "-> Med12_mutant"
    } else {
      dat_reshape$Assigned_Genotype[i] <- "-> unassigned"
    }
  }
}

dat_reshape$Geno_Med <- factor(paste(dat_reshape$Assigned_Genotype, dat_reshape$Media), levels = c("-> wt 2i", "-> Med12_mutant 2i", "-> wt N2B27", "-> Med12_mutant N2B27","-> unassigned 2i", "-> unassigned N2B27"), ordered = TRUE)


# Plot variable for one replicate (Figure S8 C lower panel)
ggplot(subset(dat_reshape, Rep == "Rep1" & !(Geno_Med %in% c("-> unassigned 2i", "-> unassigned N2B27"))), aes(y=IntDen_RNA, x=Geno_Med)) +
  geom_violin(aes(fill=Geno_Med), position = position_dodge(width = 0.9), scale = "width",width = 0.85, linewidth = 0.2) +
  geom_boxplot(position = position_dodge(width = 0.9),width = 0.2,outlier.size = 0.3,linewidth = 0.2) + scale_y_log10() +
  scale_fill_manual(values = c("#2E6A67","#E98624","#91E2DA","#F7D1A5"))+
  facet_grid(cols = vars(Labeled, Condition), rows = vars(Labeling_Time), scales = "free_x") +
  ylab("RNA") + 
  xlab("Assigned Genotype") +
  theme_bw(base_size = 10) + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(color = 'black'), 
        axis.text.x = element_text(colour = "black", angle = 90,hjust = 0.9),
        axis.text.y = element_text(colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.48),
        panel.border = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(fill = NA, linewidth = 0),
        # legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA),
        strip.background = element_rect(fill = "white", color = "white"),  # Remove the border
        strip.text = element_text(color = "black")) +
  force_panelsizes(rows = unit(4, "cm"),
                   cols = unit(2.8, "cm"))
ggsave("./Compare_Replicates_Plots/20240815_IntDenRNA_Rep1_including_negative_control.pdf",  width = 12, height = 4)

# Calculate medum of variables
dat_reshape_summerized <- dat_reshape %>% 
  group_by(Condition, Rep, Media, Labeled, Labeling_Time, Assigned_Genotype, Geno_Med) %>%
  summarize(Median_Area = median(Area), Median_IntDenRNA = median(IntDen_RNA), Median_IntDenSpry = median(IntDen_Spry4), Median_IntDenNuc = median(IntDen_Nuclei),
            Median_MeanRNA = median(Mean_RNA), Median_MeanSpry = median(Mean_Spry4), Median_MeanNuc = median(Mean_Nuclei))


# Filter for negative control and unassigned genotypes
dat_reshape_summerized_subset <-  subset(dat_reshape_summerized,Condition != "negative_control" & !(Geno_Med %in% c("-> unassigned 2i", "-> unassigned N2B27")))
       
# Plot medium of one variable (Figure S8 D)
# Scale values by dividing by -> wt 2i average
dat_reshape_norm <- dat_reshape_summerized_subset %>% 
  group_by(Labeled, Labeling_Time, Rep) %>%
  mutate(Median_Area_Norm = Median_Area/Median_Area[Geno_Med == "-> wt 2i"],
         Median_IntDenSpry_Norm = Median_IntDenSpry/Median_IntDenSpry[Geno_Med == "-> wt 2i"],
         Median_IntDenRNA_Norm = Median_IntDenRNA/Median_IntDenRNA[Geno_Med == "-> wt 2i"])

ggplot(dat_reshape_norm, aes(x=Geno_Med, y=Median_IntDenRNA_Norm)) +
  geom_bar(aes(fill = Geno_Med), position = "dodge", stat = "summary", fun = "mean", width = 0.7) +
  scale_fill_manual(values = c("#2E6A67","#E98624","#91E2DA","#F7D1A5")) +
  geom_jitter(aes(shape = Rep, color = Labeling_Time), position = position_jitter(0.2)) +
  scale_color_manual(values = c("black", "grey"))+
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = 0.3, color = "#144976") +
  scale_y_continuous(expand = c(0,0,0.1,0.1)) + 
  ylab("median RNA fluorescence\nnormalized to wild type in 2i")+
  xlab("Genotype and Medium")+
  facet_grid(cols = vars(Labeled), scales = "free_y") +
  scale_x_discrete(guide = guide_axis(angle = 90))+
  stat_compare_means(comparisons = list(c("-> Med12_mutant 2i", "-> Med12_mutant N2B27"), c("-> wt 2i", "-> wt N2B27")), 
                     method = "t.test",
                     paired = TRUE, 
                     p.adjust.method = "BH",
                     label = "p.signif") +
  theme_bw(base_size = 10) + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(color = 'black'), 
        axis.text = element_text(colour = "black"), 
        panel.border = element_rect(color = "black", linewidth = 0.48),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        plot.background = element_rect(fill = NA, linewidth = 0),
        # legend.position = "none",
        legend.background = element_rect(fill = NA, linewidth = 0),
        legend.key = element_rect(fill = NA)) +
  force_panelsizes(rows = unit(4, "cm"),
                   cols = unit(4, "cm"))
ggsave("./Compare_Replicates_Plots/20240815_Median_IntDenRNA_norm_stats_combined_times.pdf",  width = 6, height = 5)

compare_means(Median_IntDenRNA_Norm ~ Geno_Med, 
              group.by = "Labeled", 
              data = dat_reshape_norm, 
              method = "t.test",
              p.adjust.method = "BH",
              paired = TRUE)
