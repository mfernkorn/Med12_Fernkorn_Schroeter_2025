# Plot Fig S4F
# load libraries
library(ggplot2)
library(reshape2)

# Import quantifications from csv file
dat <- read.csv("./Quantification_pER_Western_all_Replicates.csv", sep = ";")
dat <- subset(dat, Cellline != "FGF KO ")
ggplot(dat, aes(x = Tubulin, y = total.ERK44+total.ERK42, shape = Cellline)) + geom_point()

dat$Erk_ratio <- (dat$pERK44+dat$pERK42)/(dat$total.ERK44+dat$total.ERK42)
ggplot(dat, aes(x = Tubulin, y = Erk_ratio, shape = Genotype)) + geom_point()

dat_ERK_ratio <- reshape2::dcast(dat, Cellline + Replicate ~ Genotype, value.var = "Erk_ratio")

dat_ERK_ratio_long <- melt(dat_ERK_ratio, value.name = "ratio", variable.name = "Genotype")
dat_ERK_ratio_long$Condition <- paste(dat_ERK_ratio_long$Cellline, dat_ERK_ratio_long$Replicate)
dat_ERK_ratio_long$Genotype <- factor(dat_ERK_ratio_long$Genotype, levels = c("wt", "Med12 KO"))
ggplot(dat_ERK_ratio_long, aes(x=Genotype, y=ratio, colour = Cellline)) + 
  geom_line(aes(group=Condition), linetype=2, colour = "black") +
  geom_point(size = 2.5) +
  scale_color_manual(values = c('#226DAE','#7DB68D')) + 
  stat_compare_means(comparisons = list(c("wt", "Med12 KO")), 
                     method = "t.test",
                     label = "p.signif") +
  ylab("phospho-ERK / total-ERK") + 
  ylim(c(0,4)) +
  theme_bw() + theme(panel.grid = element_blank())
ggsave("./20240816_plot_ERK_ratio.pdf", width = 3, height = 3)


compare_means(ratio ~ Genotype,  data = dat_ERK_ratio_long,  method = "t.test", p.adjust.method = "BH")



dat_ERK_ratio$Genotype_Erk_ratio <- dat_ERK_ratio$`Med12 KO`/dat_ERK_ratio$wt
dat_ERK_ratio$All <- "all"
ggplot(dat_ERK_ratio, aes()) + geom_bar(aes(y = mean(dat_ERK_ratio$Genotype_Erk_ratio))) +
  theme_bw() +geom_point(aes(x=All, y=Genotype_Erk_ratio))



