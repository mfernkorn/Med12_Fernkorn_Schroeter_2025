# Plot Fig S4C
# load libraries
library(ggplot2)
library("scales")
library(reshape2)

# Load FACS data, exported from FlowJo
medians_rep1 <- read.csv2("./29-Aug-2022.wsp FlowJo_medians_Rep1.csv", sep = ",", row.names = NULL)
medians_rep2 <- read.csv2("./19-Aug-2022.wsp FlowJo_medians_Rep2.csv", sep = ",", row.names = NULL)
medians_rep3and4 <- read.csv2("./12-Sep-2022.wsp FlowJo_medians.csv", sep = ",", row.names = NULL)


medians <- rbind(medians_rep1, medians_rep2, medians_rep3and4)
  
ggplot(medians, aes(x = Time, y=Median, color = Genotype, shape = factor(Replicate))) + 
  geom_point() +
  scale_shape_manual(values=c(15,16,17, 18)) + 
  scale_color_manual(values = c("black", "grey")) +
  scale_size_manual(values=c(2,2,2)) +
  theme_bw() +theme(aspect.ratio=1)
ggsave("./20220919_Spry4exp_medians.pdf", width = 4, height = 4)
