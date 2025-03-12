# Plot FigS1
#load libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library("varhandle")
library("drc")

# Load list of differentially expressed genes, produced in Seqmonk as described in the methods section
RPMs_RNA_seq = read.csv("./hitlist2.txt", header = TRUE, sep = "\t")
RPMs_RNA_seq = data.frame(RPMs_RNA_seq$Probe, RPMs_RNA_seq$X0, RPMs_RNA_seq$X0.625, RPMs_RNA_seq$X1.25, RPMs_RNA_seq$X2.5, RPMs_RNA_seq$X5, RPMs_RNA_seq$X10, RPMs_RNA_seq$X20, RPMs_RNA_seq$X40)
names(RPMs_RNA_seq) = c("Gene", 0, 0.625, 1.25, 2.5, 5, 10, 20, 40)

#Filter for minimal Expression
RPMs_expressedGenes = melt(RPMs_RNA_seq)
RPMs_expressedGenes$value = 2**RPMs_expressedGenes$value

# Calculate Foldchange
CalcFoldchange_withfilter <- function(probelist, threshold){
  #probelist$Gene = droplevels(probelist$Gene)
  Genes = levels(as.factor(probelist$Gene))
  output = data.frame()
  temp_gene = subset(probelist, Gene == Genes[1])
  temp_gene[,3] = temp_gene[,3]/temp_gene[1,3]
  if (FALSE %in% c(temp_gene$value < threshold)){ # > for upregulated, < for downregulated
    keep_not = TRUE
  } else{
    keep_not = FALSE
  }
  output = temp_gene
  for (i in 2:length(Genes)){
    temp_gene = subset(probelist, Gene == Genes[i])
    temp_gene[,3] = temp_gene[,3]/temp_gene[1,3]
    if (TRUE %in% c(temp_gene$value < threshold)){ # > for upregulated, < for downregulated
      output[(length(output[,1])+1):(length(output[,1])+length(temp_gene[,1])),1:3] = temp_gene
    }
  }
  if (keep_not == TRUE){
    output = output[length(temp_gene[,1])+1:length(output[,1]),]
  }
  return(output)
}

RPMs_upredgulated_genes = CalcFoldchange_withfilter(probelist = RPMs_expressedGenes, threshold = 35)

ScaleFoldchange <- function(probelist){
  #probelist$Gene = droplevels(probelist$Gene)
  Genes = levels(as.factor(probelist$Gene))
  output = data.frame()
  temp_gene = subset(probelist, Gene == Genes[1])
  temp_gene[,3] = (temp_gene[,3]-min(temp_gene[,3]))/(max(temp_gene[,3])-min(temp_gene[,3]))
  output = temp_gene
  for (i in 2:length(Genes)){
    temp_gene = subset(probelist, Gene == Genes[i])
    temp_gene[,3] = (temp_gene[,3]-min(temp_gene[,3]))/(max(temp_gene[,3])-min(temp_gene[,3]))
    output[(length(output[,1])+1):(length(output[,1])+length(temp_gene[,1])),1:3] = temp_gene
  }
  return(output)
}

RPMs_upredgulated_genes_scaled = ScaleFoldchange(RPMs_upredgulated_genes)


filterbyfitting <- function(probelist){
  output = data.frame()
  Genes = levels(as.factor(probelist$Gene))
  for (i in 1:length(Genes)){
    temp_gene = subset(probelist, Gene==Genes[i])
    temp_gene$variable = unfactor(temp_gene$variable)
    fit1 <- try(drm(value ~ variable, data = temp_gene, fct = LL.4()), silent=F)
    if("try-error" %in% class(fit1)){
      next
    }
    #    plot(fit1, log = "")
    q = sd(predict(fit1, temp_gene)-temp_gene$value)
    if (q < 0.1){
      if (length(output)==0){
        output = temp_gene
      } else{
        output[(length(output[,1])+1):(length(output[,1])+length(temp_gene[,1])),1:3] = temp_gene
      }
    }
  }
  return(output)
}

RPMs_upregulated_genes_scaled_filtered = filterbyfitting(RPMs_upredgulated_genes_scaled)

#plotting of most expressed genes
RPMs_upredgulated_genes = subset(RPMs_upredgulated_genes, Gene %in% RPMs_upregulated_genes_scaled_filtered$Gene)
RPMs_upredgulated_genes <- with(RPMs_upredgulated_genes,  RPMs_upredgulated_genes[order(value, decreasing = T ) , ])

RPMs_upredgulated_genes$variable = factor(RPMs_upredgulated_genes$variable)
ggplot(data = RPMs_upredgulated_genes, aes(x=variable, y=value, group=Gene, color=Gene)) + 
  geom_line(size=1) + #geom_point(size=3.5)  +
  scale_y_log10() + annotation_logticks(sides="l") +
  labs(x = "Concentration of  FGF4 [ng/mL]", y = "Expression Foldchange", size=200) + theme_bw() +
  theme(axis.text=element_text(size=16), axis.title =element_text(size = 18))
ggsave("20210407_10most_upregulated_genes.png", width = 6.5, height = 5)


