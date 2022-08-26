library(ggplot2)
library(cowplot)
library(tidyr)
library(plyr)
library(multcompView)
library(MASS)
library(readr)
library(naniar)

## locate and read data
raw_data_sample <- read_csv("Data-CT_Plot.csv")
#View(raw_data_sample)


#Clean up the data
raw_data_sample$Group <- gsub('-','_', raw_data_sample$Group)
raw_data_sample <- subset(raw_data_sample, Group != "NC" & Group != "PC_D1" & Group != "PC_D2" & Group != "PC_D3")
raw_data_sample <- raw_data_sample %>% replace_with_na(replace = list(Ct=40.00000))

p<-ggplot(raw_data_sample, aes(x=Group, y=Ct, fill=Target)) + scale_fill_manual(values=c("#d7191c","#fdae61","#ffffbf","#abdda4","#2b83ba")) +
  geom_boxplot(show.legend = TRUE) + 
  #ggtitle("Actin", subtitle = "Not normalized, transformed") +
  xlab("Treatments") + ylab("qPCR threshold cycle") +
  scale_x_discrete(labels = str_wrap(c("Antibiotic female", "Antibiotic male", "Control female", "Control male", "Colony female", "Colony male"), width = 10)) +
  theme(axis.text = element_text(size = 11)) +
  theme(axis.title = element_text(size = 12)) 
p
ggsave(filename = "Raw-Ct.tiff", path="Paper/Graphs", width = 238, height = 185, units="mm", device='tiff', dpi=300, compression = "lzw")

 
