###QC filtering
rm(list = ls())
library(ballgown)
library(biomaRt)
library(dplyr)
library(ggplot2)

setwd("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/QC/starstandard")
##ERCC cor
namelist<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/QC/starstandard/QC_starstandard_list_withoutR.csv",header=T,row.names = 1)
ERCC <- as.data.frame(read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/QC/ERCC/TPM_pearson.csv",header=T,row.names = 1))




#X
namelist_X <- namelist[grep("X",namelist[,1]),1]
filtered_list<-which(ERCC$Column %in% namelist_X)

c<-ERCC[filtered_list,2]
c<-as.numeric(c)
c<-as.data.frame(c)
ggplot(c, aes(x = c)) +geom_density() +
  labs(#title = "Spearman correlation of ERCC",
       x = "Correlation of ERCC in Plasma", y = "Samples")+  theme_bw()+
  theme(panel.grid = element_blank())
ggsave("ERCC_density_plasma.png", width = 3, height = 3, units = "in", dpi = 300)
ggsave("ERCC_density_plasma.pdf", width = 3, height = 3, units = "in", dpi = 300)

#N
namelist_N <- namelist[grep("N",namelist[,1]),1]
filtered_list<-which(ERCC$Column %in% namelist_N)

c<-ERCC[filtered_list,2]
c<-as.numeric(c)
c<-as.data.frame(c)
ggplot(c, aes(x = c)) +geom_density() +
  labs(#title = "Spearman correlation of ERCC",
    x = "Correlation of ERCC in Plasma", y = "Samples")+  theme_bw()+
  theme(panel.grid = element_blank())
ggsave("ERCC_density_urine.png", width = 3, height = 3, units = "in", dpi = 300)
ggsave("ERCC_density_urine.pdf", width = 3, height = 3, units = "in", dpi = 300)
