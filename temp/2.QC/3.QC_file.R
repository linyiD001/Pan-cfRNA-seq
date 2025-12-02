#####plot of unfiltered data (STAR mapping)

setwd("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC")
rm(list = ls())
library(readxl)
library(ggplot2)
library(ggsignif)
library(readxl)
library(ggplot2)
library(ggsignif)
library(readxl)

QC_STAR_1 <- as.data.frame(read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/240905cfRNAQC.xlsx",sheet = "Sheet1"))
QC_STAR_2 <- as.data.frame(read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/240905cfRNAQC.xlsx",sheet = "star_log"))
QC_STAR_3 <- as.data.frame(read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/240905cfRNAQC.xlsx",sheet = "samtools"))
QC_STAR_4 <- as.data.frame(read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/240905cfRNAQC.xlsx",sheet = "read_distribution"))
QC_STAR_5 <- as.data.frame(read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/ERCC/TPM_pearson.csv",header=T,row.names = 1))

QC_salmon_standard<-as.data.frame(read.table("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/multiqc_general_stats_standard.txt",header=T))
QC_salmon_repeataware<-as.data.frame(read.table("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/multiqc_general_stats_repeataware.txt",header=T))

###QC_parameter test
##mapping reads 100000(salmon)
##Spearman correlation of spike-in > 0.5
##STAR exon/intron <10

samplelist_1<-as.vector(QC_STAR_1[which(QC_STAR_1$clean>10000000),1])
samplelist_2<-as.vector(QC_STAR_2[which(QC_STAR_2$`total map`>100000),1])
samplelist_3<-as.vector(QC_STAR_5[which(QC_STAR_5$P>0.5),1])
samplelist_4<-as.vector(QC_STAR_4[which(QC_STAR_4$`Intron/exon`<10),1])
samplelist_5<-as.vector(QC_STAR_4[which(QC_STAR_4$` 3|bias `<1),1]) ##all passed
samplelist_6<-as.vector(QC_STAR_3[which(QC_STAR_3$rRNA/QC_STAR_3$hg38<0.4),1])
samplelist_7<-as.vector(QC_salmon_repeataware[which(QC_salmon_repeataware$salmon.num_mapped>0.1),1])

#QC_repeataware_STAR<-Reduce(intersect, list(samplelist_1,samplelist_2,samplelist_3,samplelist_4,samplelist_5,samplelist_6))
#write.csv(QC_repeataware_STAR,"QC_repeataware_STAR_list.csv")

QC_repeataware<-Reduce(intersect, list(samplelist_1,samplelist_2,samplelist_3,samplelist_4,samplelist_5,samplelist_6,samplelist_7))
QC_repeataware<-as.data.frame(QC_repeataware)
QC_repeataware<-as.data.frame(QC_repeataware[-grep("C107",QC_repeataware[,1]),])
QC_repeataware<-as.data.frame(QC_repeataware[-grep("C89",QC_repeataware[,1]),])
write.csv(QC_repeataware,"QC_repeataware_list_withoutR.csv")

