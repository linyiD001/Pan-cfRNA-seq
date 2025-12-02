#####plot of unfiltered data (STAR mapping)

setwd("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC")
rm(list = ls())
library(readxl)
library(ggplot2)
library(ggsignif)
library(readxl)
library(ggplot2)
library(ggsignif)


QC_STAR_1 <- as.data.frame(read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/240905cfRNAQC.xlsx",sheet = "Sheet1"))
colnames(QC_STAR_1)[1]<-"sample_name"

QC_STAR_2 <- as.data.frame(read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/240905cfRNAQC.xlsx",sheet = "star_log"))
colnames(QC_STAR_2)[1]<-"sample_name"
QC_STAR_3 <- as.data.frame(read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/240905cfRNAQC.xlsx",sheet = "samtools"))
colnames(QC_STAR_3)[1]<-"sample_name"
QC_STAR_4 <- as.data.frame(read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/240905cfRNAQC.xlsx",sheet = "read_distribution"))
colnames(QC_STAR_4)[1]<-"sample_name"
QC_STAR_5 <- as.data.frame(read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/ERCC/TPM_pearson.csv",header=T,row.names = 1))
colnames(QC_STAR_5)[1]<-"sample_name"

data_merge<-merge(QC_STAR_1,QC_STAR_2,by="sample_name",all=T)
data_merge<-merge(data_merge,QC_STAR_3,by="sample_name",all=T)
data_merge<-merge(data_merge,QC_STAR_4,by="sample_name",all=T)
data_merge<-merge(data_merge,QC_STAR_5,by="sample_name",all=T)

#QC_salmon_standard<-as.data.frame(read.table("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/multiqc_general_stats_standard.txt",header=T))
QC_salmon_repeataware<-as.data.frame(read.table("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/multiqc_general_stats_repeataware.txt",header=T))



##filtered
filtered_namelist<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/QC_repeataware_list.csv",header=T,row.names = 1)

idx_filtered_namelist<-which(data_merge$sample_name %in% filtered_namelist$x)
#QCofseq_summary<-data_merge[idx_filtered_namelist,]
QCofseq_summary<-data_merge[-277,]


idx_X<-grep("X",QCofseq_summary$sample_name)
QCofseq_summary$TYPE[idx_X]<-"Plasma"
idx_N<-grep("N",QCofseq_summary$sample_name)
QCofseq_summary$TYPE[idx_N]<-"Urine"


QCofseq_summary$raw<-as.numeric(QCofseq_summary$raw)
p1 <- ggplot(QCofseq_summary, aes(x = TYPE, y = raw, fill = TYPE)) +
  geom_jitter(aes(color = raw < 1e+7), width = 0.2, alpha = 0.6, size = 2) +  # 根据 raw 值是否小于 1e+7 控制颜色
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +  # 添加盒形图，不显示离群值
  theme_light(base_size = 10) +
  theme_classic() +
  labs(x = "Sample Type", y = "Raw Count",color="QC_pass") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     labels = c("TRUE" = "fail", "FALSE" = "pass"))  # 修改为 pass 和 fail
p1
ggsave("raw_count.png", plot = p1, width = 4, height = 4, dpi = 300,units='in')
ggsave("raw_count.pdf", plot = p1, width = 4, height = 4, dpi = 300)

QCofseq_summary$clean<-as.numeric(QCofseq_summary$clean)
p2<-ggplot(QCofseq_summary, aes(x = TYPE, y = clean, fill = TYPE)) +
  geom_jitter(aes(color = clean < 1e+7), width = 0.2, alpha = 0.6, size = 2) +  # 根据 raw 值是否小于 1e+7 控制颜色
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +  # 添加盒形图，不显示离群值
  theme_light(base_size = 10) +
  theme_classic() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     labels = c("TRUE" = "fail", "FALSE" = "pass"))  +# 修改为 pass 和 fail
  theme_classic()+labs(x = "sample type", y = "clean count",color="QC_pass")
p2
ggsave("clean_count.png", plot = p2, width = 4, height = 4, dpi = 300,units='in')
ggsave("clean_count.pdf", plot = p2, width = 4, height = 4, dpi = 300)

QCofseq_summary$`total map`<-as.numeric(QCofseq_summary$`total map`)
p3<-ggplot(QCofseq_summary, aes(x = TYPE, y = `total map`, fill = TYPE)) +
  geom_jitter(aes(color = `total map` < 1e+5), width = 0.2, alpha = 0.6, size = 2) +  # 根据 raw 值是否小于 1e+7 控制颜色
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +  # 添加盒形图，不显示离群值
  theme_light(base_size = 10) +
  theme_classic() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     labels = c("TRUE" = "fail", "FALSE" = "pass"))  +# 修改为 pass 和 fail
  theme_classic()+labs(x = "sample type", y = "total map reads",color="QC_pass")
p3
ggsave("total_map.png", plot = p3, width = 4, height = 4,  dpi = 300,units='in')
ggsave("total_map.pdf", plot = p3, width = 4, height = 4, dpi = 300)

QCofseq_summary$`gene number(ensg)`<-as.numeric(QCofseq_summary$`基因数（ENSG类型）`)
p4<-ggplot(QCofseq_summary, aes(x = TYPE, y = `gene number(ensg)`, fill = TYPE)) +
  geom_boxplot() +
  theme_light(base_size=10)+
  theme_classic()
ggsave("gene_number(ensg).png", plot = p4, width = 3, height = 3, units = "in", dpi = 300)
ggsave("gene_number(ensg).pdf", plot = p4, width = 3, height = 3, dpi = 300)


QCofseq_summary$`Intron/exon`<-as.numeric(QCofseq_summary$`Intron/exon`)

p5<-ggplot(QCofseq_summary, aes(x = TYPE, y = `Intron/exon`, fill = TYPE)) +
  geom_jitter(aes(color = `Intron/exon` > 10), width = 0.2, alpha = 0.6, size = 2) +  # 根据 raw 值是否小于 1e+7 控制颜色
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +  # 添加盒形图，不显示离群值
  theme_light(base_size = 10) +
  theme_classic() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     labels = c("TRUE" = "fail", "FALSE" = "pass"))  +# 修改为 pass 和 fail
  theme_classic()+labs(x = "sample type", y = "Intron/exon",color="QC_pass")+ylim(0,20)
ggsave("Intron_exon.png", plot = p5, width = 4, height = 4, units = "in", dpi = 300)
ggsave("Intron_exon.pdf", plot = p5, width = 4, height = 4,  dpi = 300)

QCofseq_summary$`P`<-as.numeric(QCofseq_summary$P)
p6<-ggplot(QCofseq_summary, aes(x = TYPE, y =  P, fill = TYPE)) +
  geom_jitter(aes(color = P < 0.5), width = 0.2, alpha = 0.6, size = 2) +  # 根据 raw 值是否小于 1e+7 控制颜色
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  theme_light(base_size=10)+
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                                               labels = c("TRUE" = "fail", "FALSE" = "pass"))  +# 修改为 pass 和 fail
  theme_classic()+labs(x = "sample type", y = "correlation of ERCC(pearson)",color="QC_pass")
p6
ggsave("ERCC_cor.png", plot = p6, width = 3, height = 3, units = "in", dpi = 300)
ggsave("ERCC_cor.pdf", plot = p6, width = 3, height = 3, units = "in", dpi = 300)


QCofseq_summary$`3' bias`<-as.numeric(QCofseq_summary$`3' bias`)
p7<-ggplot(QCofseq_summary, aes(x = TYPE, y =  `3' bias`, fill = TYPE)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  geom_jitter(aes(color = `3' bias`>0.3), width = 0.2, alpha = 0.6, size = 2) +  
  theme_light(base_size=10)+ scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),labels = c("TRUE" = "fail", "FALSE" = "pass"))  +# 修改为 pass 和 fail
  theme_classic()+labs(x = "sample type", y = "3' bias",color="QC_pass")+ylim(0,0.3)
p7
ggsave("3'bias.png", plot = p7, width = 4, height = 4, units = "in", dpi = 300)
ggsave("3'bias.pdf", plot = p7, width = 4, height = 4, units = "in", dpi = 300)


QCofseq_summary$`total map%`<-as.numeric(QCofseq_summary$`total map%`)
p8<-ggplot(QCofseq_summary, aes(x = TYPE, y =  `total map%`, fill = TYPE)) +
  geom_boxplot() +
  theme_light(base_size=10)+
  theme_classic()
ggsave("total_map_percentage.png", plot = p8, width = 3, height = 3, units = "in", dpi = 300)
ggsave("total_map_percentage.pdf", plot = p8, width = 3, height = 3, units = "in", dpi = 300)


p9<-ggplot(QCofseq_summary, aes(x = P)) +geom_density() +
  labs(title = "Spearman correlation of ERCC", x = "Spearman correlation of ERCC", y = "Samples")+  theme_bw()+
  theme(panel.grid = element_blank())
ggsave("density_ERCC_cor.png", plot = p9, width = 3, height = 3, units = "in", dpi = 300)
ggsave("density_ERCC_cor.pdf", plot = p9, width = 3, height = 3, units = "in", dpi = 300)


#QC_salmon_standard<-as.data.frame(read.table("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/multiqc_general_stats_standard.txt",header=T))
QC_salmon_repeataware<-as.data.frame(read.table("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/multiqc_general_stats_repeataware.txt",header=T))

filtered_namelist<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/QC_repeataware_list.csv",header=T,row.names = 1)
#idx_filtered_namelist<-which(QC_salmon_standard$Sample %in% filtered_namelist$x)
#QC_salmon_standard_filtered<-QC_salmon_standard[idx_filtered_namelist,]
idx_filtered_namelist<-which(QC_salmon_repeataware$Sample %in% filtered_namelist$x)
QC_salmon_repeataware_filtered<-QC_salmon_repeataware[idx_filtered_namelist,]


QCofseq_summary<-QC_salmon_repeataware_filtered
idx_X<-grep("X",QCofseq_summary$Sample)
QCofseq_summary$TYPE[idx_X]<-"Plasma"
idx_N<-grep("N",QCofseq_summary$Sample)
QCofseq_summary$TYPE[idx_N]<-"Urine"

p10<-ggplot(QCofseq_summary, aes(x = TYPE, y =salmon.percent_mapped, fill = TYPE)) +
  geom_boxplot() +
  theme_light(base_size=10)+
  theme_classic()+labs(x = "sample type", y = "salmonstandard_mapping_rate")
p10
ggsave("salmonrepeataware_mapping_rate.png", plot = p10, width = 3, height = 3, units = "in", dpi = 300)
ggsave("salmonrepeataware_mapping_rate.pdf", plot = p10, width = 3, height = 3, units = "in", dpi = 300)


p11<-ggplot(QCofseq_summary, aes(x = TYPE, y =salmon.num_mapped, fill = TYPE)) +
  geom_boxplot() +
  theme_light(base_size=10)+
  theme_classic()+labs(x = "sample type", y = "salmonstandard_total map reads")
p11
ggsave("salmonrepeataware_total_map.png", plot = p11, width = 3, height = 3, units = "in", dpi = 300)
ggsave("salmonrepeataware_total_map.pdf", plot = p11, width = 3, height = 3, units = "in", dpi = 300)



p12<-ggplot(QCofseq_summary, aes(x = `gene number（ENSG）`,color = TYPE)) +geom_density() +
  labs(title = "", x = "Protein Coding Gene Number", y = "density")+  theme_bw()+
  theme(panel.grid = element_blank())
p12
ggsave("density_coding_gene.png", plot = p12, width = 3, height = 3, units = "in", dpi = 300)
ggsave("density_coding_gene.pdf", plot = p12, width = 3, height = 3, units = "in", dpi = 300)


