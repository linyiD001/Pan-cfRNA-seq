######
setwd("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar")
rm(list=ls())
Normalized_log2_read<-read.csv("log2_expr_cpm.csv",header=T,row.names = 1)

#Normalized_log2_read <- cpm_matrix_from_star
data_all<-as.data.frame(Normalized_log2_read)

coldata<-matrix(ncol = 2,nrow = 223)
rownames(coldata)<-colnames(data_all)

idx_X_counts<-grep("X",colnames(data_all))
idx_N_counts<-grep("N",colnames(data_all))
idx_S_counts<-grep("S",colnames(data_all))
idx_H_counts<-grep("H",colnames(data_all))
idx_C_counts<-grep("C",colnames(data_all))

coldata[idx_X_counts,1]<-"Plasma"
coldata[idx_N_counts,1]<-"Urine"
coldata[idx_H_counts,2]<-"HC"
coldata[idx_S_counts,2]<-"BS"
coldata[idx_C_counts,2]<-"BC"

library(readxl)
colnames(coldata)<-c("category","group")

clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BC")
MIBC<-as.data.frame(cbind(clinicial_data$Patient,clinicial_data$MIBC,clinicial_data$`T`,clinicial_data$AGE,clinicial_data$SEX,clinicial_data$Tumor_size))
colnames(MIBC)<-c("individual","MIBC","Stage","Age","Sex","Tumor_size")

BS_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BS")
BS_cd<-as.data.frame(cbind(BS_clinicial_data$Patient,NA,"BS",BS_clinicial_data$AGE,BS_clinicial_data$SEX,NA))
colnames(BS_cd)<-c("individual","MIBC","Stage","Age","Sex","Tumor_size")

HC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "HC")
HC_cd<-as.data.frame(cbind(HC_clinicial_data$Patient,NA,"HC",HC_clinicial_data$AGE,HC_clinicial_data$SEX,NA))
colnames(HC_cd)<-c("individual","MIBC","Stage","Age","Sex","Tumor_size")

MIBC_all<-rbind(MIBC,HC_cd,BS_cd)
MIBC_all<-MIBC_all[-c(1:3),]
colnames(MIBC_all)<-c("Patients","MIBC","Stage","Age","Sex","Tumor_size")
MIBC_all<-as.data.frame(MIBC_all)
coldata<-as.data.frame(coldata)
coldata$sample_num<-rownames(coldata)
coldata$Patients<-rownames(coldata)
coldata$Patients <- sub(paste0("*", "_X"), "", coldata$Patients)
coldata$Patients <- sub(paste0("*", "_N"), "", coldata$Patients)


coldata_all<-merge(coldata,MIBC_all,by="Patients",all.x=T)
rownames(coldata_all)<-coldata_all$sample_num
coldata_all_N<-coldata_all[which(coldata_all$category=="Urine"),]
Normalized_log2_read_N<-Normalized_log2_read[,grep("N",colnames(Normalized_log2_read))]

library(reshape2)
library(dplyr)

# -----------------------------
# 设置你想画的基因
# -----------------------------
#genes_to_plot <- c("AQP2","BST2","CALR","CHCHD10","CIMAP3",
#                   "DCDC2","EFR3A","ENSG00000269590","ENSG00000289901","FABP4",
#                   "FXYD4","GPX2","H19","HMGN5","HNRNPA1",       
#                   "HSPD1P6","IGFBP3","LHX1","MIR4444-1","MTATP6P1",
#                   "PKP1","PLAAT4","RBM25","RNF152","RPL36AL",
#                   "RPS4Y1","RSRC2","SMG1P1","TFPT","TRUB1",     
#                   "TUBGCP5","ZNF818P")  # 替换为你的基因

genes_to_plot<- c("CD74","ENSG00000240801","ENSG00000284779","H19","IGF2",        
                   "IGFBP3","INS-IGF2","KRT17","MDK","MIR675",
                   "RARRES1","SLC2A1")



expr_sub <- Normalized_log2_read_N[genes_to_plot, ]
expr_sub$gene <- rownames(expr_sub)  # 把行名变成一列


expr_long <- melt(expr_sub, id.vars = "gene",
                  variable.name = "sample",
                  value.name = "expression")

expr_long <- expr_long %>%
  left_join(coldata_all_N, by = c("sample" = "sample_num"))  # sample_id 替换为你 coldata_all_N 中的样本名列

expr_long$Stage[grep('1',expr_long$Stage)]<-'T1_2'
expr_long$Stage[grep('2',expr_long$Stage)]<-'T1_2'
expr_long$Stage[grep('3',expr_long$Stage)]<-'T3_4'
expr_long$Stage[grep('4',expr_long$Stage)]<-'T3_4'
expr_long$Stage[which(expr_long$Stage=='a')]<-'Ta'
unique(expr_long$Stage)
expr_long$Stage <- factor(expr_long$Stage,c("HC","BS","Ta", "T1_2", "T3_4"), ordered = TRUE)
#expr_long$Stage <- factor(expr_long$Stage,c("HC","BS","Ta", "T1", "T2", "T3", "T4"), ordered = TRUE)
#expr_long$Stage <- factor(expr_long$Stage,c("HC","BS","early","late"), ordered = TRUE)

all_comparisons <- combn(levels(expr_long$Stage), 2, simplify = FALSE)
all_comparisons <- list(c("HC","BS"),c("BS","Ta"))
#all_comparisons <- list(c("HC","BS"),c("BS","early"))
all_comparisons <- list(c("HC","BS"),c("BS","Ta"),c("HC","Ta"))

ggplot(expr_long, aes(x = Stage, y = expression, fill = Stage)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  #geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  facet_wrap(~gene, scales = "free_y") +
  stat_compare_means(
    comparisons = all_comparisons,
    method = "wilcox.test",  # 如果数据符合正态可换成 "t.test"
    label = "p.signif",
    hide.ns = F
  ) +
  theme_bw() +
  labs(#title = "Expression of selected genes across stages",
       x = "Stage", y = "log2(normalied CPM+1)") +
  theme(
    text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c(
    "HC" = "#379BF4",
    "BS" = "#FCAE59",
    #"early" = "#FBD0DF",
    #"late" = "#DB498E"
    "Ta" = "#FBD0DF", 
    #"T1" = "#F3ADCA",
    #"T2" = "#EB8CB6",
    "T1_2" = "#F3ADCA",
    "T3_4" = "#DB498E"
    #"T3" = "#E36BA2",
    #"T4" = "#DB498E"
  ))

ggsave("gene_for_multi_N.png",height=10,width=10,dpi=300)
ggsave("gene_for_multi_N.pdf",height=10,width=10,dpi=300)
