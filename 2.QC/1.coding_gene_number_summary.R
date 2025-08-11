library(DESeq2)
library(apeglm)
library(regionReport)
library(ggplot2)
library(biomaRt)
library(ggrepel)
library(clusterProfiler)
setwd("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC")


count<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome/salmon.repeataware_simple_genename.csv",header = T,row.names=1)
genelist<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome/all_genelist_annotation.csv",header = T,row.names=1)

protein_coding<-genelist[which(genelist$genetype=="protein_coding"),4]

CEG<-count[which(count$Group.1 %in% protein_coding),]


filtered_namelist<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/QC_repeataware_list.csv",header=T,row.names = 1)

CEG<-CEG[,c(1,which(colnames(CEG) %in% filtered_namelist$x))]
  
#X N CEG
idx_X<-grep("X",colnames(CEG))
idx_N<-grep("N",colnames(CEG))

data_X<-CEG[,idx_X]#121
data_N<-CEG[,idx_N]#118

a<-apply(data_X,1,function(x)(length(which(x>1))))
length(which(a>117*0.8))#47
length(which(a>117*0.5))#2902
length(which(a>117*0.3))#10272

b<-apply(data_N,1,function(x)(length(which(x>1))))
length(which(b>105*0.8))#6883
length(which(b>105*0.5))#12288
length(which(b>105*0.3))#15132


data_X
c<-as.numeric(apply(data_X,2,function(x)(length(which(x>1)))))
c<-data.frame(colnames(data_X),c)
idx_S<-grep("S",c$colnames.data_X.)
c$group[idx_S]<-"S"
idx_H<-grep("H",c$colnames.data_X.)
c$group[idx_H]<-"H"
idx_C<-grep("C",c$colnames.data_X.)
c$group[idx_C]<-"C"


ggplot(c, aes(x = group, y = c, fill = group)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(), alpha = 0.7) +  # 柱状图
  geom_jitter(aes(color = group), width = 0.2,size = 1.5,alpha = 0.8) +  # 抖动散点
  theme_light(base_size = 10) +
  theme_classic()

ggplot(c, aes(x = c)) +geom_density() +
  labs(title = "", x = "Plasma Protein Coding Gene number", y = "Density")+  theme_bw()+
  theme(axis.title.y = element_blank(),  # 去掉纵坐标标题
        axis.text.y = element_blank(),   # 去掉纵坐标刻度文本
        axis.ticks.y = element_blank()) +  # 去掉纵坐标刻度
  scale_y_continuous(breaks = NULL)  # 移除纵坐标刻度
ggsave("Density Plot of Plasma Protein Coding Genes.png",height=4,width=4,dpi=1200)


data_N
d<-as.numeric(apply(data_N,2,function(x)(length(which(x>0)))))
d<-data.frame(colnames(data_N),d)
idx_S<-grep("S",d$colnames.data_N.)
d$group[idx_S]<-"S"
idx_H<-grep("H",d$colnames.data_N.)
d$group[idx_H]<-"H"
idx_C<-grep("C",d$colnames.data_N.)
d$group[idx_C]<-"C"

ggplot(d, aes(x = d)) +geom_density() +
  labs(title = "", x = "Urine Protein Coding Gene number", y = "Density")+  theme_bw()+
  theme(axis.title.y = element_blank(),  # 去掉纵坐标标题
        axis.text.y = element_blank(),   # 去掉纵坐标刻度文本
        axis.ticks.y = element_blank()) +  # 去掉纵坐标刻度
  scale_y_continuous(breaks = NULL)  # 移除纵坐标刻度
ggsave("Density Plot of Urine Protein Coding Genes.png",height=4,width=4,dpi=1200)

