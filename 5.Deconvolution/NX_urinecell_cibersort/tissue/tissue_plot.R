library(pheatmap)
setwd("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/cell_type_marker/NX_urinecell_cibersort/tissue")
data<-read.table("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/cell_type_marker/NX_urinecell_cibersort/tissue/CIBERSORTx_Adjusted.txt",header=T,row.names = 1)
data<-as.data.frame(data)
data<-data[-which(rownames(data)=='C103P'),]
data<-data[-which(rownames(data)=='C89P'),]
data<-data[-which(rownames(data)=='C89T'),]

pheatmap(data[,1:14], cluster_rows=T, show_rownames=T,cluster_cols=T,show_colnames = T, 
         #annotation_col=coldata_all,
         #annotation_colors = ann_colors,na_col = "gray",
         scale = "none")

coldata<-matrix(ncol = 1,nrow = 106)
rownames(coldata)<-rownames(data)

idx_T_counts<-grep("T",rownames(data))
idx_P_counts<-grep("P",rownames(data))

coldata[idx_T_counts,1]<-"Tumor"
coldata[idx_P_counts,1]<-"NAT"

colnames(coldata)<-c("category")

library(readxl)

BC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BC")
BC_cd<-as.data.frame(cbind(BC_clinicial_data$Patient,BC_clinicial_data$MIBC,BC_clinicial_data$AGE,BC_clinicial_data$SEX,
                           BC_clinicial_data$Tumor_size,BC_clinicial_data$`T`,BC_clinicial_data$Tumor_number,BC_clinicial_data$Grade))
colnames(BC_cd)<-c("Patients","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade')


clinicial_design<-rbind(BC_cd)

coldata<-as.data.frame(coldata)
coldata$sample_num<-rownames(coldata)
coldata$Patients<-rownames(coldata)
coldata$Patients <- sub(paste0("*", "T"), "", coldata$Patients)
coldata$Patients <- sub(paste0("*", "P"), "", coldata$Patients)


coldata_all<-merge(coldata,clinicial_design,by="Patients",all.x=T)
rownames(coldata_all)<-coldata_all$sample_num
coldata_all<-coldata_all[,-1]
coldata_all<-coldata_all[,-3]
coldata_all<-coldata_all[,-2]

t_coldata_all<-as.data.frame(t(coldata_all))
t_coldata_all <- t_coldata_all %>% dplyr::select(rownames(data))

coldata_all<-as.data.frame(t(t_coldata_all))
coldata_all$Tumor_size<-as.numeric(coldata_all$Tumor_size)
coldata_all$Age<-as.numeric(coldata_all$Age)

idx_2<-grep('2',coldata_all$StageT)
coldata_all$StageT[idx_2]<-'2'


idx_4<-grep('4',coldata_all$StageT)
coldata_all$StageT[idx_4]<-'4'

idx_single<-grep('1',coldata_all$Tumor_number)
coldata_all$Tumor_number[idx_single]<-'single'

idx_multiple<-grep('2',coldata_all$Tumor_number)
coldata_all$Tumor_number[idx_multiple]<-'multiple'

ann_colors <- list(
  Sex = c("Male" = "#406A93", "Female" = "#EFA9AE"),
  MIBC = c("MIBC" = "#65AE00", "NMIBC" = "#D3B356"),
  category = c("Tumor"="#DB498E","NAT"="#379FB4"),
  Tumor_size=colorRampPalette(c("gray", "#2C1714"))(100),
  Age = colorRampPalette(c("#EEE0CC", "#91191E"))(100),
  StageT=c('4'="#FF3300",'3'='#FF6633','2'='#FF9966','1'='#FFCC66','a'='#FFFF00'),
  Grade=c('Low'='#7A84D0','High'="#D78203"),
  sample=c('T'='#8E5D2A','NAT'="#64BE52"),
  Tumor_number=c('single'='#90D4A4','multiple'="#C49406")
)

data_heatmap<-as.matrix(t(data))
data_heatmap<-data_heatmap[1:14,]

library(Cairo)
#CairoPNG(file = "tissue_cibersort.png", width = 12, height = 5,res=600,units='in',family='Arial')
CairoPDF(file = "tissue_cibersort.pdf", width = 12, height = 5,family='Arial')

pheatmap(data_heatmap,
         cluster_rows = TRUE,
         show_rownames = T,
         cluster_cols = T,
         show_colnames = F,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = coldata_all, 
         annotation_colors = ann_colors, 
         fontsize_row =10,
         fontsize = 10,
         na_col = "gray",
         scale = "none",
         #color = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
         color = colorRampPalette(c("blue", "white", "red"))(100),
         legend=TRUE,
         name = 'Z-score')

dev.off()

data<-data[,1:14]
# 给 data 添加 category 信息
data$category <- coldata_all$category[match(rownames(data), rownames(coldata_all))]

# 转为长表
data_long <- melt(data, id.vars = "category", variable.name = "CellType", value.name = "Proportion")

# 创建 PDF 保存所有图
pdf("CellType_Proportion_Tumor_vs_NAT.pdf", width = 5, height = 5)

cell_types <- unique(data_long$CellType)

library(reshape2)
library(ggplot2)
library(ggpubr)
library(dplyr)

# 给 data 添加 category 信息
data$category <- coldata_all$category[match(rownames(data), rownames(coldata_all))]

# 转为长表
data_long <- melt(data, id.vars = "category", variable.name = "CellType", value.name = "Proportion")

cell_types <- unique(data_long$CellType)

for (cell in cell_types) {
  # 子集数据
  df <- data_long %>% filter(CellType == cell)
  # Wilcoxon 测试
  p_val <- wilcox.test(Proportion ~ category, data = df)$p.value
  p_val_label <- ifelse(p_val < 0.001, "***",
                        ifelse(p_val < 0.01, "**",
                               ifelse(p_val < 0.05, "*", "ns")))
  # 文件名安全处理
  file_name <- paste0("CellType_", gsub("[^A-Za-z0-9]", "_", cell), ".pdf")
  # 创建 PDF
  pdf(file_name, width = 3, height = 3)
  # 绘图
  p <- ggplot(df, aes(x = category, y = Proportion, fill = category)) +
    geom_boxplot(outliers = F,width = 0.6) +
    #geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
    scale_fill_manual(values = c("Tumor"="#DB498E", "NAT"="#379FB4")) +
    theme_classic() +
    ylab(cell) +
    xlab("") +
    annotate("text", x = 1.5, y = max(df$Proportion)*1.05, label = p_val_label, size = 6)
  print(p)
  dev.off()
}
