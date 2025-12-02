
load("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/DEG_starstandard.RData")


data <-count_N

coldata<-matrix(ncol = 2,nrow = 106)
rownames(coldata)<-colnames(data)

idx_N_counts<-grep("N",colnames(data))
idx_S_counts<-grep("S",colnames(data))
idx_H_counts<-grep("H",colnames(data))
idx_C_counts<-grep("C",colnames(data))

coldata[idx_N_counts,1]<-"Urine"
coldata[idx_H_counts,2]<-"HC"
coldata[idx_S_counts,2]<-"BS"
coldata[idx_C_counts,2]<-"BC"


colnames(coldata)<-c("category","Group")

library(readxl)

BC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BC")
BC_cd<-as.data.frame(cbind(BC_clinicial_data$Patient,BC_clinicial_data$MIBC,BC_clinicial_data$AGE,BC_clinicial_data$SEX,
                           BC_clinicial_data$Tumor_size,BC_clinicial_data$`T`,BC_clinicial_data$Tumor_number,BC_clinicial_data$Grade))
colnames(BC_cd)<-c("Patients","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade')

BS_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BS")
BS_cd<-as.data.frame(cbind(BS_clinicial_data$Patient,NA,BS_clinicial_data$AGE,BS_clinicial_data$SEX,NA,NA,NA,NA))
colnames(BS_cd)<-c("Patients","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade')

HC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "HC")
HC_cd<-as.data.frame(cbind(HC_clinicial_data$Patient,NA,HC_clinicial_data$AGE,HC_clinicial_data$SEX,NA,NA,NA,NA))
colnames(HC_cd)<-c("Patients","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade')


clinicial_design<-rbind(BC_cd,BS_cd,HC_cd)

coldata<-as.data.frame(coldata)
coldata$sample_num<-rownames(coldata)
coldata$Patients<-rownames(coldata)
coldata$Patients <- sub(paste0("*", "_X"), "", coldata$Patients)
coldata$Patients <- sub(paste0("*", "_N"), "", coldata$Patients)


coldata_all<-merge(coldata,clinicial_design,by="Patients",all.x=T)
rownames(coldata_all)<-coldata_all$sample_num
coldata_all<-coldata_all[,-1]
coldata_all<-coldata_all[,-3]

coldata_all <- coldata_all[colnames(data), , drop = FALSE]
dds <- DESeqDataSetFromMatrix(countData = round(data),
                              colData = coldata_all,
                              design = ~ Sex)
dds

smallestGroupSize <- 10
keep <- rowSums(counts(dds) > 5) >= smallestGroupSize
dds <- dds[keep,]
#dds$type <- relevel(dds$type, ref = "HC")

dds <- DESeq(dds)
resultsNames(dds)

ntd <- normTransform(dds)
library("pheatmap")
df <- as.data.frame(colData(dds)[,c("Sex")])
vsd <- vst(dds, blind=FALSE)
#rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Sex
colnames(sampleDistMatrix) <- colnames(count_analysis_N)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
#nor_count_N<-vst(dds, blind=T)


plotPCA(vsd,ntop=1000, intgroup=c("Sex"))

# 创建 PCA 数据
pcaData <- plotPCA(vsd,ntop=1000,  intgroup = c("Sex"), returnData = TRUE)

# 自定义颜色
colors <- c("Female" = "pink", "Male" = "lightblue")

# 绘制 PCA 图
ggplot(pcaData, aes(x = PC1, y = PC2, color = Sex)) +
  geom_point(size = 3) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  stat_ellipse(aes(fill = Sex, color = Sex), geom = "polygon", alpha = 0.15) +  # 添加半透明圆圈
  theme_bw() +
  labs(title = "PCA Plot",
       x = "PC1 49%",
       y = "PC2 9%") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank())
ggsave("PCA_for_gender.png",height=6,width=6,unit='in',dpi=300)
ggsave("PCA_for_gender.pdf",height=6,width=6)


library(umap)
# 提取归一化后的表达矩阵
mat <- assay(vsd)

# 运行 UMAP
set.seed(123)
topVarGenes <- head(order(rowVars(mat), decreasing=TRUE), 1000)
map_data <- umap(t(mat[topVarGenes, ]),
                 n_neighbors = 30,
                 min_dist = 0.3,
                 metric = "euclidean",
                 random_state = 123)
# 合并结果和样本信息
umap_df <- data.frame(map_data$layout, colData(vsd))
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")
umap_df$type <- coldata_all$Sex

# 画图
ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=type)) +
  geom_point(size=3, alpha=0.8) +
  theme_classic(base_size = 14) +
  labs(title="RNA-seq UMAP", color="Group") +
  theme(plot.title = element_text(hjust = 0.5))


data <-count_X

coldata<-matrix(ncol = 2,nrow = 117)
rownames(coldata)<-colnames(data)

idx_X_counts<-grep("X",colnames(data))
idx_S_counts<-grep("S",colnames(data))
idx_H_counts<-grep("H",colnames(data))
idx_C_counts<-grep("C",colnames(data))

coldata[idx_X_counts,1]<-"Plasma"
coldata[idx_H_counts,2]<-"HC"
coldata[idx_S_counts,2]<-"BS"
coldata[idx_C_counts,2]<-"BC"


colnames(coldata)<-c("category","Group")

library(readxl)

BC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BC")
BC_cd<-as.data.frame(cbind(BC_clinicial_data$Patient,BC_clinicial_data$MIBC,BC_clinicial_data$AGE,BC_clinicial_data$SEX,
                           BC_clinicial_data$Tumor_size,BC_clinicial_data$`T`,BC_clinicial_data$Tumor_number,BC_clinicial_data$Grade))
colnames(BC_cd)<-c("Patients","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade')

BS_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BS")
BS_cd<-as.data.frame(cbind(BS_clinicial_data$Patient,NA,BS_clinicial_data$AGE,BS_clinicial_data$SEX,NA,NA,NA,NA))
colnames(BS_cd)<-c("Patients","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade')

HC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "HC")
HC_cd<-as.data.frame(cbind(HC_clinicial_data$Patient,NA,HC_clinicial_data$AGE,HC_clinicial_data$SEX,NA,NA,NA,NA))
colnames(HC_cd)<-c("Patients","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade')


clinicial_design<-rbind(BC_cd,BS_cd,HC_cd)

coldata<-as.data.frame(coldata)
coldata$sample_num<-rownames(coldata)
coldata$Patients<-rownames(coldata)
coldata$Patients <- sub(paste0("*", "_X"), "", coldata$Patients)
coldata$Patients <- sub(paste0("*", "_N"), "", coldata$Patients)


coldata_all<-merge(coldata,clinicial_design,by="Patients",all.x=T)
rownames(coldata_all)<-coldata_all$sample_num
coldata_all<-coldata_all[,-1]
coldata_all<-coldata_all[,-3]

# 去除 Sex 为 NA 的样本
coldata_all <- coldata_all[!is.na(coldata_all$Sex), ]

# 同时在 countData 中也去掉对应列
data <- data[, rownames(coldata_all)]


coldata_all <- coldata_all[colnames(data), , drop = FALSE]
dds <- DESeqDataSetFromMatrix(countData = round(data),
                              colData = coldata_all,
                              design = ~ Sex)
dds

smallestGroupSize <- 10
keep <- rowSums(counts(dds) > 5) >= smallestGroupSize
dds <- dds[keep,]
#dds$type <- relevel(dds$type, ref = "HC")

dds <- DESeq(dds)
resultsNames(dds)

ntd <- normTransform(dds)
library("pheatmap")
df <- as.data.frame(colData(dds)[,c("Sex")])
vsd <- vst(dds, blind=FALSE)
#rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Sex
colnames(sampleDistMatrix) <- colnames(count_analysis_X)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
#nor_count_N<-vst(dds, blind=T)


plotPCA(vsd,ntop=1000, intgroup=c("Sex"))

# 创建 PCA 数据
pcaData <- plotPCA(vsd,ntop=1000,  intgroup = c("Sex"), returnData = TRUE)

# 自定义颜色
colors <- c("Female" = "pink", "Male" = "lightblue")

# 绘制 PCA 图
ggplot(pcaData, aes(x = PC1, y = PC2, color = Sex)) +
  geom_point(size = 3) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  stat_ellipse(aes(fill = Sex, color = Sex), geom = "polygon", alpha = 0.15) +  # 添加半透明圆圈
  theme_bw() +
  labs(title = "PCA Plot",
       x = "PC1 76%",
       y = "PC2 3%") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank())
ggsave("X_PCA_for_gender.png",height=6,width=6,unit='in',dpi=300)
ggsave("X_PCA_for_gender.pdf",height=6,width=6)
