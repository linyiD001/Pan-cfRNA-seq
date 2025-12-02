# 安装和加载包

#if(!require(GSVA)) install.packages("BiocManager"); BiocManager::install("GSVA")
library(GSVA)
library(readxl)
CPM_cfRNA <- read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/cpm_matrix_from_star.csv",header=T,row.names = 1)
CPM_tissue <- read.csv("/dssg/home/acct-dahan/share/BC/tissue/RESULT/NORMALIZATION/cpm_matrix_from_star.csv",header=T,row.names = 1)
BC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BC")
clinical <- as.data.frame(cbind(BC_clinicial_data$Patient,BC_clinicial_data$MIBC,BC_clinicial_data$AGE,BC_clinicial_data$SEX,
                                BC_clinicial_data$Tumor_size,BC_clinicial_data$`T`,BC_clinicial_data$Tumor_number,BC_clinicial_data$Grade))
colnames(clinical)<-c("Sample","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade')

clinical$StageT[grep('1',clinical$StageT)]<-'T1_2'
clinical$StageT[grep('2',clinical$StageT)]<-'T1_2'
clinical$StageT[grep('3',clinical$StageT)]<-'T3_4'
clinical$StageT[grep('4',clinical$StageT)]<-'T3_4'
clinical$StageT[which(clinical$StageT=='a')]<-'Ta'

idx_T<-grep("T",colnames(CPM_tissue))
T_CPM_tissue<-CPM_tissue[,idx_T]
rownames(T_CPM_tissue)<-rownames(CPM_tissue)

idx_P<-grep("P",colnames(CPM_tissue))
P_CPM_tissue<-CPM_tissue[,idx_P]
rownames(P_CPM_tissue)<-rownames(CPM_tissue)

idx_N<-grep("N",colnames(CPM_cfRNA))
N_CPM_cfRNA<-CPM_cfRNA[,idx_N]
rownames(N_CPM_cfRNA)<-rownames(CPM_cfRNA)

idx_X<-grep("X",colnames(CPM_cfRNA))
X_CPM_cfRNA<-CPM_cfRNA[,idx_X]
rownames(X_CPM_cfRNA)<-rownames(CPM_cfRNA)

colnames(T_CPM_tissue) <- sub(paste0("*", "T"), "", colnames(T_CPM_tissue))
colnames(P_CPM_tissue) <- sub(paste0("*", "P"), "", colnames(P_CPM_tissue))
colnames(N_CPM_cfRNA) <- sub(paste0("*", "_N"), "", colnames(N_CPM_cfRNA))
#colnames(X_CPM_cfRNA) <- sub(paste0("*", "_X"), "", colnames(X_CPM_cfRNA))

all_paired_list<-Reduce(intersect, list(colnames(T_CPM_tissue),
                                        colnames(P_CPM_tissue),
                                        colnames(N_CPM_cfRNA)
                                        #colnames(X_CPM_cfRNA)
))


all_paired_list<-as.data.frame(all_paired_list)[,1]
#write.table(all_paired_list,"all_4_paired_list.txt",col.names = F,row.names=F,quote = F)

###24 all paired samples
#share_gene<-intersect(rownames(T_CPM_tissue),rownames(N_CPM_cfRNA))
#share_gene<-intersect(share_gene,rownames(P_CPM_tissue))

paired_T_CPM_tissue_log2<-log2(T_CPM_tissue[,which(colnames(T_CPM_tissue) %in% all_paired_list)]+1)
paired_P_CPM_tissue_log2<-log2(P_CPM_tissue[,which(colnames(P_CPM_tissue) %in% all_paired_list)]+1)
#paired_X_CPM_cfRNA_log2<-log2(X_CPM_cfRNA[,which(colnames(X_CPM_cfRNA) %in% all_paired_list)]+1)
paired_N_CPM_cfRNA_log2<-log2(N_CPM_cfRNA[,which(colnames(N_CPM_cfRNA) %in% all_paired_list)]+1)

#if(!require(msigdbr)) install.packages("msigdbr")
library(msigdbr)
library(dplyr)

hallmark <- msigdbr(species = "Homo sapiens", category = "H")
gene_sets_hallmark <- split(x = hallmark$gene_symbol, f = hallmark$gs_name)
immune <- msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB")
gene_sets_immune <- split(immune$gene_symbol, immune$gs_name)

expr_matrix <- as.matrix(paired_N_CPM_cfRNA_log2)
gsvaPar_N_H <- ssgseaParam(expr_matrix,gene_sets_hallmark,  minSize = 1,
                           maxSize = Inf,
                           alpha = 0.25,
                           normalize = TRUE)
gsva_N_H<-gsva(gsvaPar_N_H, verbose=FALSE)

#expr_matrix <- as.matrix(paired_X_CPM_cfRNA_log2)
#gsvaPar_X_H <- ssgseaParam(expr_matrix,gene_sets_hallmark,  minSize = 1,
#                           maxSize = Inf,
#                           alpha = 0.25,
#                           normalize = TRUE)
#gsva_X_H<-gsva(gsvaPar_X_H, verbose=FALSE)

expr_matrix <- as.matrix(paired_T_CPM_tissue_log2)
gsvaPar_T_H <- ssgseaParam(expr_matrix, gene_sets_hallmark,  minSize = 1,
                           maxSize = Inf,
                           alpha = 0.25,
                           normalize = TRUE)
gsva_T_H<-gsva(gsvaPar_T_H, verbose=FALSE)

expr_matrix <- as.matrix(paired_P_CPM_tissue_log2)
gsvaPar_P_H <- ssgseaParam(expr_matrix, gene_sets_hallmark,  minSize = 1,
                           maxSize = Inf,
                           alpha = 0.25,
                           normalize = TRUE)
gsva_P_H<-gsva(gsvaPar_P_H, verbose=FALSE)


expr_matrix <- as.matrix(paired_N_CPM_cfRNA_log2)
gsvaPar_N_I <- ssgseaParam(expr_matrix,gene_sets_immune,  minSize = 1,
                           maxSize = Inf,
                           alpha = 0.25,
                           normalize = TRUE)
gsva_N_I<-gsva(gsvaPar_N_I, verbose=FALSE)

#expr_matrix <- as.matrix(paired_X_CPM_cfRNA_log2)
#gsvaPar_X_I <- ssgseaParam(expr_matrix,gene_sets_immune,  minSize = 1,
#                           maxSize = Inf,
#                           alpha = 0.25,
#                           normalize = TRUE)
#gsva_X_I<-gsva(gsvaPar_X_I, verbose=FALSE)

expr_matrix <- as.matrix(paired_T_CPM_tissue_log2)
gsvaPar_T_I <- ssgseaParam(expr_matrix, gene_sets_immune,  minSize = 1,
                           maxSize = Inf,
                           alpha = 0.25,
                           normalize = TRUE)
gsva_T_I<-gsva(gsvaPar_T_I, verbose=FALSE)

expr_matrix <- as.matrix(paired_P_CPM_tissue_log2)
gsvaPar_P_I <- ssgseaParam(expr_matrix,gene_sets_immune,  minSize = 1,
                           maxSize = Inf,
                           alpha = 0.25,
                           normalize = TRUE)
gsva_P_I<-gsva(gsvaPar_P_I, verbose=FALSE)

# 所有通路名称
pathways <- names(gene_sets_hallmark)  # 或 names(gene_sets_hallmark)

# 找出公共样本
common_samples <- intersect(colnames(gsva_T_H), colnames(gsva_N_H))

# 初始化存储结果的数据框
cor_results <- data.frame(
  Pathway = character(),
  SpearmanR = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

# 循环每条通路
for (pw in pathways) {
  df_tmp <- data.frame(
    Tumor = gsva_T_H[pw, common_samples],
    Urine_cfRNA = gsva_N_H[pw, common_samples]
  )
  
  cor_res <- cor.test(df_tmp$Tumor, df_tmp$Urine_cfRNA, method = "spearman")
  
  cor_results <- rbind(cor_results, data.frame(
    Pathway = pw,
    SpearmanR = round(cor_res$estimate, 2),
    P_value = signif(cor_res$p.value, 3)
  ))
}

# 查看前几行
head(cor_results)
cor_Tumor_N <- cor_results

pathways <- names(gene_sets_hallmark)  # 或 names(gene_sets_hallmark)
# 找出公共样本
common_samples <- intersect(colnames(gsva_P_H), colnames(gsva_N_H))
cor_results <- data.frame(
  Pathway = character(),
  SpearmanR = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)
# 循环每条通路
for (pw in pathways) {
  df_tmp <- data.frame(
    NAT = gsva_P_H[pw, common_samples],
    Urine_cfRNA = gsva_N_H[pw, common_samples]
  )
  
  cor_res <- cor.test(df_tmp$NAT, df_tmp$Urine_cfRNA, method = "spearman")
  
  cor_results <- rbind(cor_results, data.frame(
    Pathway = pw,
    SpearmanR = round(cor_res$estimate, 2),
    P_value = signif(cor_res$p.value, 3)
  ))
}
cor_NAT_N <- cor_results


#write.csv(cor_NAT_N,"cor_NAT_N.csv")
#write.csv(cor_Tumor_N,"cor_Tumor_N.csv")

###以每个样本进行
# ===== Step 1: 定义比较函数 =====
compare_sources <- function(source_mat, urine_mat, label) {
  # overlap patients
  Patients <- intersect(colnames(source_mat), colnames(urine_mat))
  # 计算 Cor 与 Spearman correlation
  Cor <- sapply(Patients, function(gs) {
    cor(source_mat[, gs], urine_mat[, gs], method = "spearman")
  })
  Pvalue <- sapply(Patients, function(gs) {
    cor.test(source_mat[, gs], urine_mat[, gs], method = "spearman")$p.value
  })
  # 返回三个变量：Patient, R2, Spearman
  data.frame(
    Patients = Patients,
    Cor = Cor,
    Pvalue =Pvalue,
    Source = label
  )
}
# ===== Step 2: 确保矩阵是 numeric matrix =====
gsva_N_H <- as.matrix(gsva_N_H)
#gsva_X_H <- as.matrix(gsva_X_H)
gsva_T_H <- as.matrix(gsva_T_H)
gsva_P_H <- as.matrix(gsva_P_H)
gsva_N_I <- as.matrix(gsva_N_I)
#gsva_X_I <- as.matrix(gsva_X_I)
gsva_T_I <- as.matrix(gsva_T_I)
gsva_P_I <- as.matrix(gsva_P_I)

# ===== Step 3: 对不同来源计算 R² =====
#res_blood_H <- compare_sources(gsva_X_H, gsva_N_H, "Blood_H")       # 血液 -> 尿液
res_tumor_H <- compare_sources(gsva_T_H, gsva_N_H, "Tumor_H")       # 肿瘤 -> 尿液
res_normal_H <- compare_sources(gsva_P_H, gsva_N_H, "NAT_H")   # 癌旁 -> 尿液

#res_blood_I <- compare_sources(gsva_X_I, gsva_N_I, "Blood_I")       # 血液 -> 尿液
res_tumor_I <- compare_sources(gsva_T_I, gsva_N_I, "Tumor_I")       # 肿瘤 -> 尿液
res_normal_I <- compare_sources(gsva_P_I, gsva_N_I, "NAT_I")   # 癌旁 -> 尿液

res_all<-rbind(res_tumor_H,res_normal_H,res_tumor_I,res_normal_I)
res_all_mean <-cbind(res_tumor_H[,2],res_normal_H[,2],res_tumor_I[,2],res_normal_I[,2])
cor_list_mean <- colMeans(res_all_mean,na.rm = TRUE)
cor_list_mean
# ===== Step 5: 设置横轴顺序 =====
# 假设 Source 列中有 "Tumor_H", "NAT_H", "Tumor_I", "NAT_I"
# 这里可以按你想要的顺序设置
res_all$Source <- factor(res_all$Source, levels = c("Tumor_H", "NAT_H", "Tumor_I", "NAT_I"))


library(ggpubr)
# ===== Step 6: 可视化 R² 分布 + 配对比较 =====
ggplot(res_all, aes(x = Source, y = Cor, fill = Source)) +
  geom_boxplot(alpha = 0.7,outliers = F) +
  #geom_jitter(width = 0.2, alpha = 0.8, size = 1.5) +
  theme_bw(base_size = 14) +
  labs(
    y = "Correlation between bladder tissue and urine",
    x = "",
    #title = "Explained variance of urinary ssGSEA by different sources"
  ) +
  scale_fill_manual(values = c(
    "Tumor_H"="pink", "NAT_H"="lightblue",
    "Tumor_I"="pink", "NAT_I"="lightblue"
  )) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(
    comparisons = list(
      c("Tumor_H", "NAT_H"),
      c("Tumor_I", "NAT_I")
    ),
    test = "wilcox.test",
    map_signif_level = TRUE
  )+ylim(0,1)

ggsave("Tumor_NAT_Urine.png",height=6,width=8,unit='in')
ggsave("Tumor_NAT_Urine.pdf",height=6,width=8,unit='in')

# ===== Step 1: 指定样本名 =====
sample_id <- "C108"  # ← 改成你要看的样本名

# ===== Step 2: 提取该样本在不同来源的所有 ssGSEA 分数 =====
df_plot <- data.frame(
  Geneset = rownames(gsva_N_H),
  Urine = gsva_N_H[,sample_id],
  Blood = gsva_X_H[,sample_id],
  Tumor = gsva_T_H[,sample_id],
  NAT = gsva_P_H[,sample_id]
)

# ===== Step 3: 转换为长格式方便绘图 =====
library(reshape2)
df_long <- melt(df_plot, id.vars = "Geneset", variable.name = "Source", value.name = "Score")

# ===== Step 4: 绘制散点图：每个来源 vs 尿液 =====
library(ggplot2)

ggplot(df_long, aes(x = Score, y = Geneset)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  facet_wrap(~ Source, nrow = 1) +
  theme_bw(base_size = 14) +
  labs(y = "ssGSEA score",
       title = paste("Sample:", sample_id, "- Relationship between urine and sources")) +
  scale_color_manual(values = c("Blood"="#DB498E", "Tumor"="#379FB4", "NAT"="#FCAE59")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_blank(),     # 隐藏y轴文字
    axis.ticks.y = element_blank()     # 同时去掉y轴刻度
  )

ggsave("C108_ssGESA.png",height=6,width=6)
ggsave("C108_ssGESA.pdf",height=6,width=6)

