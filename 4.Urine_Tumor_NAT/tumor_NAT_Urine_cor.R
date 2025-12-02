##paired raw count cor gene level

setwd("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/cell_type_marker/paired_cor")

CPM_cfRNA<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/cpm_matrix_from_star.csv",header=T,row.names = 1)
CPM_tissue<-read.csv("/dssg/home/acct-dahan/share/BC/tissue/RESULT/NORMALIZATION/cpm_matrix_from_star.csv",header=T,row.names = 1)

idx_T<-grep("T",colnames(CPM_tissue))
T_CPM_tissue<-CPM_tissue[,idx_T]
rownames(T_CPM_tissue)<-rownames(CPM_tissue)

idx_P<-grep("P",colnames(CPM_tissue))
P_CPM_tissue<-CPM_tissue[,idx_P]
rownames(P_CPM_tissue)<-rownames(CPM_tissue)

idx_N<-grep("N",colnames(CPM_cfRNA))
N_CPM_cfRNA<-CPM_cfRNA[,idx_N]
rownames(N_CPM_cfRNA)<-rownames(CPM_cfRNA)

colnames(T_CPM_tissue) <- sub(paste0("*", "T"), "", colnames(T_CPM_tissue))
colnames(P_CPM_tissue) <- sub(paste0("*", "P"), "", colnames(P_CPM_tissue))
colnames(N_CPM_cfRNA) <- sub(paste0("*", "_N"), "", colnames(N_CPM_cfRNA))

all_paired_list<-Reduce(intersect, list(colnames(T_CPM_tissue),
                                        colnames(P_CPM_tissue),
                                        colnames(N_CPM_cfRNA)))

#all_paired_list<-as.data.frame(all_paired_list)
#write.table(all_paired_list,"all_3_paired_list.txt",col.names = F,row.names=F,quote = F)

# 计算每个基因在 N 组中表达的样本比例
# 计算每个基因在 T、P 中是否完全有值（无 NA）
#expr_T <- rowMeans(!is.na(T_CPM_tissue))
#expr_P <- rowMeans(!is.na(P_CPM_tissue))

# 计算 N 组中高表达比例
expr_10 <- rowMeans(N_CPM_cfRNA > 10, na.rm = F)

# 保留在尿液中至少 80% 样本表达，且 T、P 组中无 NA 的基因
keep_genes <- (expr_10 >= 0.5)
keep_genes_names <- names(keep_genes[keep_genes == TRUE])

all_paired_gene<-Reduce(intersect, list(keep_genes_names,
                                        rownames(T_CPM_tissue),
                                        rownames(P_CPM_tissue)))


#keep_genes_names <- up_N_C_H_DEG
#keep_genes_names <- keep_genes_names[-grep("ENSG",keep_genes_names)]
#keep_genes_names
###24 all paired samples
paired_T_CPM_tissue<-T_CPM_tissue[all_paired_gene,which(colnames(T_CPM_tissue) %in% all_paired_list)]
paired_P_CPM_tissue<-P_CPM_tissue[all_paired_gene,which(colnames(P_CPM_tissue) %in% all_paired_list)]

#keep_genes_names <- rownames(paired_T_CPM_tissue)
paired_N_CPM_cfRNA<-N_CPM_cfRNA[all_paired_gene,which(colnames(N_CPM_cfRNA) %in% all_paired_list)]

all.equal(rownames(paired_T_CPM_tissue),
          rownames(paired_P_CPM_tissue),
          rownames(paired_N_CPM_cfRNA))

all.equal(colnames(paired_T_CPM_tissue),
          colnames(paired_P_CPM_tissue),
          colnames(paired_N_CPM_cfRNA))

log2_paired_T_CPM_tissue <- log2(paired_T_CPM_tissue + 1)
log2_paired_P_CPM_tissue <- log2(paired_P_CPM_tissue + 1)
log2_paired_N_CPM_cfRNA <- log2(paired_N_CPM_cfRNA + 1)

#log2_paired_P_CPM_tissue <- na.omit(log2_paired_P_CPM_tissue)  # 去掉含NA的列
#log2_paired_T_CPM_tissue <- na.omit(log2_paired_T_CPM_tissue)  # 去掉含NA的列

j=5
# 数据准备
data_combine <- cbind(log2_paired_P_CPM_tissue[,j], log2_paired_N_CPM_cfRNA[,j])
colnames(data_combine) <- c("NAT", "Urine")
data_combine <- as.data.frame(data_combine)

# 计算相关性统计值
cor_test <- cor.test(data_combine$NAT, data_combine$Urine)
r_val <- cor_test$estimate
R2_val <- r_val^2
p_val <- cor_test$p.value

# 绘图
plot_test <- ggplot(data_combine, aes(x = NAT, y = Urine)) +
  geom_point(size = 2, shape = 16, alpha = 0.8, color = "#1f77b4") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("R² = ", round(R2_val, 3), 
                          "\nP = ", formatC(p_val, format = "e", digits = 2)),
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  theme_bw(base_size = 14) +
  labs(x = "NAT expression (log2 CPM)",
       y = "Urine cfRNA expression (log2 CPM)",
       title = "Tumor–Urine Expression Correlation")

plot_test


# 数据准备
data_combine <- cbind(log2_paired_T_CPM_tissue[,j], log2_paired_N_CPM_cfRNA[,j])
colnames(data_combine) <- c("Tumor", "Urine")
data_combine <- as.data.frame(data_combine)

# 计算相关性统计值
cor_test <- cor.test(data_combine$Tumor, data_combine$Urine)
r_val <- cor_test$estimate
R2_val <- r_val^2
p_val <- cor_test$p.value

# 绘图
plot_test <- ggplot(data_combine, aes(x = Tumor, y = Urine)) +
  geom_point(size = 2, shape = 16, alpha = 0.8, color = "#1f77b4") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("R² = ", round(R2_val, 3), 
                          "\nP = ", formatC(p_val, format = "e", digits = 2)),
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  theme_bw(base_size = 14) +
  labs(x = "Tumor expression (log2 CPM)",
       y = "Urine cfRNA expression (log2 CPM)",
       title = "Tumor–Urine Expression Correlation")

plot_test
# 假设三个矩阵都有相同的列名（基因名）且样本配对顺序一致
# log2_paired_T_CPM_tissue: 肿瘤组织
# log2_paired_P_CPM_tissue: 癌旁组织 (NAT)
# log2_paired_N_CPM_cfRNA:  尿液 cfRNA

# 计算函数
calc_R2_P <- function(x, y) {
  if (length(unique(x)) <= 2 || length(unique(y)) <= 2) {
    return(c(R2 = NA,r=NA ,P = NA))
  }
  test <- suppressWarnings(cor.test(x, y))
  r <- as.numeric(test$estimate)
  R2 <- r^2
  P <- test$p.value
  return(c(R2 = R2,r=r, P = P))
}

# 初始化结果矩阵
genes <- intersect(rownames(log2_paired_T_CPM_tissue), keep_genes_names)
res_Tumor <- data.frame(Gene = genes, R2 = NA,r=NA,  P = NA)
res_NAT   <- data.frame(Gene = genes, R2 = NA,r=NA , P = NA)

# 遍历基因逐一计算
for (i in seq_along(genes)) {
  g <- genes[i]
  
  # Tumor vs Urine
  res_Tumor[i, 2:4] <- calc_R2_P(
    as.numeric(log2_paired_T_CPM_tissue[g,]),
    as.numeric(log2_paired_N_CPM_cfRNA[g,])
  )
  
  # NAT vs Urine
  res_NAT[i, 2:4] <- calc_R2_P(
    as.numeric(log2_paired_P_CPM_tissue[g,]),
    as.numeric(log2_paired_N_CPM_cfRNA[g,])
  )
}

# 保存结果
write.csv(res_Tumor, "Tumor_vs_Urine_R2_P.csv", row.names = FALSE)
write.csv(res_NAT, "NAT_vs_Urine_R2_P.csv", row.names = FALSE)

# 查看显著的相关基因
sig_Tumor <- subset(res_Tumor, P < 0.05)
sig_NAT <- subset(res_NAT, P < 0.05)

cat("Tumor–Urine 显著基因数:", nrow(sig_Tumor), "\n")
cat("NAT–Urine 显著基因数:", nrow(sig_NAT), "\n")

# 如果想看相关性分布
library(ggplot2)
ggplot(res_Tumor, aes(x = R2)) +
  geom_histogram(bins = 30, fill = "#1f77b4", alpha = 0.8) +
  theme_bw() + ggtitle("Tumor–Urine cfRNA R² 分布")

ggplot(res_NAT, aes(x = R2)) +
  geom_histogram(bins = 30, fill = "#ff7f0e", alpha = 0.8) +
  theme_bw() + ggtitle("NAT–Urine cfRNA R² 分布")

data_combine_2<- as.data.frame(cbind(res_Tumor[,3],res_NAT[,3]))
# 计算相关性统计值
colnames(data_combine_2) <- c("Tumor_Urine","NAT_Urine")
cor_test <- cor.test(as.numeric(data_combine_2$Tumor_Urine),as.numeric(data_combine_2$NAT_Urine))
r_val <- cor_test$estimate
R2_val <- r_val^2
p_val <- cor_test$p.value

# 绘图
plot_test <- ggplot(data_combine_2, aes(x = Tumor_Urine, y = NAT_Urine)) +
  geom_point(size = 2, shape = 16, alpha = 0.8, color = "#1f77b4") +
  geom_smooth(method = "lm", se = F, color = "red", linetype = "dashed") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("R² = ", round(R2_val, 3), 
                          "\nP = ", formatC(p_val, format = "e", digits = 2)),
           hjust = 1.1, vjust = 1.1, size = 5, color = "black") +
  theme_bw(base_size = 14) +
  labs(x = "Correlation Tumor_Urine",
       y = "Correlation NAT_Urine",
       title = "Correlation")
plot_test

#tissue <- log2(paired_N_CPM_cfRNA + 1)

cor_list<-matrix(ncol = 3,nrow = 9742)
rownames(cor_list)<-rownames(paired_T_CPM_tissue)
colnames(cor_list)<-c("Tumor_NAT","Tumor_Urine","NAT_Urine")
p_list<-matrix(ncol = 3,nrow = 9742)
rownames(p_list)<-rownames(paired_T_CPM_tissue)
colnames(p_list)<-c("Tumor_NAT","Tumor_Urine","NAT_Urine")


for (k in 1:9742) {
  # 取出每个基因在不同组的表达
  t_T <- as.numeric(log2_paired_T_CPM_tissue[k, ])
  t_P <- as.numeric(log2_paired_P_CPM_tissue[k, ])
  t_N <- as.numeric(log2_paired_N_CPM_cfRNA[k, ])
  
  # T vs P
  if (length(unique(t_T)) > 3 && length(unique(t_P)) > 3) {
    result_TP <- suppressWarnings(cor.test(t_T, t_P, method = 'pearson'))
    cor_list[k, 1] <- result_TP$estimate
    p_list[k, 1] <- result_TP$p.value
  } else {
    cor_list[k, 1] <- NA
    p_list[k, 1] <- NA
  }
  
  # T vs N
  if (length(unique(t_T)) > 3 && length(unique(t_N)) > 3) {
    result_TN <- suppressWarnings(cor.test(t_T, t_N, method = 'pearson'))
    cor_list[k, 2] <- result_TN$estimate
    p_list[k, 2] <- result_TN$p.value
  } else {
    cor_list[k, 2] <- NA
    p_list[k, 2] <- NA
  }
  
  # P vs N
  if (length(unique(t_P)) > 3 && length(unique(t_N)) > 3) {
    result_PN <- suppressWarnings(cor.test(t_P, t_N, method = 'pearson'))
    cor_list[k, 3] <- result_PN$estimate
    p_list[k, 3] <- result_PN$p.value
  } else {
    cor_list[k, 3] <- NA
    p_list[k, 3] <- NA
  }
}

write.csv(cor_list,"NTP_gene_cor_list.csv")
write.csv(p_list,"NTP_gene_p_list.csv")

colSums(!is.na(cor_list))


cor_list<-read.csv("NTP_gene_cor_list.csv",header=T,row.names=1)
p_list<-read.csv("NTP_gene_p_list.csv",header=T,row.names=1)

data_for_plot <- cor_list
data_for_plot <- na.omit(data_for_plot)
p1 <- ggplot(data_for_plot, aes(x =Tumor_Urine, y = NAT_Urine)) +
  geom_point(size = 2, shape = 16, alpha = 0.8)    # <<<<<< 实心点，shape=16

p1

select_T <- cor_list[which(cor_list$Tumor_Urine>0.5&p_list$Tumor_Urine<0.05&(cor_list$NAT_Urine<0.3|p_list$NAT_Urine>0.05)),]
select_T <- na.omit(select_T)


log2_N <- log2_paired_N_CPM_cfRNA
log2_T <- log2_paired_T_CPM_tissue
head(log2_N)

genes_to_plot=''
log2_N_sub <- t(log2_N[genes_to_plot, ])
log2_T_sub <- t(log2_T[genes_to_plot, ])

df <- as.data.frame(cbind(log2_N_sub,log2_T_sub))
colnames(df)<- c("Urine","Tumor")

# 计算 Pearson 相关系数和 p 值
cor_test <- cor.test(df$Urine, df$Tumor, method = "pearson")
R_val <- round(cor_test$estimate, 2)
p_val <- signif(cor_test$p.value, 3)
label_text <- paste0("R = ", R_val, ", p = ", p_val)

# 绘图
p <- ggplot(df, aes(x = Urine, y = Tumor)) +
  geom_point(color = "#379FB4", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  annotate("text", x = min(df$Urine), y = max(df$Tumor), 
           label = label_text, hjust = 0, vjust = 1, size = 5) +
  theme_bw(base_size = 14) +
  labs(
    x = "log2(Urine cfRNA expression + 1)",
    y = "log2(Tumor expression + 1)",
    title = genes_to_plot
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  )


print(p)

#ggsave(paste0(genes_to_plot,"_Urine_Tumor.png"),height=6,width = 6,dpi=600,unit='in')
#ggsave(paste0(genes_to_plot,"_Urine_Tumor.pdf"),height=6,width = 6)

select_P <- cor_list[which(cor_list$NAT_Urine>0.5&p_list$NAT_Urine<0.05&(cor_list$Tumor_Urine<0.3|p_list$Tumor_Urine>0.05)),]
select_P <- na.omit(select_P)

intersect(rownames(select_P),N_C_H_up_rownames)
log2_N <- log2_paired_N_CPM_cfRNA
log2_P <- log2_paired_P_CPM_tissue
head(log2_N)

genes_to_plot='MYO9A'
log2_N_sub <- t(log2_N[genes_to_plot, ])
log2_P_sub <- t(log2_P[genes_to_plot, ])

df <- as.data.frame(cbind(log2_N_sub,log2_P_sub))
colnames(df)<- c("Urine","NAT")

# 计算 Pearson 相关系数和 p 值
cor_test <- cor.test(df$Urine, df$NAT, method = "pearson")
R_val <- round(cor_test$estimate, 2)
p_val <- signif(cor_test$p.value, 3)
label_text <- paste0("R = ", R_val, ", p = ", p_val)

# 绘图
p <- ggplot(df, aes(x = Urine, y = NAT)) +
  geom_point(color = "#379FB4", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  annotate("text", x = min(df$Urine), y = max(df$NAT), 
           label = label_text, hjust = 0, vjust = 1, size = 5) +
  theme_bw(base_size = 14) +
  labs(
    x = "log2(Urine cfRNA expression + 1)",
    y = "log2(NAT expression + 1)",
    title = genes_to_plot
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  )


print(p)

ggsave(paste0(genes_to_plot,"_Urine_NAT.png"),height=6,width = 6,dpi=600,unit='in')
ggsave(paste0(genes_to_plot,"_Urine_NAT.pdf"),height=6,width = 6)

N_C_H_DEG<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/N_type_C_vs_H_result.csv",header=T,row.names = 1)
N_S_H_DEG<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/N_type_S_vs_H_result.csv",header=T,row.names = 1)
N_C_S_DEG<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/N_type_C_vs_S_result.csv",header=T,row.names = 1)
#T_P_paired_DEG<-read.csv("/dssg/home/acct-dahan/share/BC/tissue/RESULT/NORMALIZATION/paired_type_T_vs_P_result.csv",header=T,row.names = 1)
#T_P_unpaired_DEG<-read.csv("/dssg/home/acct-dahan/share/BC/tissue/RESULT/NORMALIZATION/unpaired_type_T_vs_P_result.csv",header=T,row.names = 1)

#idx_noname<-which(T_P_unpaired_DEG$hgnc_symbol=='')
#T_P_unpaired_DEG[idx_noname,7]<-T_P_unpaired_DEG[idx_noname,1]
#idx_noname<-which(T_P_paired_DEG$hgnc_symbol=='')
#T_P_paired_DEG[idx_noname,7]<-T_P_paired_DEG[idx_noname,1]

padjthreshold<-0.05
log2FoldChangethreshold<-1
basemeanthreshold<-0
tissuelog2FoldChangethreshold<-log2(2)
N_C_H_up_rownames<-rownames(N_C_H_DEG[which(N_C_H_DEG$padj<padjthreshold&(N_C_H_DEG$log2FoldChange)>log2FoldChangethreshold&N_C_H_DEG$baseMean>basemeanthreshold),])
N_C_S_up_rownames<-rownames(N_C_S_DEG[which(N_C_S_DEG$padj<padjthreshold&(N_C_S_DEG$log2FoldChange)>log2FoldChangethreshold&N_C_S_DEG$baseMean>basemeanthreshold),])
N_C_H_down_rownames<-rownames(N_C_H_DEG[which(N_C_H_DEG$padj<padjthreshold&(N_C_H_DEG$log2FoldChange)<(-log2FoldChangethreshold)&N_C_H_DEG$baseMean>basemeanthreshold),])
N_C_S_down_rownames<-rownames(N_C_S_DEG[which(N_C_S_DEG$padj<padjthreshold&(N_C_S_DEG$log2FoldChange)<(-log2FoldChangethreshold)&N_C_S_DEG$baseMean>basemeanthreshold),])


#genelist<-rownames(N_C_H_DEG[which(N_C_H_DEG$baseMean>basemeanthreshold),])

data_for_plot<-cbind(cor_list[,c(2)],p_list[,c(2)])
rownames(data_for_plot)<-rownames(cor_list)
#data_for_plot<-data_for_plot[which(rownames(data_for_plot) %in% genelist),]
data_for_plot<-as.data.frame(data_for_plot)
colSums(!is.na(data_for_plot))

colnames(data_for_plot)<-c("correlation","p")
data_for_plot$ID<-rownames(data_for_plot)

data_for_plot$Significance <- 'NS'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_H_up_rownames &!(ID %in% N_C_S_up_rownames))] <- 'N_C_H_only_up'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_S_up_rownames &!(ID %in% N_C_H_up_rownames))] <- 'N_C_S_only_up'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_H_up_rownames & ID %in% N_C_S_up_rownames)] <- 'both_up'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_H_down_rownames &!(ID %in% N_C_S_up_rownames))] <- 'N_C_H_only_down'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_S_down_rownames &!(ID %in% N_C_H_up_rownames))] <- 'N_C_S_only_down'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_H_down_rownames & ID %in% N_C_S_up_rownames)] <- 'both_down'

dim(data_for_plot[which(data_for_plot$p<0.05&data_for_plot$correlation>0.3),])
dim(data_for_plot[which(data_for_plot$p<0.05&data_for_plot$correlation<0.3),])
#dim(data_for_plot[which(data_for_plot$p<0.05&data_for_plot$correlation>0.3&data_for_plot$Significance=='N_C_H_only_DEG'),])
#dim(data_for_plot[which(data_for_plot$p<0.05&data_for_plot$correlation>0.3&data_for_plot$Significance=='N_C_S_only_DEG'),])
#dim(data_for_plot[which(data_for_plot$p<0.05&data_for_plot$correlation>0.5&data_for_plot$Significance=='both_DEG'),])
#data_for_plot[which(data_for_plot$p<0.05&data_for_plot$correlation>0.3&data_for_plot$Significance=='N_C_H_only_DEG'),]
#data_for_plot[which(data_for_plot$p<0.05&data_for_plot$correlation>0.3&data_for_plot$Significance=='N_C_S_only_DEG'),]
#data_for_plot[which(data_for_plot$p<0.05&data_for_plot$correlation>0.3&data_for_plot$Significance=='both_DEG'),]

cols <- c(
  "NS" = "#B0B0B0",              # 灰色
  "N_C_H_only_up" = "#F39B7FB2",   # 红色
  "N_C_S_only_up" = "#DB7093",   # 橙色
  "both_up" = "#FF7F00",         # 深橙/橘色
  "N_C_H_only_down" = "#5E3C99", # 紫色
  "N_C_S_only_down" = "#3288BD", # 蓝色
  "both_down" = "darkblue"        # 深蓝/紫色
)
data_for_plot<- na.omit(data_for_plot)
p1 <- ggplot(data_for_plot, aes(x = correlation, y = -log10(p), color = Significance)) +
  labs(x = expression('correlation'), 
       y = expression('-Log'['10']*'(pvalue)'), 
       title = NULL) +
  scale_alpha_manual(values = c(
    "NS" = 0.3,               # 浅灰透明
    "N_C_H_only_up" = 1,
    "N_C_S_only_up" = 1,
    "both_up" = 1,
    "N_C_H_only_down" = 1,
    "N_C_S_only_down" = 1,
    "both_down" = 1
  )) +
  geom_point(size = 2, shape = 16, alpha = 0.8) +   # <<<<<< 实心点，shape=16
  scale_color_manual(values = cols) +
  #geom_vline(xintercept = c(-1, 1), color = 'darkgreen', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), color = 'darkgreen', linetype = 'dashed') +
  geom_text_repel(data = subset(data_for_plot, (-log10(p) > log10(0.05) & correlation > 0.3 & Significance != "NS")),
                  aes(label = ID), nudge_x = 0.1,
                  segment.alpha = 0.5, min.segment.length = 0.3,max.overlaps = 10,
                  size = 3) +
  geom_text_repel(data = subset(data_for_plot, correlation < (-0.3) & -log10(p) > log10(0.05) & Significance != "NS"),
                  aes(label = ID),
                  segment.alpha = 0.5, min.segment.length = 0.3, max.overlaps = 10,
                  size = 3) +
  theme_minimal()+theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10,family='Arial'),
    legend.text = element_text(size = 10,family='Arial'),
    axis.title = element_text(size = 10, margin = margin(t = 5),family='Arial'),
    axis.title.y = element_text(margin = margin(r = 0),family='Arial'),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   # 设置图例项的水平间距
    legend.spacing.y = unit(1, 'cm'),   # 设置图例项的垂直间距
    legend.key.spacing.y = unit(0.1, "cm"),    # 设置图例色块的高度
    legend.key.height = unit(0.5, "cm"),    # 设置图例色块的高度
    legend.key.width = unit(0.5, "cm"),     # 设置图例色块的宽度
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    # 图例背景和边框
    #legend.background = element_rect(fill = "white", color = "black", size = 0.5), # 图例背景和边框
    #legend.key = element_rect(color = "black", size = 0.5)  # 每个图例项的边框
  )
p1
ggsave("NTP_T_N.png",height=5,width=6,units='in',dpi=1200)
ggsave("NTP_T_N.pdf",height=5,width=6,units='in',device = cairo_pdf,family='Arial')


data_for_plot<-cbind(cor_list[,c(3)],p_list[,c(3)])
rownames(data_for_plot)<-rownames(cor_list)
#data_for_plot<-data_for_plot[which(rownames(data_for_plot) %in% genelist),]
colSums(!is.na(data_for_plot))

data_for_plot<-as.data.frame(data_for_plot)
colnames(data_for_plot)<-c("correlation","p")
data_for_plot$ID<-rownames(data_for_plot)
data_for_plot$Significance <- 'NS'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_H_up_rownames &!(ID %in% N_C_S_up_rownames))] <- 'N_C_H_only_up'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_S_up_rownames &!(ID %in% N_C_H_up_rownames))] <- 'N_C_S_only_up'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_H_up_rownames & ID %in% N_C_S_up_rownames)] <- 'both_up'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_H_down_rownames &!(ID %in% N_C_S_up_rownames))] <- 'N_C_H_only_down'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_S_down_rownames &!(ID %in% N_C_H_up_rownames))] <- 'N_C_S_only_down'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_H_down_rownames & ID %in% N_C_S_up_rownames)] <- 'both_down'

dim(data_for_plot[which(data_for_plot$p<0.05&data_for_plot$correlation>0.5),])
dim(data_for_plot[which(data_for_plot$p<0.05&data_for_plot$correlation<0.5),])
dim(data_for_plot[which(data_for_plot$p<0.05&data_for_plot$correlation>0.5&data_for_plot$Significance=='N_C_H_only_DEG'),])
dim(data_for_plot[which(data_for_plot$p<0.05&data_for_plot$correlation>0.5&data_for_plot$Significance=='N_C_S_only_DEG'),])
dim(data_for_plot[which(data_for_plot$p<0.05&data_for_plot$correlation>0.5&data_for_plot$Significance=='both_DEG'),])

cols <- c(
  "NS" = "#B0B0B0",              # 灰色
  "N_C_H_only_up" = "#F39B7FB2",   # 红色
  "N_C_S_only_up" = "#DB7093",   # 橙色
  "both_up" = "#FF7F00",         # 深橙/橘色
  "N_C_H_only_down" = "#5E3C99", # 紫色
  "N_C_S_only_down" = "#3288BD", # 蓝色
  "both_down" = "darkblue"        # 深蓝/紫色
)
p2 <- ggplot(data_for_plot, aes(x = correlation, y = -log10(p), color = Significance)) +
  labs(x = expression('correlation'), 
       y = expression('-Log'['10']*'(pvalue)'), 
       title = NULL) +
  scale_alpha_manual(values = c(
    "NS" = 0.3,               # 浅灰透明
    "N_C_H_only_up" = 1,
    "N_C_S_only_up" = 1,
    "both_up" = 1,
    "N_C_H_only_down" = 1,
    "N_C_S_only_down" = 1,
    "both_down" = 1
  )) +
  geom_point(size = 2, shape = 16, alpha = 0.8) +   # <<<<<< 实心点，shape=16
  scale_color_manual(values = cols) +
  #geom_vline(xintercept = c(-1, 1), color = 'darkgreen', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), color = 'darkgreen', linetype = 'dashed') +
  geom_text_repel(data = subset(data_for_plot, (-log10(p) > log10(0.05) & correlation > 0.3 & Significance != "NS")),
                  aes(label = ID), nudge_x = 0.1,
                  segment.alpha = 0.5, min.segment.length = 0.3,max.overlaps = 10,
                  size = 3, ) +
  geom_text_repel(data = subset(data_for_plot, correlation < (-0.3) & -log10(p) > log10(0.05) & Significance != "NS"),
                  aes(label = ID),
                  segment.alpha = 0.5, min.segment.length = 0.3, max.overlaps = 10,
                  size = 3) +
  theme_minimal()+theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10,family='Arial'),
    legend.text = element_text(size = 10,family='Arial'),
    axis.title = element_text(size = 10, margin = margin(t = 5),family='Arial'),
    axis.title.y = element_text(margin = margin(r = 0),family='Arial'),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   # 设置图例项的水平间距
    legend.spacing.y = unit(1, 'cm'),   # 设置图例项的垂直间距
    legend.key.spacing.y = unit(0.1, "cm"),    # 设置图例色块的高度
    legend.key.height = unit(0.5, "cm"),    # 设置图例色块的高度
    legend.key.width = unit(0.5, "cm"),     # 设置图例色块的宽度
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    # 图例背景和边框
    #legend.background = element_rect(fill = "white", color = "black", size = 0.5), # 图例背景和边框
    #legend.key = element_rect(color = "black", size = 0.5)  # 每个图例项的边框
  )
p2
ggsave("NTX_NAT_N.png",height=5,width=6,units='in',dpi=1200)
ggsave("NTX_NAT_N.pdf",height=5,width=6,units='in',device = cairo_pdf,family='Arial')


data_for_plot<-cbind(cor_list[,c(2)],p_list[,c(2)])
# 转换成数据框
rownames(data_for_plot) <- rownames(cor_list)
cor_df <- as.data.frame(data_for_plot)
colnames(cor_df) <- c("r", "p")
sel_genes <- cor_df[which(cor_df$r > 0.3 & cor_df$p < 0.05),]
# ===== 1️⃣ 筛选相关性显著的基因 =====
sig_genes <- rownames(sel_genes)
cat("筛选到显著正相关基因数量：", length(sig_genes), "\n")

goAll_up <- enrichGO(sig_genes,
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "ALL",
                     pAdjustMethod= 'BH',
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.1)  # 以上三种 

top10_go <- goAll_up@result %>%
  dplyr::arrange(p.adjust) %>%      
  dplyr::slice_head(n = 10) %>%      
  dplyr::mutate(GeneRatio = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))) %>%  # 转为数值
  dplyr::arrange(desc(GeneRatio))  
go_for_plot <- top10_go
go_for_plot$Description <- str_wrap(go_for_plot$Description, width = 30)
p2_ggplot <- ggplot(go_for_plot, aes(x = reorder(Description, GeneRatio), y = GeneRatio, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.7) +  # 使用点绘制
  coord_flip() +  # 翻转坐标轴
  scale_size_continuous(range = c(3, 10)) +  # 设置点的大小范围
  scale_color_gradient(low = "blue", high = "red") +  # p-value 渐变色
  theme_minimal() +
  labs(x = "Gene Ontology enrichment analysis", y = "Gene Ratio", size = "Count", color = "Adjusted p-value") +
  theme(
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size = 20, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   # 设置图例项的水平间距
    legend.spacing.y = unit(1, 'cm'),   # 设置图例项的垂直间距
    legend.key.spacing.y = unit(0.1, "cm"),    # 设置图例色块的高度
    legend.key.height = unit(0.5, "cm"),    # 设置图例色块的高度
    legend.key.width = unit(0.5, "cm"),     # 设置图例色块的宽度
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    # 图例背景和边框
    #legend.background = element_rect(fill = "white", color = "black", size = 0.5), # 图例背景和边框
    #legend.key = element_rect(color = "black", size = 0.5)  # 每个图例项的边框
  ) +
  theme(legend.position = "right")

p2_ggplot
ggsave("Top_Tumor_correlated_KEGG.png",height=6,width=6,dpi=600,unit='in')
ggsave("Top_Tumor_correlated_KEGG.pdf",height=6,width=6)


# ===== 2️⃣ 计算这些基因在尿液中的平均表达 =====
expr_urine <- paired_N_CPM_cfRNA[sig_genes, , drop = FALSE]
mean_expr_urine <- rowMeans(expr_urine)

# ===== 3️⃣ 计算这些基因在尿液样本中的占比 =====
total_urine_expr <- colSums(paired_N_CPM_cfRNA)  # 每个样本的总表达

# ===== 4️⃣ 整理结果表 =====
result_df <- data.frame(
  Gene = sig_genes,
  Correlation = cor_df[sig_genes, "r"],
  Pvalue = cor_df[sig_genes, "p"],
  Mean_CPM_in_Urine = mean_expr_urine
)

# 排序（相关性高的在前）
result_df <- result_df[order(-result_df$Correlation), ]

# 输出结果
write.csv(result_df, "NTP_Tumor_Significant_correlated_genes_urine_expression.csv", row.names = FALSE)

library(ggplot2)
library(reshape2)
library(dplyr)

# 选出 top 10 显著正相关的基因
top10 <- head(result_df, 15)
top10_genes <- top10$Gene
# 把行名转成一列（基因名）
top10_expr <- paired_N_CPM_cfRNA[rownames(paired_N_CPM_cfRNA) %in% top10_genes, ]

top10_expr_df <- as.data.frame(top10_expr)
top10_expr_df$Gene <- rownames(top10_expr_df)

# 转换为长格式
plot_df <- melt(top10_expr_df, id.vars = "Gene", variable.name = "Sample", value.name = "CPM")

# 合并相关性信息
plot_df <- merge(plot_df, top10[, c("Gene", "Correlation")], by = "Gene", all.x = TRUE)

# 按相关性排序
plot_df$Gene <- factor(plot_df$Gene, levels = top10$Gene[order(top10$Correlation)])

# 绘图
ggplot(plot_df, aes(x = Gene, y = log2(CPM + 1), fill = Correlation)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5, alpha = 0.9, color = "grey30") +
  scale_fill_gradient2(
    low = "#4575b4", mid = "white", high = "#d73027",
    midpoint = 0.5,  # 颜色中点，可根据相关性范围调整
    name = "Correlation (r)"
  ) +
  coord_flip() +
  labs(
    title = "Tumor–Urine Correlated Genes",
    #subtitle = "Urinary cfRNA expression (log2 CPM + 1)\nColor represents correlation with plasma cfRNA",
    x = "Gene",
    y = "Expression in Urine(log2(CPM + 1))"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "right"
  )
ggsave("TPN_top_cor_T_in_urine.png",height=6,width=6,dpi=600,unit='in')
ggsave("TPN_top_cor_T_in_urine.pdf",height=6,width=6)


# log2转换
log2_N <- log2(paired_N_CPM_cfRNA + 1)
log2_T <- log2(paired_T_CPM_tissue + 1)

# 初始化向量
sample_corr <- numeric(ncol(log2_N))

for (i in 1:ncol(log2_N)) {
  sample_corr[i] <- cor(log2_N[, i], log2_T[, i], method = "pearson")
}

# 查看每对样本相关性
sample_corr
summary(sample_corr)

# 可画直方图
# 直方图
hist(sample_corr, breaks = 10, main = "Paired Tumor-Urine cfRNA Correlation", 
     xlab = "Pearson r", col = "#377EB8")

# 保存为 PNG
png("NTP_Tumor_Urine_cfRNA_correlation_hist.png", width = 6, height = 6, units = "in", res = 300)
hist(sample_corr, breaks = 10, main = "Paired Tumor-Urine cfRNA Correlation", 
     xlab = "Pearson r", col = "#377EB8")
dev.off()

# 保存为 PDF
pdf("NTP_Tumor_Urine_cfRNA_correlation_hist.pdf", width = 6, height = 6)
hist(sample_corr, breaks = 10, main = "Paired Tumor-Urine cfRNA Correlation", 
     xlab = "Pearson r", col = "#377EB8")
dev.off()



data_for_plot<-cbind(cor_list[,c(3)],p_list[,c(3)])
# 转换成数据框
rownames(data_for_plot) <- rownames(cor_list)
cor_df <- as.data.frame(data_for_plot)
colnames(cor_df) <- c("r", "p")
sel_genes <- cor_df[which(cor_df$r > 0.3 & cor_df$p < 0.05),]
# ===== 1️⃣ 筛选相关性显著的基因 =====
sig_genes <- rownames(sel_genes)
cat("筛选到显著正相关基因数量：", length(sig_genes), "\n")
sel_genes
goAll_up <- enrichGO(sig_genes,
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "ALL",
                     pAdjustMethod= 'BH',
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.1)  # 以上三种 
top10_go <- goAll_up@result %>%
  dplyr::arrange(p.adjust) %>%      
  dplyr::slice_head(n = 10) %>%      
  dplyr::mutate(GeneRatio = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))) %>%  # 转为数值
  dplyr::arrange(desc(GeneRatio))  
go_for_plot <- top10_go
go_for_plot$Description <- str_wrap(go_for_plot$Description, width = 30)
p2_ggplot <- ggplot(go_for_plot, aes(x = reorder(Description, GeneRatio), y = GeneRatio, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.7) +  # 使用点绘制
  coord_flip() +  # 翻转坐标轴
  scale_size_continuous(range = c(3, 10)) +  # 设置点的大小范围
  scale_color_gradient(low = "blue", high = "red") +  # p-value 渐变色
  theme_minimal() +
  labs(x = "Gene Ontology enrichment analysis", y = "Gene Ratio", size = "Count", color = "Adjusted p-value") +
  theme(
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size = 20, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   # 设置图例项的水平间距
    legend.spacing.y = unit(1, 'cm'),   # 设置图例项的垂直间距
    legend.key.spacing.y = unit(0.1, "cm"),    # 设置图例色块的高度
    legend.key.height = unit(0.5, "cm"),    # 设置图例色块的高度
    legend.key.width = unit(0.5, "cm"),     # 设置图例色块的宽度
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    # 图例背景和边框
    #legend.background = element_rect(fill = "white", color = "black", size = 0.5), # 图例背景和边框
    #legend.key = element_rect(color = "black", size = 0.5)  # 每个图例项的边框
  ) +
  theme(legend.position = "right")

p2_ggplot
ggsave("Top_NAT_correlated_KEGG.png",height=6,width=6,dpi=600,unit='in')
ggsave("Top_NAT_correlated_KEGG.pdf",height=6,width=6)


# ===== 2️⃣ 计算这些基因在尿液中的平均表达 =====
expr_urine <- paired_N_CPM_cfRNA[sig_genes, , drop = FALSE]
mean_expr_urine <- rowMeans(expr_urine)

# ===== 3️⃣ 计算这些基因在尿液样本中的占比 =====
total_urine_expr <- colSums(paired_N_CPM_cfRNA)  # 每个样本的总表达

# ===== 4️⃣ 整理结果表 =====
result_df <- data.frame(
  Gene = sig_genes,
  Correlation = cor_df[sig_genes, "r"],
  Pvalue = cor_df[sig_genes, "p"],
  Mean_CPM_in_Urine = mean_expr_urine
)

# 排序（相关性高的在前）
result_df <- result_df[order(-result_df$Correlation), ]

# 输出结果
write.csv(result_df, "NTP_NAT_Significant_correlated_genes_urine_expression.csv", row.names = FALSE)

library(ggplot2)
library(reshape2)
library(dplyr)

# 选出 top 10 显著正相关的基因
top10 <- head(result_df, 15)
top10_genes <- top10$Gene
# 把行名转成一列（基因名）
top10_expr <- paired_N_CPM_cfRNA[rownames(paired_N_CPM_cfRNA) %in% top10_genes, ]

top10_expr_df <- as.data.frame(top10_expr)
top10_expr_df$Gene <- rownames(top10_expr_df)

# 转换为长格式
plot_df <- melt(top10_expr_df, id.vars = "Gene", variable.name = "Sample", value.name = "CPM")

# 合并相关性信息
plot_df <- merge(plot_df, top10[, c("Gene", "Correlation")], by = "Gene", all.x = TRUE)

# 按相关性排序
plot_df$Gene <- factor(plot_df$Gene, levels = top10$Gene[order(top10$Correlation)])

# 绘图
ggplot(plot_df, aes(x = Gene, y = log2(CPM + 1), fill = Correlation)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5, alpha = 0.9, color = "grey30") +
  scale_fill_gradient2(
    low = "#4575b4", mid = "white", high = "#d73027",
    midpoint = 0.5,  # 颜色中点，可根据相关性范围调整
    name = "Correlation (r)"
  ) +
  coord_flip() +
  labs(
    title = "Tumor–Urine Correlated Genes",
    #subtitle = "Urinary cfRNA expression (log2 CPM + 1)\nColor represents correlation with plasma cfRNA",
    x = "Gene",
    y = "Expression in Urine(log2(CPM + 1))"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "right"
  )
ggsave("TPN_top_cor_P_in_urine.png",height=6,width=6,dpi=600,unit='in')
ggsave("TPN_top_cor_P_in_urine.pdf",height=6,width=6)


# log2转换
log2_N <- log2(paired_N_CPM_cfRNA + 1)
log2_P <- log2(paired_P_CPM_tissue + 1)

# 初始化向量
sample_corr <- numeric(ncol(log2_N))

for (i in 1:ncol(log2_N)) {
  sample_corr[i] <- cor(log2_N[, i], log2_P[, i], method = "pearson")
}

# 查看每对样本相关性
sample_corr
summary(sample_corr)

# 可画直方图
# 直方图
hist(sample_corr, breaks = 10, main = "Paired NAT-Urine cfRNA Correlation", 
     xlab = "Pearson r", col = "#377EB8")

# 保存为 PNG
png("NTP_NAT_Urine_cfRNA_correlation_hist.png", width = 6, height = 6, units = "in", res = 300)
hist(sample_corr, breaks = 10, main = "Paired NAT-Urine cfRNA Correlation", 
     xlab = "Pearson r", col = "#377EB8")
dev.off()

# 保存为 PDF
pdf("NTP_NAT_Urine_cfRNA_correlation_hist.pdf", width = 6, height = 6)
hist(sample_corr, breaks = 10, main = "Paired NAT-Urine cfRNA Correlation", 
     xlab = "Pearson r", col = "#377EB8")
dev.off()

