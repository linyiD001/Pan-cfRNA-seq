# 安装和加载包
setwd("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/cell_type_marker/urine_plasma")
#if(!require(GSVA)) install.packages("BiocManager"); BiocManager::install("GSVA")
library(GSVA)

CPM_cfRNA <- read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/cpm_matrix_from_star.csv",header=T,row.names = 1)

idx_N<-grep("N",colnames(CPM_cfRNA))
N_CPM_cfRNA<-CPM_cfRNA[,idx_N]
rownames(N_CPM_cfRNA)<-rownames(CPM_cfRNA)

idx_X<-grep("X",colnames(CPM_cfRNA))
X_CPM_cfRNA<-CPM_cfRNA[,idx_X]
rownames(X_CPM_cfRNA)<-rownames(CPM_cfRNA)

colnames(N_CPM_cfRNA) <- sub(paste0("*", "_N"), "", colnames(N_CPM_cfRNA))
colnames(X_CPM_cfRNA) <- sub(paste0("*", "_X"), "", colnames(X_CPM_cfRNA))

all_paired_list<-Reduce(intersect, list(colnames(N_CPM_cfRNA),
                                        colnames(X_CPM_cfRNA)
))

all_paired_gene<-Reduce(intersect, list(colnames(N_CPM_cfRNA),
                                        colnames(X_CPM_cfRNA)
))


all_paired_list<-as.data.frame(all_paired_list)[,1]
#write.table(all_paired_list,"all_4_paired_list.txt",col.names = F,row.names=F,quote = F)


# 计算 N 组中高表达比例
expr_N <- rowMeans(N_CPM_cfRNA > 50, na.rm = F)
expr_X <- rowMeans(X_CPM_cfRNA > 50, na.rm = F)

# 保留在尿液中至少 80% 样本表达，且 T、P 组中无 NA 的基因
keep_genes <- (expr_N > 0.3) & (expr_X > 0.3)
keep_genes_names <- names(keep_genes[keep_genes == TRUE])

paired_X_CPM_cfRNA_log2<-log2(X_CPM_cfRNA[keep_genes_names,which(colnames(X_CPM_cfRNA) %in% all_paired_list)]+1)
paired_N_CPM_cfRNA_log2<-log2(N_CPM_cfRNA[keep_genes_names,which(colnames(N_CPM_cfRNA) %in% all_paired_list)]+1)

# 计算函数
calc_R2_P <- function(x, y) {
  if (length(unique(x)) <= 1 || length(unique(y)) <= 1) {
    return(c(R2 = NA, P = NA))
  }
  test <- suppressWarnings(cor.test(x, y))
  r <- as.numeric(test$estimate)
  R2 <- r^2
  P <- test$p.value
  return(c(R2 = R2,r=r,P = P))
}

# 初始化结果矩阵
genes <- keep_genes_names
res_NX<- data.frame(Gene = genes, R2 = NA, r=NA, P = NA)

# 遍历基因逐一计算
for (i in seq_along(genes)) {
  g <- genes[i]
  
  # Urine vs Plasma
  res_NX[i, 2:4] <- calc_R2_P(
    as.numeric(paired_N_CPM_cfRNA_log2[g, ]),
    as.numeric(paired_X_CPM_cfRNA_log2[g, ])
  )

}

head(res_NX)
# 保存结果
#write.csv(res_NX, "Urine_vs_Plasma_R2_P.csv", row.names = FALSE)

# 查看显著的相关基因
sig_NX <- subset(res_NX, P < 0.05)

# 如果想看相关性分布
library(ggplot2)
ggplot(res_NX, aes(x = R2)) +
  geom_histogram(bins = 30, fill = "#1f77b4", alpha = 0.8) +
  theme_bw() + ggtitle("Tumor–Urine cfRNA R² per genes")

ggplot(res_NX, aes(x = R2)) +
  geom_density(fill = "#1f77b4", alpha = 0.6) +
  theme_bw() +
  ggtitle("Urine–Plasma") +
  xlab("R²") +
  ylab("Density")

#ggsave("Tumor–Urine cfRNA R² per genes.png",height=6,width=6,unit='in')
#ggsave("Tumor–Urine cfRNA R² per genes.pdf",height=6,width=6,unit='in')


p1 <- ggplot(res_NX, aes(x = r, y = -log10(P))) +
  labs(x = expression('correlation'), 
       y = expression('-Log'['10']*'(pvalue)'), 
       title = NULL) +
  
  ## 在这里按条件设颜色
  geom_point(aes(color = r < 0.5),
             size = 2, shape = 16, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "skyblue", "TRUE" = "grey70")) +
  
  geom_hline(yintercept = -log10(0.05), color = 'darkgreen', linetype = 'dashed') +
  geom_vline(xintercept = 0.5, color = 'darkgreen', linetype = 'dashed') +
  
  geom_text_repel(
    data = subset(res_NX, (-log10(P) > -log10(0.05)) & abs(r) > 0.5),
    aes(label = Gene),
    nudge_x = 0.1, segment.alpha = 0.5,
    min.segment.length = 0.3, max.overlaps = 10, size = 3
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 10, margin = margin(t = 5), family='Arial'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )


p1

#ggsave("NX_N_X.png",height=6,width=6,units='in',dpi=1200)
#ggsave("NX_N_X.pdf",height=6,width=6,units='in',device = cairo_pdf,family='Arial')



library(ggplot2)

paired_X_CPM_cfRNA_log2<-log2(X_CPM_cfRNA[,which(colnames(X_CPM_cfRNA) %in% all_paired_list)]+1)
paired_N_CPM_cfRNA_log2<-log2(N_CPM_cfRNA[,which(colnames(N_CPM_cfRNA) %in% all_paired_list)]+1)

# 指定你感兴趣的基因名称
gene_of_interest <- "CCDC7"  # ← 这里改成你要画的基因名

####这里介入性别的信息
#
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

# 提取该基因在尿液与血浆中的 log2 表达量
df_plot <- data.frame(
  Urine = as.numeric(paired_N_CPM_cfRNA_log2[gene_of_interest, all_paired_list]),
  Plasma = as.numeric(paired_X_CPM_cfRNA_log2[gene_of_interest, all_paired_list])
)
df_plot$Patients <- sub(paste0("*", "_N"), "", colnames(paired_N_CPM_cfRNA_log2[, all_paired_list]))
df_plot<-merge(df_plot,clinicial_design,by="Patients")
p_gene <- ggplot(df_plot, aes(x = Plasma, y = Urine, color = Sex)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  theme_bw() +
  labs(
    title = paste0("Expression Correlation of ", gene_of_interest),
    x = "Plasma (log2 CPM + 1)",
    y = "Urine (log2 CPM + 1)",
    color = "Sex"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(size = 10, family = "Arial"),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  scale_color_manual(values = c("Male" = "#1f77b4", "Female" = "#e377c2"))

# 计算 Pearson 相关
cor_val <- cor(df_plot$Plasma, df_plot$Urine, method = "pearson")

p_gene <- p_gene + annotate(
  "text", x = Inf, y = Inf,
  label = paste0("r = ", round(cor_val, 2)),
  hjust = 1.2, vjust = 1.5,
  size = 4, color = "black"
)

p_gene

ggsave("CCDC7_gender.png",height=6,width=6)
ggsave("CCDC7_gender.pdf",height=6,width=6,device=cairo_pdf)

# 计算平均表达量
mean_N <- rowMeans(N_CPM_cfRNA, na.rm = TRUE)
mean_X <- rowMeans(X_CPM_cfRNA, na.rm = TRUE)

# 各自取前100个基因
top100_N <- names(sort(mean_N, decreasing = TRUE))[1:1000]
top100_N
top100_X <- names(sort(mean_X, decreasing = TRUE))[1:1000]
top100_X

# 查看结果
head(top100_N)
head(top100_X)

# 取并集或交集（可选）
top100_union <- union(top100_N, top100_X)   # 血+尿中任一高表达
top100_intersect <- intersect(top100_N, top100_X)  # 血和尿中都高表达

library(VennDiagram)

# 定义集合
venn.plot <- draw.pairwise.venn(
  area1 = length(top100_N),
  area2 = length(top100_X),
  cross.area = length(intersect(top100_N, top100_X)),
  category = c("Urine Top 1000 genes", "Plasma Top 1000 genes"),
  fill = c("#1f77b4", "#ff7f0e"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05)
)
grid.draw(venn.plot)
