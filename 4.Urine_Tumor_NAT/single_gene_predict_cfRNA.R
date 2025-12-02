##paired raw count cor gene level
library(dplyr)
library(VariableScreening)
library(biomaRt)
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

all_paired_gene<-Reduce(intersect, list(rownames(T_CPM_tissue),
                                        rownames(P_CPM_tissue),
                                        rownames(N_CPM_cfRNA)))

#all_paired_list<-as.data.frame(all_paired_list)
#write.table(all_paired_list,"all_3_paired_list.txt",col.names = F,row.names=F,quote = F)

#expr_10 <- rowMeans(N_CPM_cfRNA > 10)
#keep_genes <- (expr_10 >= 0.5)
#keep_genes_names <- names(keep_genes[keep_genes == TRUE])

N_C_H_DEG<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/N_type_C_vs_H_result.csv",header=T,row.names = 1)
padjthreshold<-0.05
log2FoldChangethreshold<-4
basemeanthreshold<-5
up_N_C_H_DEG<-rownames(N_C_H_DEG[which(N_C_H_DEG$padj<padjthreshold&N_C_H_DEG$log2FoldChange>log2FoldChangethreshold&N_C_H_DEG$baseMean>basemeanthreshold),])

#keep_genes_names <- c("IGF2","H19","IVL")
keep_genes_names <- up_N_C_H_DEG
keep_genes_names <- keep_genes_names[-grep("ENSG",keep_genes_names)]
keep_genes_names


gene_fit_df <- data.frame(
  Gene = character(),
  R2 = numeric(),
  P_overall = numeric(),       # 新增总体显著性P值
  Coef_Tumor = numeric(),
  Coef_Para = numeric(),
  P_Tumor = numeric(),
  P_Para = numeric(),
  stringsAsFactors = FALSE
)

paired_T_CPM_tissue<- T_CPM_tissue[which(rownames(T_CPM_tissue) %in% keep_genes_names),which(colnames(T_CPM_tissue) %in% all_paired_list)]
paired_P_CPM_tissue<- P_CPM_tissue[which(rownames(P_CPM_tissue) %in% keep_genes_names),which(colnames(P_CPM_tissue) %in% all_paired_list)]
paired_N_CPM_cfRNA <- N_CPM_cfRNA[which(rownames(N_CPM_cfRNA) %in% keep_genes_names),which(colnames(N_CPM_cfRNA) %in% all_paired_list)]
log2_paired_T_CPM_tissue <- log2(paired_T_CPM_tissue + 1)
log2_paired_P_CPM_tissue <- log2(paired_P_CPM_tissue + 1)
log2_paired_N_CPM_cfRNA <- log2(paired_N_CPM_cfRNA + 1)


for (g in keep_genes_names) {
  
  y <- as.numeric(log2_paired_N_CPM_cfRNA[g, ])
  x_T <- as.numeric(log2_paired_T_CPM_tissue[g, ])
  x_P <- as.numeric(log2_paired_P_CPM_tissue[g, ])
  
  df <- data.frame(Urine = y, Tumor = x_T, Para = x_P)
  
  # 忽略 NA 计算标准差
  # if (sd(y, na.rm = TRUE) == 0 || sd(x_T, na.rm = TRUE) == 0 || sd(x_P, na.rm = TRUE) == 0) next
  
  model <- lm(Urine ~ Tumor + Para, data = df)
  smry <- summary(model)
  
  coef_Tumor <- coef_Para <- p_Tumor <- p_Para <- P_overall <- NA
  
  if ("Tumor" %in% rownames(smry$coefficients)) {
    coef_Tumor <- smry$coefficients["Tumor", "Estimate"]
    p_Tumor <- smry$coefficients["Tumor", "Pr(>|t|)"]
  }
  if ("Para" %in% rownames(smry$coefficients)) {
    coef_Para <- smry$coefficients["Para", "Estimate"]
    p_Para <- smry$coefficients["Para", "Pr(>|t|)"]
  }
  
  if (!is.null(smry$fstatistic)) {
    fval <- smry$fstatistic[1]
    df1 <- smry$fstatistic[2]
    df2 <- smry$fstatistic[3]
    P_overall <- pf(fval, df1, df2, lower.tail = FALSE)
  }
  smry$r.squared
  R2 <- smry$r.squared
  
  gene_fit_df <- rbind(gene_fit_df, data.frame(
    Gene = g,
    R2 = R2,
    P_overall = P_overall,
    Coef_Tumor = coef_Tumor,
    Coef_Para = coef_Para,
    P_Tumor = p_Tumor,
    P_Para = p_Para
  ))
}
#gene_fit_df_select <- gene_fit_df[which(gene_fit_df$P_overall<0.05),]

ggplot(gene_fit_df, aes(x = R2)) +
  geom_density(fill = "#377EB8", alpha = 0.6) +
  geom_vline(aes(xintercept = mean(R2, na.rm = TRUE)), 
             color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = mean(gene_fit_df$R2, na.rm = TRUE), 
           y = 0, label = "Mean", color = "red", vjust = -1, hjust = 0) +
  labs(
    title = "Paired Tissue–Urine cfRNA Prediction",
    x = "R²",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

ggsave("paired_tissue_urine_single_cfRNA_prediction.png",height=6,width=12,dpi=600,unit='in')
ggsave("paired_tissue_urine_single_cfRNA_prediction.pdf",height=6,width=12)

#gene_fit_df<-read.csv("Paired_Prediction.csv",header=T,row.names = 1)
write.csv(gene_fit_df,"Paired_Prediction.csv")
