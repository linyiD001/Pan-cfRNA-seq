####separate dataset
rm(list=ls())

setwd("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/model/star_standard")
count_N <- read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/nor_N_count.csv",row.names = 1)
count_X <- read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/nor_X_count.csv", row.names = 1)

#提取样本主体名（去掉_N 或 _X）
get_base_name <- function(x) gsub("(_N|_X)$", "", x)
base_N <- sapply(colnames(count_N), get_base_name)
base_X <- sapply(colnames(count_X), get_base_name)

#取交集
common_base <- intersect(base_N, base_X)
cat("共有对应样本数:", length(common_base), "\n")

keep_N <- colnames(count_N)[base_N %in% common_base]
keep_X <- colnames(count_X)[base_X %in% common_base]
count_N <- count_N[, keep_N]
count_X <- count_X[, keep_X]

count_analysis<-as.data.frame.matrix(count_N)

coldata = data.frame(row.names = colnames(count_analysis))
idx_S<-grep('S',colnames(count_analysis))
coldata$type[idx_S]<-'2'
idx_H<-grep('H',colnames(count_analysis))
coldata$type[idx_H]<-'1'
idx_C<-grep('C',colnames(count_analysis))
coldata$type[idx_C]<-'3'


###cfRNA_filter
TPM<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/cpm_matrix_from_star.csv",header=T,row.names = 1)
N_TPM<-TPM[,grep("N",colnames(TPM))]
gene_mean <- rowMeans(N_TPM,na.rm = TRUE)
filtered_genes <- names(gene_mean[gene_mean > 10])
filtered_genes
count_analysis_DEG<-count_analysis[which(rownames(count_analysis) %in% filtered_genes),]

t_data <- as.data.frame(t(count_analysis_DEG))
data_for_model <- base::transform(t_data,outcome=coldata$type)

data_for_model$outcome<-as.factor(data_for_model$outcome)
set.seed(123)
##split data
############
library(caret)

train_indices <- createDataPartition(data_for_model$outcome, p = 0.8, list = FALSE)
train_data <- data_for_model[train_indices, ]
test_data <- data_for_model[-train_indices, ]

write.csv(train_data, "train_data_N.csv")
write.csv(test_data,  "test_data_N.csv")

#rm(list=ls())

####separate dataset
count_analysis<-as.data.frame.matrix(count_X)

coldata = data.frame(row.names = colnames(count_analysis))
idx_S<-grep('S',colnames(count_analysis))
coldata$type[idx_S]<-'2'
idx_H<-grep('H',colnames(count_analysis))
coldata$type[idx_H]<-'1'
idx_C<-grep('C',colnames(count_analysis))
coldata$type[idx_C]<-'3'


###cfRNA_filter
TPM<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/cpm_matrix_from_star.csv",header=T,row.names = 1)
X_TPM<-TPM[,grep("X",colnames(TPM))]
gene_mean <- rowMeans(X_TPM,na.rm = TRUE)

# 计算每个基因在多少个样本中表达为 0
zero_prop <- rowMeans(X_TPM == 0, na.rm = TRUE)

# 只保留在 ≤50% 样本中为 0 的基因
filtered_genes <- names(zero_prop[zero_prop <= 0.3])
count_analysis_DEG<-count_analysis[which(rownames(count_analysis) %in% filtered_genes),]

t_data <- as.data.frame(t(count_analysis_DEG))
data_for_model <- base::transform(t_data,outcome=coldata$type)

data_for_model$outcome<-as.factor(data_for_model$outcome)
train_indices <- createDataPartition(data_for_model$outcome, p = 0.8, list = FALSE)

train_data <- data_for_model[train_indices, ]
test_data <- data_for_model[-train_indices, ]


write.csv(train_data, "train_data_X.csv")
write.csv(test_data,  "test_data_X.csv")

