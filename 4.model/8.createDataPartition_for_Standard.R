####separate dataset
rm(list=ls())

setwd("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/model/salmon_standard")
count_N<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_salmonstandard/nor_count_N.csv",row.names = 1)
#count_N<-count_N[,-1]
count_X<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_salmonstandard/nor_count_X.csv",row.names = 1)
#count_X<-count_X[,-1]

count_analysis<-as.data.frame.matrix(count_N)

coldata = data.frame(row.names = colnames(count_analysis))
idx_S<-grep('S',colnames(count_analysis))
coldata$type[idx_S]<-'2'
idx_H<-grep('H',colnames(count_analysis))
coldata$type[idx_H]<-'1'
idx_C<-grep('C',colnames(count_analysis))
coldata$type[idx_C]<-'3'


###cfRNA_filter
NCH_DEG<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_salmonstandard/N_type_C_vs_H_result.csv",header=T,row.names = 1)
gene_list_DEG<-subset(NCH_DEG,baseMean>10)
count_analysis_DEG<-count_analysis[which(rownames(count_analysis) %in% rownames(gene_list_DEG)),]

t_data <- as.data.frame(t(count_analysis_DEG))
data_for_model <- base::transform(t_data,outcome=coldata$type)

data_for_model$outcome<-as.factor(data_for_model$outcome)
set.seed(250805)
##split data
############
# Use createDataPartition to carry out stratified sampling
train_indices <- createDataPartition(data_for_model$outcome, p=0.8, list=FALSE)
train_data <- data_for_model[train_indices, ]
test_data <- data_for_model[-train_indices, ]

write.csv(train_data,"train_data_N.csv")
write.csv(test_data,"test_data_N.csv")

rm(list=ls())

####separate dataset


count_N<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_salmonstandard/nor_count_N.csv",row.names = 1)
#count_N<-count_N[,-1]
count_X<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_salmonstandard/nor_count_X.csv",row.names = 1)
#count_X<-count_X[,-1]

count_analysis<-as.data.frame.matrix(count_X)

coldata = data.frame(row.names = colnames(count_analysis))
idx_S<-grep('S',colnames(count_analysis))
coldata$type[idx_S]<-'2'
idx_H<-grep('H',colnames(count_analysis))
coldata$type[idx_H]<-'1'
idx_C<-grep('C',colnames(count_analysis))
coldata$type[idx_C]<-'3'


###cfRNA_filter
NCH_DEG<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_salmonstandard/salmon_repeataware_stranded_allgenome/X_type_C_vs_H_result.csv",header=T,row.names = 1)
gene_list_DEG<-subset(NCH_DEG,baseMean>10)
count_analysis_DEG<-count_analysis[which(rownames(count_analysis) %in% rownames(gene_list_DEG)),]

t_data <- as.data.frame(t(count_analysis_DEG))
data_for_model <- base::transform(t_data,outcome=coldata$type)

data_for_model$outcome<-as.factor(data_for_model$outcome)
set.seed(250805)
##split data
############
# Use createDataPartition to carry out stratified sampling
train_indices <- createDataPartition(data_for_model$outcome, p=0.8, list=FALSE)
train_data <- data_for_model[train_indices, ]
test_data <- data_for_model[-train_indices, ]

write.csv(train_data,"train_data_X.csv")
write.csv(test_data,"test_data_X.csv")
