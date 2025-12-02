#####GTEx data 


setwd("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/cell_type_marker/database/GTEx")

folder_path <- "/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/cell_type_marker/database/GTEx"
list <-list.files(folder_path, pattern = "\\.gct$", full.names = FALSE)
list<-list[grep('tpm',list)]
all_mean_TPM<-as.data.frame(matrix(ncol=1))
colnames(all_mean_TPM)<-'ENSG'
for (file in list){
#file=list[1]
GTEx_tpm_data <- read.table(paste0("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/cell_type_marker/database/GTEx/",file,sep=''),skip=2,header=T,row.names = 1)
modified_data<-apply(GTEx_tpm_data[,-c(1:2)],1,mean)
modified_data_1<-cbind(GTEx_tpm_data[,c(1)],modified_data)
tissue_name <- sub(".*_v[0-9]+_(.*)\\.gct", "\\1",file )
colnames(modified_data_1)<-c("ENSG",tissue_name)
modified_data_1<-as.data.frame(modified_data_1)
all_mean_TPM<-merge(all_mean_TPM,modified_data_1,by='ENSG',all=T)
print(tissue_name)
}

all_mean_TPM_1<-all_mean_TPM[-56201,]
write.csv(all_mean_TPM_1,"GTEx_mean_TPM.csv")

all_mean_TPM_1<-read.csv("GTEx_mean_TPM.csv",header=T,row.names = 1)
#GTEx_tpm_data <- read.table(paste0("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/cell_type_marker/database/GTEx/",file,sep=''),skip=2,header=T,row.names = 1)


dataforanalysis<-all_mean_TPM_1
list <-list.files(folder_path, pattern = "\\.gct$", full.names = FALSE)
list<-list[grep('tpm',list)]
file=list[1]
GTEx_tpm_data <- read.table(paste0("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/cell_type_marker/database/GTEx/",file,sep=''),skip=2,header=T,row.names = 1)

symbols<-GTEx_tpm_data[,c(1:2)]

colnames(symbols)[1]<-'ENSG'
data_annotation<-merge(dataforanalysis,symbols,by='ENSG',all.x=T)
write.csv(data_annotation,"GTEx_mean_TPM_anno.csv")


colnames(data_annotation)
data_annotation[, 2:55] <- apply(data_annotation[, 2:55], 2, as.numeric)
Anno_count <- aggregate(data_annotation[,2:55], by=list(data_annotation$Description), FUN=sum)
write.csv(Anno_count,"GTEx_mean_TPM_anno_by_genename.csv")
