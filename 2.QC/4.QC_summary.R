####organize all the sample list

library(dplyr) 
library(readxl)

setwd("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC")

BC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BC")
BS_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BS")
HC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "HC")
clinic <- bind_rows(BC_clinicial_data,BS_clinicial_data,HC_clinicial_data)

clinic$Plasma<-NA
clinic$Urine<-NA
clinic$Tumor<-NA
clinic$NAT<-NA

cfRNA_count<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/STAR/count.csv",header=T,row.names = 1)
tissue_count<-read.csv("/dssg/home/acct-dahan/share/tissue/RESULT/MAPPING/star/unanno_count.csv",header=T,row.names = 2)
tissue_count<-tissue_count[,-1]

idx_T<-grep("T",colnames(tissue_count))
T_tissue_count<-tissue_count[,idx_T]
colnames(T_tissue_count)<-sub(paste0("*", "T"), "", colnames(T_tissue_count))

idx_P<-grep("P",colnames(tissue_count))
P_tissue_count<-tissue_count[,idx_P]
colnames(P_tissue_count)<-sub(paste0("*", "P"), "", colnames(P_tissue_count))

idx_N<-grep("N",colnames(cfRNA_count))
N_tissue_count<-cfRNA_count[,idx_N]
colnames(N_tissue_count)<-sub(paste0("*", "_N"), "", colnames(N_tissue_count))

idx_X<-grep("X",colnames(cfRNA_count))
X_tissue_count<-cfRNA_count[,idx_X]
colnames(X_tissue_count)<-sub(paste0("*", "_X"), "", colnames(X_tissue_count))


clinic[which(clinic$Patient %in% colnames(X_tissue_count)),15]<-'plamsa_cfRNA'
clinic[which(clinic$Patient %in% colnames(N_tissue_count)),16]<-'urine_cfRNA'
clinic[which(clinic$Patient %in% colnames(T_tissue_count)),17]<-'Tumor_bulkRNA'
clinic[which(clinic$Patient %in% colnames(P_tissue_count)),18]<-'NAT_bulkRNA'

clinic<-as.data.frame(clinic)
idx_H<-grep("H",clinic$Patient)
clinic$group[idx_H]<-"Healthy Control"
idx_S<-grep("S",clinic$Patient)
clinic$group[idx_S]<-"Bladder Stone"
idx_C<-grep("C",clinic$Patient)
clinic$group[idx_C]<-"Bladder Cancer"

idx_R<-grep("C89",clinic$Patient)
clinic$group[idx_R]<-"R"
idx_R<-grep("C103",clinic$Patient)
clinic$group[idx_R]<-"R"
idx_R<-grep("C107",clinic$Patient)
clinic$group[idx_R]<-"R"

####add QC info
cfRNA_QC_list<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/QC_repeataware_list_withoutR.csv",header=T,row.names = 1)

idx_N<-grep("N",cfRNA_QC_list[,1])
passed_N_tissue<-cfRNA_QC_list[idx_N,1]
passed_N_tissue<-sub(paste0("*", "_N"), "", cfRNA_QC_list[idx_N,1])

idx_X<-grep("X",cfRNA_QC_list[,1])
passed_X_tissue<-cfRNA_QC_list[idx_X,1]
passed_X_tissue<-sub(paste0("*", "_X"), "", cfRNA_QC_list[idx_X,1])

clinic$Urine_cfRNA_QC<-NA
clinic$Plasma_cfRNA_QC<-NA

clinic[which(clinic$Patient %in% passed_N_tissue),20]<-'urine_cfRNA_QCpassed'
clinic[which(clinic$Patient %in% passed_X_tissue),21]<-'plasma_cfRNA_QCpassed'

write.csv(clinic,"sample_info_with_QC.csv")
