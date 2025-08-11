##paired raw count cor gene level


setwd("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/cell_type_marker/paired_cor")

TPM_cfRNA<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.standard/TPM_all_anno_by_gene_name.csv",header=T,row.names = 2)
TPM_tissue<-read.csv("/dssg/home/acct-dahan/share/tissue/RESULT/MAPPING/salmon.standard/TPM_all_anno_by_gene_name.csv",header=T,row.names = 2)
filtered_name<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/QC_repeataware_list.csv",header=T,row.names = 1)

idx_T<-grep("T",colnames(TPM_tissue))
T_TPM_tissue<-TPM_tissue[,idx_T]
rownames(T_TPM_tissue)<-rownames(TPM_tissue)

idx_P<-grep("P",colnames(TPM_tissue))
P_TPM_tissue<-TPM_tissue[,idx_P]
rownames(P_TPM_tissue)<-rownames(TPM_tissue)

idx_N<-grep("N",colnames(TPM_cfRNA))
N_TPM_cfRNA<-TPM_cfRNA[,idx_N]
rownames(N_TPM_cfRNA)<-rownames(TPM_cfRNA)

idx_X<-grep("X",colnames(TPM_cfRNA))
X_TPM_cfRNA<-TPM_cfRNA[,idx_X]
rownames(X_TPM_cfRNA)<-rownames(TPM_cfRNA)

colnames(T_TPM_tissue) <- sub(paste0("*", "T"), "", colnames(T_TPM_tissue))
colnames(P_TPM_tissue) <- sub(paste0("*", "P"), "", colnames(P_TPM_tissue))
colnames(N_TPM_cfRNA) <- sub(paste0("*", "_N"), "", colnames(N_TPM_cfRNA))
colnames(X_TPM_cfRNA) <- sub(paste0("*", "_X"), "", colnames(X_TPM_cfRNA))


all_paired_list<-Reduce(intersect, list(colnames(T_TPM_tissue),
                                        colnames(P_TPM_tissue),
                                        colnames(N_TPM_cfRNA),
                                        colnames(X_TPM_cfRNA)))

#all_paired_list<-as.data.frame(all_paired_list)
#write.table(all_paired_list,"all_3_paired_list.txt",col.names = F,row.names=F,quote = F)

###24 all paired samples
paired_T_TPM_tissue<-T_TPM_tissue[,which(colnames(T_TPM_tissue) %in% all_paired_list)]
paired_P_TPM_tissue<-P_TPM_tissue[,which(colnames(P_TPM_tissue) %in% all_paired_list)]
paired_X_TPM_cfRNA<-X_TPM_cfRNA[,which(colnames(X_TPM_cfRNA) %in% all_paired_list)]
paired_N_TPM_cfRNA<-N_TPM_cfRNA[,which(colnames(N_TPM_cfRNA) %in% all_paired_list)]


all.equal(colnames(paired_T_TPM_tissue),colnames(paired_P_TPM_tissue),
          colnames(paired_X_TPM_cfRNA),colnames(paired_N_TPM_cfRNA))

all.equal(rownames(paired_T_TPM_tissue),rownames(paired_P_TPM_tissue),
          rownames(paired_X_TPM_cfRNA),rownames(paired_N_TPM_cfRNA))

paired_T_TPM_tissue_rank_list<-as.data.frame(apply(-paired_T_TPM_tissue,2,rank))
paired_P_TPM_tissue_rank_list<-as.data.frame(apply(-paired_P_TPM_tissue,2,rank))
paired_X_TPM_cfRNA_rank_list<-as.data.frame(apply(-paired_X_TPM_cfRNA,2,rank))
paired_N_TPM_cfRNA_rank_list<-as.data.frame(apply(-paired_N_TPM_cfRNA,2,rank))

cor_list<-matrix(ncol = 6,nrow = 24)
rownames(cor_list)<-all_paired_list
colnames(cor_list)<-c("Tumor_NAT","Tumor_Plasma","Tumor_Urine","NAT_Plasma","NAT_Urine","Plasma_Urine")
for (k in 1:24){
  #top10000_N_idx <- order(paired_N_tissue_count_rank_list[,k], decreasing = F)[1:10000]
  #top10000_X_idx <- order(paired_X_tissue_count_rank_list[,k], decreasing = F)[1:10000]
  cor_list[k,1]<-cor(paired_T_TPM_tissue_rank_list[,k],paired_P_TPM_tissue_rank_list[,k],method = 'spearman')
  cor_list[k,2]<-cor(paired_T_TPM_tissue_rank_list[,k],paired_X_TPM_cfRNA_rank_list[,k],method = 'spearman')
  cor_list[k,3]<-cor(paired_T_TPM_tissue_rank_list[,k],paired_N_TPM_cfRNA_rank_list[,k],method = 'spearman')
  cor_list[k,4]<-cor(paired_P_TPM_tissue_rank_list[,k],paired_X_TPM_cfRNA_rank_list[,k],method = 'spearman')
  cor_list[k,5]<-cor(paired_P_TPM_tissue_rank_list[,k],paired_N_TPM_cfRNA_rank_list[,k],method = 'spearman')
  cor_list[k,6]<-cor(paired_X_TPM_cfRNA_rank_list[,k],paired_N_TPM_cfRNA_rank_list[,k],method = 'spearman')
  #cor_list[k,7]<-cor(paired_T_tissue_count_rank_list[top10000_N_idx,k],paired_N_tissue_count_rank_list[top10000_N_idx,k],method = 'spearman')
  #cor_list[k,8]<-cor(paired_P_tissue_count_rank_list[top10000_N_idx,k],paired_N_tissue_count_rank_list[top10000_N_idx,k],method = 'spearman')
}
cor_list<-as.data.frame(cor_list)

write.csv(cor_list,"cor_list_all.csv")
library(pheatmap)

#png("heatmap.png", width = 1600, height = 1800,res = 300)
###annotation

cor_list_plot<-t(cor_list)


coldata<-matrix(ncol = 0,nrow = 24)
rownames(coldata)<-colnames(cor_list_plot)

clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BC")
MIBC<-as.data.frame(cbind(clinicial_data$Patient,clinicial_data$MIBC,clinicial_data$AGE,clinicial_data$SEX,clinicial_data$Tumor_size))
colnames(MIBC)<-c("Patients","MIBC","Age","Sex","Tumor_size")

MIBC<-as.data.frame(MIBC)
coldata<-as.data.frame(coldata)
coldata$sample_num<-rownames(coldata)
coldata$Patients<-rownames(coldata)

coldata_all<-merge(coldata,MIBC,by="Patients",all.x=T)
rownames(coldata_all)<-coldata_all$sample_num
coldata_all<-coldata_all[,-2]
coldata_all<-coldata_all[,-1]

t_coldata_all<-as.data.frame(t(coldata_all))
t_coldata_all <- t_coldata_all %>% dplyr::select(colnames(cor_list_plot))

coldata_all<-as.data.frame(t(t_coldata_all))
coldata_all$Tumor_size<-as.numeric(coldata_all$Tumor_size)
coldata_all$Age<-as.numeric(coldata_all$Age)

ann_colors <- list(
  Sex = c("Male" = "#406A93", "Female" = "#EFA9AE"),
  MIBC = c("MIBC" = "#65AE00", "NMIBC" = "#D3B356"),
  Group = c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4"),
  Tumor_size=colorRampPalette(c("gray", "#2C1714"))(100),
  Age = colorRampPalette(c("#EEE0CC", "#91191E"))(100)
)

#CairoPDF(file = "heatmap_cor_4_samples.pdf", width = 10, height = 5,family='Arial')
CairoPNG(file = "heatmap_cor_4_samples.png", width = 10, height = 5,res=600,units='in',family='Arial')

pheatmap(cor_list_plot,
         cluster_rows = TRUE,
         show_rownames = T,
         cluster_cols = TRUE,
         show_colnames = F,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = coldata_all, 
         annotation_colors = ann_colors, 
         fontsize_row =8,
         fontsize = 10,
         na_col = "gray",
         scale = "column",
         #color = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
         color = colorRampPalette(c("blue", "white", "red"))(100),
         legend=TRUE,
         name = 'Z-score')
dev.off()
