library(readxl)


count_cfRNA<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.standard/count_all_anno_by_gene_name.csv",header=T,row.names = 2)
count_cfRNA<-count_cfRNA[,-1]
TSP_data<-read.table("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/cell_type_marker/database/tsp_v1_basisMatrix.txt",header=T,sep='\t',row.names=1)

count_cfRNA_select<-count_cfRNA[which(rownames(count_cfRNA) %in%rownames(TSP_data)),]
TSP_data_select<-TSP_data[which(rownames(TSP_data) %in%rownames(count_cfRNA)),]
TSP_data_select <- TSP_data_select[match(rownames(count_cfRNA_select), rownames(TSP_data_select)), ]

all.equal(rownames(TSP_data_select),rownames(count_cfRNA_select))

basis_matrix_rank_list<-apply(-TSP_data_select,2,rank)
counts_rank_list<-apply(-count_cfRNA_select,2,rank)


cor_table<-matrix(nrow=62,ncol=222)
rownames(cor_table)<-colnames(basis_matrix_rank_list)
colnames(cor_table)<-colnames(counts_rank_list)

for (k in 1:222){
  for (i in 1:62){
    cor_table[i,k]<-cor(counts_rank_list[,k],basis_matrix_rank_list[,i],method = 'spearman')
    #print(i)
  }
  print(k)
}

write.csv(cor_table,"spearman_cor_between_cf_basismatrix.csv")

N_cor_table<-cor_table[,grep('N',colnames(cor_table))]

data_all<-cor_table

idx_R<-grep('C89',colnames(data_all))
data_all<-data_all[,-idx_R]
#idx_R<-grep('C103',colnames(data_all))
#data_all<-data_all[-idx_R,]
idx_R<-grep('C107',colnames(data_all))
data_all<-data_all[,-idx_R]

cor_table_plot<-as.data.frame(data_all)

coldata<-matrix(ncol = 2,nrow = 219)
rownames(coldata)<-colnames(data_all)

idx_X_counts<-grep("X",colnames(data_all))
idx_N_counts<-grep("N",colnames(data_all))
idx_S_counts<-grep("S",colnames(data_all))
idx_H_counts<-grep("H",colnames(data_all))
idx_C_counts<-grep("C",colnames(data_all))

coldata[idx_X_counts,1]<-"Plasma"
coldata[idx_N_counts,1]<-"Urine"
coldata[idx_H_counts,2]<-"HC"
coldata[idx_S_counts,2]<-"BS"
coldata[idx_C_counts,2]<-"BC"

library(readxl)
colnames(coldata)<-c("category","group")

clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BC")
MIBC<-as.data.frame(cbind(clinicial_data$Patient,clinicial_data$MIBC,clinicial_data$AGE,clinicial_data$SEX,clinicial_data$Tumor_size))
colnames(MIBC)<-c("individual","MIBC","Age","Sex","Tumor_size")

BS_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BS")
BS_cd<-as.data.frame(cbind(BS_clinicial_data$Patient,NA,BS_clinicial_data$AGE,BS_clinicial_data$SEX,NA))
colnames(BS_cd)<-c("individual","MIBC","Age","Sex","Tumor_size")

HC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "HC")
HC_cd<-as.data.frame(cbind(HC_clinicial_data$Patient,NA,HC_clinicial_data$AGE,HC_clinicial_data$SEX,NA))
colnames(HC_cd)<-c("individual","MIBC","Age","Sex","Tumor_size")

MIBC_all<-rbind(MIBC,HC_cd,BS_cd)

colnames(MIBC_all)<-c("Patients","MIBC","Age","Sex","Tumor_size")
MIBC_all<-as.data.frame(MIBC_all)
coldata<-as.data.frame(coldata)
coldata$sample_num<-rownames(coldata)
coldata$Patients<-rownames(coldata)
coldata$Patients <- sub(paste0("*", "_X"), "", coldata$Patients)
coldata$Patients <- sub(paste0("*", "_N"), "", coldata$Patients)


coldata_all<-merge(coldata,MIBC_all,by="Patients",all.x=T)
rownames(coldata_all)<-coldata_all$sample_num
coldata_all<-coldata_all[,-1]
coldata_all<-coldata_all[,-3]

t_coldata_all<-as.data.frame(t(coldata_all))
t_coldata_all <- t_coldata_all %>% dplyr::select(colnames(cor_table_plot))

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

cor_table_plot_2<-as.matrix(cor_table_plot)
library(pheatmap)

png("heatmap.png", width = 8000, height = 3600,res = 300)

pheatmap(cor_table_plot_2, cluster_rows=T, show_rownames=T,cluster_cols=T,show_colnames = F, 
         annotation_col=coldata_all,annotation_colors = ann_colors,na_col = "gray",scale = "column",fontsize_row = 5)
dev.off()

N_cor_table<-cor_table[,grep('N',colnames(cor_table))]

data_all<-N_cor_table

#idx_R<-grep('C89',colnames(data_all))
#data_all<-data_all[,-idx_R]
#idx_R<-grep('C103',colnames(data_all))
#data_all<-data_all[-idx_R,]
idx_R<-grep('C107',colnames(data_all))
data_all<-data_all[,-idx_R]

cor_table_plot<-as.data.frame(data_all)

coldata<-matrix(ncol = 2,nrow = 104)
rownames(coldata)<-colnames(data_all)

idx_X_counts<-grep("X",colnames(data_all))
idx_N_counts<-grep("N",colnames(data_all))
idx_S_counts<-grep("S",colnames(data_all))
idx_H_counts<-grep("H",colnames(data_all))
idx_C_counts<-grep("C",colnames(data_all))

coldata[idx_X_counts,1]<-"Plasma"
coldata[idx_N_counts,1]<-"Urine"
coldata[idx_H_counts,2]<-"HC"
coldata[idx_S_counts,2]<-"BS"
coldata[idx_C_counts,2]<-"BC"

library(readxl)
colnames(coldata)<-c("category","Group")

clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BC")
MIBC<-as.data.frame(cbind(clinicial_data$Patient,clinicial_data$MIBC,clinicial_data$AGE,clinicial_data$SEX,clinicial_data$Tumor_size))
colnames(MIBC)<-c("individual","MIBC","Age","Sex","Tumor_size")

BS_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BS")
BS_cd<-as.data.frame(cbind(BS_clinicial_data$Patient,NA,BS_clinicial_data$AGE,BS_clinicial_data$SEX,NA))
colnames(BS_cd)<-c("individual","MIBC","Age","Sex","Tumor_size")

HC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "HC")
HC_cd<-as.data.frame(cbind(HC_clinicial_data$Patient,NA,HC_clinicial_data$AGE,HC_clinicial_data$SEX,NA))
colnames(HC_cd)<-c("individual","MIBC","Age","Sex","Tumor_size")

MIBC_all<-rbind(MIBC,HC_cd,BS_cd)

colnames(MIBC_all)<-c("Patients","MIBC","Age","Sex","Tumor_size")
MIBC_all<-as.data.frame(MIBC_all)
coldata<-as.data.frame(coldata)
coldata$sample_num<-rownames(coldata)
coldata$Patients<-rownames(coldata)
coldata$Patients <- sub(paste0("*", "_X"), "", coldata$Patients)
coldata$Patients <- sub(paste0("*", "_N"), "", coldata$Patients)


coldata_all<-merge(coldata,MIBC_all,by="Patients",all.x=T)
rownames(coldata_all)<-coldata_all$sample_num
coldata_all<-coldata_all[,-1]
coldata_all<-coldata_all[,-3]

t_coldata_all<-as.data.frame(t(coldata_all))
t_coldata_all <- t_coldata_all %>% dplyr::select(colnames(cor_table_plot))

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
cor_table_plot_2<-as.matrix(cor_table_plot)
library(pheatmap)

png("heatmap_N.png", width = 8000, height = 3600,res = 300)

pheatmap(cor_table_plot_2, cluster_rows=T, show_rownames=T,cluster_cols=T,show_colnames = F, 
         annotation_col=coldata_all,annotation_colors = ann_colors,na_col = "gray",scale = "column",fontsize_row = 6)
dev.off()


# Separate samples by group
idx_C <- which(coldata_all$Group == "BC")
idx_H <- which(coldata_all$Group == "HC")
idx_S <- which(coldata_all$Group == "BS")

# Cluster each group separately
hc_C <- hclust(dist(t(cor_table[, idx_C])))
hc_H <- hclust(dist(t(cor_table[, idx_H])))
hc_S <- hclust(dist(t(cor_table[, idx_S])))

# Combine the clustering results by group
ordered_indices <- c(idx_C[hc_C$order], idx_H[hc_H$order], idx_S[hc_S$order])

# Reorder cor_table and coldata_all based on combined clustering order
cor_table_ordered <- cor_table_plot_2[, ordered_indices]
coldata_all_ordered <- coldata_all[ordered_indices, ]

#CairoPDF(file = "heatmap_N_group.pdf", width = 20, height = 10,family='Arial')
CairoPNG(file = "heatmap_N_group.png", width = 20, height = 10,res=600,units='in',family='Arial')

split_name <- function(name, n = 90) {
  # 每 n 个字符插入一个换行符
  gsub(paste0("(.{", n, "})"), "\\1\n", name)
}
cor_table_ordered_new <- cor_table_ordered
rownames(cor_table_ordered_new) <- sapply(rownames(cor_table_ordered), split_name)


# Plot the heatmap with the ordered data
pheatmap(cor_table_ordered_new, 
         cluster_rows = TRUE, 
         show_rownames = TRUE, 
         cluster_cols = FALSE,  # Disable automatic clustering since we pre-ordered columns
         show_colnames = FALSE, 
         annotation_col = coldata_all_ordered, 
         annotation_colors = ann_colors, 
         na_col = "gray", 
         scale = "column",
         fontsize_row = 4.4,
         fontsize = 10,
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         #color = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
         color = colorRampPalette(c("blue", "white", "red"))(100),
         legend=TRUE,
         name = 'Z-score')
dev.off()
