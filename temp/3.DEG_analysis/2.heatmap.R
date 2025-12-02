library(edgeR)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(readxl)
library(dplyr)
library(tximport)
library(DESeq2)

sampleTable <- data.frame(row.names=additional_list)
idx_S<-grep('S',rownames(sampleTable))
sampleTable$type[idx_S]<-'BS'
idx_C<-grep('C',rownames(sampleTable))
sampleTable$type[idx_C]<-'BC'
idx_H<-grep('H',rownames(sampleTable))
sampleTable$type[idx_H]<-'HC'

dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~type)

#dds_counts <- estimateSizeFactors(dds); 
#dds_counts_N<-counts(dds_counts, normalized=TRUE)
#write.csv(dds_counts_N,"DeSeq_2_nor_dds_counts_N.csv")

smallestGroupSize <- 10
keep <- rowSums(counts(dds) > 0) >= smallestGroupSize
dds <- dds[keep,]
dds$type <- relevel(dds$type, ref = "HC")

dds <- DESeq(dds)
resultsNames(dds)
N_dds<-assay(dds) 
#write.csv(N_dds,"N_dds.csv")

ntd <- normTransform(dds)
library("pheatmap")
df <- as.data.frame(colData(dds)[,c("type")])
vsd <- vst(dds, blind=FALSE)


log2_cpm_matrix <- assay(ntd)

DGE_C_H<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_repeataware/salmon_repeataware_stranded_allgenome/N_type_C_vs_H_result.csv",header=T,row.names = 1)
DGE_C_S<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_repeataware/salmon_repeataware_stranded_allgenome/N_type_C_vs_S_result.csv",header=T,row.names = 1)
DGE_S_H<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_repeataware/salmon_repeataware_stranded_allgenome/N_type_S_vs_H_result.csv",header=T,row.names = 1)

logFCThreshold = 2
adjPThreshold = 0.01

select_DGE_C_H <- DGE_C_H %>% filter(padj < adjPThreshold, baseMean > 10, log2FoldChange > logFCThreshold)
select_DGE_C_S <- DGE_C_S %>% filter(padj < adjPThreshold, baseMean > 10, log2FoldChange > logFCThreshold)
select_DGE_S_H <- DGE_S_H %>% filter(padj < adjPThreshold, baseMean > 10, log2FoldChange > logFCThreshold)


all_genes <- rownames(log2_cpm_matrix)

Significance <- rep(NA, length(all_genes))
names(Significance) <- all_genes

Significance[all_genes %in% setdiff(rownames(select_DGE_C_H), rownames(select_DGE_C_S))] <- "N_C_H_only_DEG"
Significance[all_genes %in% setdiff(rownames(select_DGE_C_S), rownames(select_DGE_C_H))] <- "N_C_S_only_DEG"
Significance[all_genes %in% intersect(rownames(select_DGE_C_H), rownames(select_DGE_C_S))] <- "both_DEG"

keep_genes <- !is.na(Significance)
scaled_matrix <- t(scale(t(log2_cpm_matrix[keep_genes, ])))
Significance <- data.frame(Significance = Significance[keep_genes])

genelist<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome/all_genelist_annotation.csv",header=T,row.names = 1)
genelist<-genelist[,4:5]

genelist[grep("DNA",genelist$genetype),2]<-"Other"
genelist[grep("IG",genelist$genetype),2]<-"Other"
genelist[grep("TR",genelist$genetype),2]<-"Other"
genelist[grep("LTR",genelist$genetype),2]<-"LTR"
genelist[grep("rRNA",genelist$genetype),2]<-"Other"
genelist[grep("tRNA",genelist$genetype),2]<-"tRNA"
genelist[grep("decay",genelist$genetype),2]<-"Other"
genelist[grep("ribozyme",genelist$genetype),2]<-"Other"
genelist[grep("lncRNA",genelist$genetype),2]<-"lncRNA"
genelist[grep("protein_coding",genelist$genetype),2]<-"protein_coding"
genelist[grep("RC",genelist$genetype),2]<-"RC"
genelist[grep("pseudogene",genelist$genetype),2]<-"Other"
genelist[grep("SINE",genelist$genetype),2]<-"SINE"
genelist[grep("artifact",genelist$genetype),2]<-"Other"
genelist[grep("scaRNA",genelist$genetype),2]<-"Other"
genelist[grep("scRNA",genelist$genetype),2]<-"Other"
genelist[grep("snoRNA",genelist$genetype),2]<-"Other"
genelist[grep("snRNA",genelist$genetype),2]<-"Other"
genelist[grep("sRNA",genelist$genetype),2]<-"Other"
genelist[grep("srpRNA",genelist$genetype),2]<-"Other"
genelist[grep("Unknown",genelist$genetype),2]<-"Other"
genelist[grep("vault_RNA",genelist$genetype),2]<-"Other"
genelist[grep("TEC",genelist$genetype),2]<-"Other"
genelist[grep("retained_intron",genelist$genetype),2]<-"Other"
genelist[grep("processed_transcript",genelist$genetype),2]<-"Other"
genelist[grep("misc_RNA",genelist$genetype),2]<-"Other"

genelist_clean <- genelist %>%
  group_by(gene_symbol) %>%
  arrange(genetype == "Other") %>%  # 
  ungroup()

genelist_clean<-genelist_clean[!duplicated(genelist_clean$gene_symbol),]

Significance$gene_symbol<-rownames(Significance)
Significance<-merge(Significance,genelist_clean,by='gene_symbol',all.x=T)
rownames(Significance)<-Significance[,1]
Significance<-Significance[,-1]

library(readxl)
clinicial_data_BC <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BC")
clinicial_data_BS <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BS")
clinicial_data_HC <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "HC")

MIBC <- data.frame(Patient = clinicial_data_BC$Patient,
                   MIBC = clinicial_data_BC$MIBC,
                   Age = clinicial_data_BC$AGE,
                   Sex = clinicial_data_BC$SEX,
                   Tumor_size = clinicial_data_BC$Tumor_size)

BS_cd <- data.frame(Patient = clinicial_data_BS$Patient,
                    MIBC = NA,
                    Age = clinicial_data_BS$AGE,
                    Sex = clinicial_data_BS$SEX,
                    Tumor_size = NA)

HC_cd <- data.frame(Patient = clinicial_data_HC$Patient,
                    MIBC = NA,
                    Age = clinicial_data_HC$AGE,
                    Sex = clinicial_data_HC$SEX,
                    Tumor_size = NA)

MIBC_all <- rbind(MIBC, BS_cd, HC_cd)

coldata <- data.frame(sample_num = colnames(scaled_matrix))
coldata$Patients <- sub("_X|_N", "", coldata$sample_num)

coldata_all <- left_join(coldata, MIBC_all, by = c("Patients" = "Patient"))
rownames(coldata_all) <- coldata_all$sample_num

coldata_all$Group <- ifelse(grepl("^C", coldata_all$sample_num), "BC",
                           ifelse(grepl("^S", coldata_all$sample_num), "BS", "HC"))

coldata_all <- coldata_all[colnames(scaled_matrix), ]

coldata_all$Tumor_size <- as.numeric(coldata_all$Tumor_size)
coldata_all$Age <- as.numeric(coldata_all$Age)
coldata_all<-coldata_all[,-c(1:2)]

ann_colors <- list(
  Sex = c("Male" = "#406A93", "Female" = "#EFA9AE"),
  MIBC = c("MIBC" = "#65AE00", "NMIBC" = "#D3B356"),
  Group = c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4"),
  Tumor_size=colorRampPalette(c("gray", "#2C1714"))(100),
  Age = colorRampPalette(c("#EEE0CC", "#91191E"))(100)
)
row_anno_colors <- list(Significance = c(
  'N_C_H_only_DEG' = '#64BE52',
  'N_C_S_only_DEG' = '#D78203',
  'both_DEG' = '#373F89'),
genetype = c(
  'LINE' = 'darkred', 
  'SINE' = 'blue', 
  'protein_coding' = 'darkgreen', 
  'Simple_repeat' = 'pink', 
  'LTR' = 'orange', 
  'Low_complexity' = 'purple',
  'lncRNA'='lightblue',
  'Other'='grey',
  'Satellite'='yellow'
))
#install.packages("Cairo", type = "source")

library(Cairo)
#CairoPNG(file="heatmap.png",width=15,height=10,res=900,units = "in",family='Arial')
CairoPDF(file = "heatmap.pdf", width = 15, height = 10,family='Arial')

pheatmap(
  scaled_matrix,
  cluster_rows = TRUE,
  show_rownames = FALSE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = coldata_all,
  annotation_row = Significance,
  annotation_colors = c(ann_colors, row_anno_colors),
  na_col = "gray",
  scale = "none",
  #color = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  legend=TRUE,
  name = 'Z-score'
)

dev.off()
setdiff(unique(Significance$genetype), names(row_anno_colors$genetype))
