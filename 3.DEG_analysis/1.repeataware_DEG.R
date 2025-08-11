library(Biostrings)
library(BSgenome)
library(eisaR)
library(GenomicFeatures)
library(SummarizedExperiment)
library(tximeta)
library(tximport)
library(rjson)
library(reticulate)
library(SingleCellExperiment)
library(scater)
library(BiocFileCache)
library(tximportData)
library(DESeq2)
library(ggrepel)

####TRY USE TXIMPORT


#setwd("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING")

#test tximeta on process
name<-read.table("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/cfRNA.txt",header = F,sep = "\t",fill=T)
name<-as.vector(name$V1)

files <- file.path("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome", name, "quant.sf") 
file.exists(files)
additional_list<-name[which(file.exists(files)==TRUE)]
#repeataware_all_count <- read.csv("salmon.repeataware/repeataware_all_count.csv",row.names = 1)

QC_list<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/QC_repeataware_list.csv",row.names = 1)
additional_list<-additional_list[which(additional_list %in% QC_list[,1])]

idx_R<-grep('C89',additional_list)
additional_list<-additional_list[-idx_R]
#idx_R<-grep('C103',colnames(count_analysis))
#count_analysis<-count_analysis[,-idx_R]
idx_R<-grep('C107',additional_list)
additional_list_all<-additional_list[-idx_R]


########### N or X 
idx<-grep('N',additional_list_all)
additional_list<-additional_list_all[idx]

files <- file.path("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome", additional_list, "quant.sf") 

#coldata <- data.frame(files, names=files, condition="A", stringsAsFactors=FALSE)
#coldata
tx2gene<-read.csv("/dssg/home/acct-dahan/share/references/gencode/RE_allgenome/gencode.v44.ucsc.rmsk.tx.to.gene.csv",header=F)


library(tximport)

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
#save.image(file='/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_repeataware/salmon_repeataware_stranded_allgenome/N_data.RData')

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
#rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$type
colnames(sampleDistMatrix) <- rownames(sampleTable)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

dds <- estimateSizeFactors(dds)

nor_count_N <- assay(vst(dds,blind=F))
#write.csv(nor_count_N,"nor_count_N.csv")

plotPCA(vsd, intgroup=c("type"))
# Assuming 'vsd' is your variance-stabilizing transformation object
pca_plot <- plotPCA(vsd, intgroup = c("type"), returnData = TRUE)

# Create the plot and add sample names as labels
ggplot(pca_plot, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 3) +
  #geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
  scale_color_manual(values = c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4"))+
  stat_ellipse(aes(group = type), linetype = "dashed", size = 1) +  
  theme_minimal()+
  theme_minimal() +
  labs(x = "PC1 (64%)", y = "PC2 (7%)", color = "Group") +  
  theme(
    legend.title = element_text(size = 10),  
    legend.text = element_text(size = 10)    
  )
#ggsave("PCA_N_all.png",height=4,width=5,dpi=600)
#ggsave("PCA_N_all.pdf",width=5,height=4,dpi=1200, units = "in", device = cairo_pdf)


res_type_C_vs_H <- results(dds, name="type_BC_vs_HC")
#write.csv(res_type_C_vs_H,"N_type_C_vs_H_result.csv")
hist(res_type_C_vs_H$pvalue, col = "grey50", border = "white", main = "", xlab = "p-values")

type_C_vs_H_result<- lfcShrink(dds, coef ="type_BC_vs_HC",type = 'apeglm',res=res_type_C_vs_H)
#write.csv(type_C_vs_H_result,"lfc_N_type_C_vs_H_result.csv")


res_type_S_vs_H <- results(dds, name="type_BS_vs_HC")
#write.csv(res_type_S_vs_H,"N_type_S_vs_H_result.csv")
hist(res_type_S_vs_H$pvalue, col = "grey50", border = "white", main = "", xlab = "p-values")

type_S_vs_H_result<- lfcShrink(dds, coef ="type_BS_vs_HC",type = 'apeglm',res=res_type_S_vs_H)
#write.csv(type_S_vs_H_result,"lfc_N_type_S_vs_H_result.csv")

dds$type <- relevel(dds$type, ref = "BS")

dds <- DESeq(dds)
resultsNames(dds)

res_type_C_vs_S <- results(dds, name="type_BC_vs_BS")
#write.csv(res_type_C_vs_S,"N_type_C_vs_S_result.csv")
hist(res_type_C_vs_S$pvalue, col = "grey50", border = "white", main = "", xlab = "p-values")
type_C_vs_S_result<- lfcShrink(dds, coef ="type_BC_vs_BS",type = 'apeglm',res=res_type_C_vs_S)
#write.csv(type_C_vs_S_result,"lfc_N_type_C_vs_S_result.csv")

dataForVolcanoPlot <- as.data.frame(type_C_vs_H_result)

logFcThreshold <- log2(2)
adjPvalThreshold <- 0.05
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange < logFcThreshold | padj> adjPvalThreshold)] <- 'NS'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange >= logFcThreshold & padj <= adjPvalThreshold)] <- 'UP'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange <= -logFcThreshold & padj <= adjPvalThreshold)] <- 'DOWN'


dataForVolcanoPlot
dataForVolcanoPlot$ID = rownames(dataForVolcanoPlot)

cols <- c('UP'="red", 'NS'='gray','DOWN'="blue")

dev.new()
p <- ggplot(dataForVolcanoPlot, aes(x = log2FoldChange, y = -log10(padj))) +
  #xlim(-2,2) +
  labs(x=expression(bold('Log'['2']*'(Fold Change)')), 
       y=(expression(bold('-Log'['10']*'(FDR)'))), 
       title=NULL) +
  geom_point(color='gray', alpha=1, size=2) +
  geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
             color='darkgreen', linetype='dashed') +
  geom_hline(yintercept = -log10(adjPvalThreshold), 
             color='darkgreen',linetype='dashed')+
  scale_color_manual(values = cols) +
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'UP'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="red", 
             size=1, color="red",shape=21)+
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'DOWN'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="blue", 
             size=1, color="blue",shape=21) +
  scale_fill_manual(values = cols) +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'UP' & (-log10(padj) > 2 & log2FoldChange > 5)),
                  aes(label = ID), nudge_x = 0.1,
                  segment.alpha = 0.3, #segment.size = 0.5,
                  min.segment.length = 0.2,max.overlaps = 15,
                  size = 3, color='black', segment.color = 'black') +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'DOWN' & log2FoldChange < -5 & -log10(padj) > 2),
                  aes(label = ID),
                  segment.alpha = 0.3, #segment.size = 0.5,
                  min.segment.length = 0.2, max.overlaps = 15,
                  size = 3, color='black', segment.color = 'black') +
  theme_minimal()+theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),     
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    #legend.background = element_rect(fill = "white", color = "black", size = 0.5), 
    legend.key = element_rect(color = "black", size = 0.5) 
  )
p

ggsave("N_type_C_vs_H_volcano.png", width = 8, height = 8, units = "in", dpi = 300)
ggsave("N_type_C_vs_H_volcano.pdf", width = 8, height = 8, units = "in", dpi = 300)


dataForVolcanoPlot <- as.data.frame(type_C_vs_S_result)

logFcThreshold <- log2(2)
adjPvalThreshold <- 0.05
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange < logFcThreshold | padj> adjPvalThreshold)] <- 'NS'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange >= logFcThreshold & padj <= adjPvalThreshold)] <- 'UP'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange <= -logFcThreshold & padj <= adjPvalThreshold)] <- 'DOWN'

dataForVolcanoPlot
dataForVolcanoPlot$ID = rownames(dataForVolcanoPlot)

cols <- c('UP'="red", 'NS'='gray','DOWN'="blue")

dev.new()
p <- ggplot(dataForVolcanoPlot, aes(x = log2FoldChange, y = -log10(padj))) +
  #xlim(-2,2) +
  labs(x=expression(bold('Log'['2']*'(Fold Change)')), 
       y=(expression(bold('-Log'['10']*'(FDR)'))), 
       title=NULL) +
  geom_point(color='gray', alpha=1, size=2) +
  geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
             color='darkgreen', linetype='dashed') +
  geom_hline(yintercept = -log10(adjPvalThreshold), 
             color='darkgreen',linetype='dashed')+
  scale_color_manual(values = cols) +
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'UP'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="red", 
             size=1, color="red",shape=21)+
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'DOWN'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="blue", 
             size=1, color="blue",shape=21) +
  scale_fill_manual(values = cols) +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'UP' & (-log10(padj) > 2 & log2FoldChange > 5)),
                  aes(label = ID), nudge_x = 0.1,
                  segment.alpha = 0.3, #segment.size = 0.5,
                  min.segment.length = 0.2,max.overlaps = 15,
                  size = 3, color='black', segment.color = 'black') +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'DOWN' & log2FoldChange < -5 & -log10(padj) > 2),
                  aes(label = ID),
                  segment.alpha = 0.3, #segment.size = 0.5,
                  min.segment.length = 0.2, max.overlaps = 15,
                  size = 3, color='black', segment.color = 'black') +
  theme_minimal()+theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),   
    legend.key.spacing.y = unit(0.1, "cm"),   
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),     
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    #legend.background = element_rect(fill = "white", color = "black", size = 0.5), 
    legend.key = element_rect(color = "black", size = 0.5)  
  )
p

ggsave("N_type_C_vs_S_volcano.png", width = 8, height = 8, units = "in", dpi = 300)
ggsave("N_type_C_vs_S_volcano.pdf", width = 8, height = 8, units = "in", dpi = 300)

dataForVolcanoPlot <- as.data.frame(type_S_vs_H_result)

logFcThreshold <- log2(2)
adjPvalThreshold <- 0.05
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange < logFcThreshold | padj> adjPvalThreshold)] <- 'NS'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange >= logFcThreshold & padj <= adjPvalThreshold)] <- 'UP'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange <= -logFcThreshold & padj <= adjPvalThreshold)] <- 'DOWN'


dataForVolcanoPlot
dataForVolcanoPlot$ID = rownames(dataForVolcanoPlot)

cols <- c('UP'="red", 'NS'='gray','DOWN'="blue")

dev.new()
p <- ggplot(dataForVolcanoPlot, aes(x = log2FoldChange, y = -log10(padj))) +
  #xlim(-2,2) +
  labs(x=expression(bold('Log'['2']*'(Fold Change)')), 
       y=(expression(bold('-Log'['10']*'(FDR)'))), 
       title=NULL) +
  geom_point(color='gray', alpha=1, size=2) +
  geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
             color='darkgreen', linetype='dashed') +
  geom_hline(yintercept = -log10(adjPvalThreshold), 
             color='darkgreen',linetype='dashed')+
  scale_color_manual(values = cols) +
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'UP'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="red", 
             size=1, color="red",shape=21)+
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'DOWN'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="blue", 
             size=1, color="blue",shape=21) +
  scale_fill_manual(values = cols) +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'UP' & (-log10(padj) > 2 & log2FoldChange > 5)),
                  aes(label = ID), nudge_x = 0.1,
                  segment.alpha = 0.3, #segment.size = 0.5,
                  min.segment.length = 0.2,max.overlaps = 15,
                  size = 3, color='black', segment.color = 'black') +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'DOWN' & log2FoldChange < -5 & -log10(padj) > 2),
                  aes(label = ID),
                  segment.alpha = 0.3, #segment.size = 0.5,
                  min.segment.length = 0.2, max.overlaps = 15,
                  size = 3, color='black', segment.color = 'black') +
  theme_minimal()+theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),  
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),     
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    #legend.background = element_rect(fill = "white", color = "black", size = 0.5), # 图例背景和边框
    legend.key = element_rect(color = "black", size = 0.5)  
  )
p

ggsave("N_type_S_vs_H_volcano.png", width = 8, height = 8, units = "in", dpi = 300)
ggsave("N_type_S_vs_H_volcano.pdf", width = 8, height = 8, units = "in", dpi = 300)


library(clusterProfiler)
DEGresult<-as.data.frame(res_type_C_vs_H)
result_up_gene<-subset(DEGresult,padj<0.05&log2FoldChange>1)
result_down_gene<-subset(DEGresult,padj<0.05& (-log2FoldChange)>1)

#up
goAll_up <- enrichGO(rownames(result_up_gene),
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "ALL",
                     pAdjustMethod= 'BH',
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.2)   

p2<-dotplot(goAll_up, font.size = 8,
            #title= "GO analysis for Urine BC vs HC",
            label_format= 60,
            showCategory=20)+
  ggtitle("GO analysis for upregulated genes in Urine BC vs HC") +
  theme(plot.title = element_text(size = 8, hjust = 0.5))+
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),    
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),) +
  theme(legend.position = "right")

p2
ggsave("N_type_C_vs_H_up_KEGG.pdf", width = 8, height = 8, units = "in", dpi = 300)
ggsave("N_type_C_vs_H_up_KEGG.png", width = 8, height = 8, units = "in", dpi = 300)

#down
goAll_down <- enrichGO(rownames(result_down_gene),
                       keyType = "SYMBOL",
                       OrgDb = "org.Hs.eg.db",
                       ont = "ALL",
                       pAdjustMethod= 'BH',
                       pvalueCutoff= 0.05,
                       qvalueCutoff= 0.2)   

p2<-dotplot(goAll_down, font.size = 8,
            title= "GO Dot plot",
            label_format= 60,
            showCategory=20)+
  ggtitle("GO analysis for downregulated genes in Urine BC vs HC") +
  theme(plot.title = element_text(size = 8, hjust = 0.5))+
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),    
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),) +
  theme(legend.position = "right")
p2
ggsave("N_type_C_vs_H_down_KEGG.png", width = 8, height = 8, units = "in", dpi = 300)
ggsave("N_type_C_vs_H_down_KEGG.pdf", width = 8, height = 8, units = "in", dpi = 300)



DEGresult<-as.data.frame(res_type_C_vs_S)
result_up_gene<-subset(DEGresult,padj<0.05&log2FoldChange>1)
result_down_gene<-subset(DEGresult,padj<0.05& (-log2FoldChange)>1)

#up
goAll_up <- enrichGO(rownames(result_up_gene),
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "ALL",
                     pAdjustMethod= 'BH',
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.2)   

p2<-dotplot(goAll_up, font.size = 8,
            label_format= 60,
            showCategory=20)+
  ggtitle("GO analysis for upregulated genes in Urine BC vs BS") +
  theme(plot.title = element_text(size = 8, hjust = 0.5))+
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),    
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),) +
  theme(legend.position = "right")
p2
ggsave("N_type_C_vs_S_up_KEGG.png", width = 8, height = 8, units = "in", dpi = 300)
ggsave("N_type_C_vs_S_up_KEGG.pdf", width = 8, height = 8, units = "in", dpi = 300)



#down
goAll_down <- enrichGO(rownames(result_down_gene),
                       keyType = "SYMBOL",
                       OrgDb = "org.Hs.eg.db",
                       ont = "ALL",
                       pAdjustMethod= 'BH',
                       pvalueCutoff= 0.05,
                       qvalueCutoff= 0.2)   

p2<-dotplot(goAll_down, font.size = 8,
            label_format= 60,
            showCategory=20)+
  ggtitle("GO analysis for downregulated genes in Urine BC vs BS") +
  theme(plot.title = element_text(size = 8, hjust = 0.5))+
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),    
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),) +
  theme(legend.position = "right")
p2
#ggsave("N_type_C_vs_S_down_KEGG.png", width = 8, height = 8, units = "in", dpi = 300)
#ggsave("N_type_C_vs_S_down_KEGG.pdf", width = 8, height = 8, units = "in", dpi = 300)


DEGresult<-as.data.frame(res_type_S_vs_H)
result_up_gene<-subset(DEGresult,padj<0.05&log2FoldChange>1)
result_down_gene<-subset(DEGresult,padj<0.05& (-log2FoldChange)>1)

#up
goAll_up <- enrichGO(rownames(result_up_gene),
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "ALL",
                     pAdjustMethod= 'BH',
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.2)   
p2<-dotplot(goAll_up, font.size = 8,
            label_format= 60,
            showCategory=20)+
  ggtitle("GO analysis for upregulated genes in Urine BS vs HC") +
  theme(plot.title = element_text(size = 8, hjust = 0.5))+
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),    
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),) +
  theme(legend.position = "right")
p2
ggsave("N_type_S_vs_H_up_KEGG.png", width = 8, height = 8, units = "in", dpi = 300)
ggsave("N_type_S_vs_H_up_KEGG.pdf", width = 8, height = 8, units = "in", dpi = 300)


#down
goAll_down <- enrichGO(rownames(result_down_gene),
                       keyType = "SYMBOL",
                       OrgDb = "org.Hs.eg.db",
                       ont = "ALL",
                       pAdjustMethod= 'BH',
                       pvalueCutoff= 0.05,
                       qvalueCutoff= 0.2)   

p2<-dotplot(goAll_down, font.size = 8,
            label_format= 60,
            showCategory=20)+
  ggtitle("GO analysis for downregulated genes in Urine BS vs HC") +
  theme(plot.title = element_text(size = 8, hjust = 0.5))+
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),    
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),) +
  theme(legend.position = "right")
p2
ggsave("N_type_S_vs_H_down_KEGG.png", width = 8, height = 8, units = "in", dpi = 300)
ggsave("N_type_S_vs_H_down_KEGG.pdf", width = 8, height = 8, units = "in", dpi = 300)




########### N or X 
idx<-grep('X',additional_list_all)
additional_list<-additional_list_all[idx]
files <- file.path("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome", additional_list, "quant.sf") 

#coldata <- data.frame(files, names=files, condition="A", stringsAsFactors=FALSE)
#coldata
tx2gene<-read.csv("/dssg/home/acct-dahan/share/references/gencode/RE_allgenome/gencode.v44.ucsc.rmsk.tx.to.gene.csv",header=F)


library(tximport)

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
#save.image(file='/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_repeataware/salmon_repeataware_stranded_allgenome/X_data.RData')

sampleTable <- data.frame(row.names=additional_list)
idx_S<-grep('S',rownames(sampleTable))
sampleTable$type[idx_S]<-'BS'
idx_C<-grep('C',rownames(sampleTable))
sampleTable$type[idx_C]<-'BC'
idx_H<-grep('H',rownames(sampleTable))
sampleTable$type[idx_H]<-'HC'

dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~type)

#dds_counts <- estimateSizeFactors(dds); 
#dds_counts_X<-counts(dds_counts, normalized=TRUE)
#write.csv(dds_counts_X,"DeSeq_2_nor_dds_counts_X.csv")

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
#rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$type
colnames(sampleDistMatrix) <- rownames(sampleTable)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

nor_count_X <- assay(vst(dds,blind=F))
#write.csv(nor_count_X,"nor_count_X.csv")

plotPCA(vsd, intgroup=c("type"))
# Assuming 'vsd' is your variance-stabilizing transformation object
pca_plot <- plotPCA(vsd, intgroup = c("type"), returnData = TRUE)
# Create the plot and add sample names as labels
# Create the plot and add sample names as labels
ggplot(pca_plot, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 3) +
  #geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
  scale_color_manual(values = c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4"))+
  stat_ellipse(aes(group = type), linetype = "dashed", size = 1) +  
  theme_minimal()+
  theme_minimal() +
  labs(x = "PC1 (35%)", y = "PC2 (17%)", color = "Group") +  
  theme(
    legend.title = element_text(size = 10),  
    legend.text = element_text(size = 10)    
  )

#ggsave("PCA_X_all.png",height=4,width=5,dpi=600)
#ggsave("PCA_X_all.pdf",width=5,height=4,dpi=1200, units = "in", device = cairo_pdf)

res_type_C_vs_H <- results(dds, name="type_BC_vs_HC")
#write.csv(res_type_C_vs_H,"X_type_C_vs_H_result.csv")
hist(res_type_C_vs_H$pvalue, col = "grey50", border = "white", main = "", xlab = "p-values")
type_C_vs_H_result<- lfcShrink(dds, coef ="type_BC_vs_HC",type = 'apeglm',res=res_type_C_vs_H)
#write.csv(type_C_vs_H_result,"lfc_X_type_C_vs_H_result.csv")


res_type_S_vs_H <- results(dds, name="type_BS_vs_HC")
#write.csv(res_type_S_vs_H,"X_type_S_vs_C_result.csv")
hist(res_type_S_vs_H$pvalue, col = "grey50", border = "white", main = "", xlab = "p-values")
type_S_vs_H_result<- lfcShrink(dds, coef ="type_BS_vs_HC",type = 'apeglm',res=res_type_S_vs_H)
#write.csv(type_S_vs_H_result,"lfc_X_type_S_vs_H_result.csv")

dds$type <- relevel(dds$type, ref = "BS")

dds <- DESeq(dds)
resultsNames(dds)

res_type_C_vs_S <- results(dds, name="type_BC_vs_BS")
#write.csv(res_type_C_vs_S,"X_type_C_vs_S_result.csv")
hist(res_type_C_vs_S$pvalue, col = "grey50", border = "white", main = "", xlab = "p-values")
type_C_vs_S_result<- lfcShrink(dds, coef ="type_BC_vs_BS",type = 'apeglm',res=res_type_C_vs_S)
#write.csv(type_C_vs_S_result,"lfc_X_type_C_vs_S_result.csv")


dataForVolcanoPlot <- as.data.frame(type_C_vs_H_result)

logFcThreshold <- log2(2)
adjPvalThreshold <- 0.05
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange < logFcThreshold | padj> adjPvalThreshold)] <- 'NS'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange >= logFcThreshold & padj <= adjPvalThreshold)] <- 'UP'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange <= -logFcThreshold & padj <= adjPvalThreshold)] <- 'DOWN'


dataForVolcanoPlot
dataForVolcanoPlot$ID = rownames(dataForVolcanoPlot)

cols <- c('UP'="red", 'NS'='gray','DOWN'="blue")


dev.new()

p <- ggplot(dataForVolcanoPlot, aes(x = log2FoldChange, y = -log10(padj))) +
  #xlim(-2,2) +
  labs(x=expression(bold('Log'['2']*'(Fold Change)')), 
       y=(expression(bold('-Log'['10']*'(FDR)'))), 
       title=NULL) +
  geom_point(color='gray', alpha=1, size=2) +
  geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
             color='darkgreen', linetype='dashed') +
  geom_hline(yintercept = -log10(adjPvalThreshold), 
             color='darkgreen',linetype='dashed')+
  scale_color_manual(values = cols) +
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'UP'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="red", 
             size=1, color="red",shape=21)+
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'DOWN'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="blue", 
             size=1, color="blue",shape=21) +
  scale_fill_manual(values = cols) +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'UP' & (-log10(padj) > 2 & log2FoldChange > 5)),
                  aes(label = ID), nudge_x = 0.1,
                  segment.alpha = 0.5, #segment.size = 0.5,
                  min.segment.length = 0.2,max.overlaps = 20,
                  size = 3, color='black', segment.color = 'black') +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'DOWN' & log2FoldChange < -5 & -log10(padj) > 2),
                  aes(label = ID),
                  segment.alpha = 0.3, #segment.size = 0.5,
                  min.segment.length = 0.2, max.overlaps = 20,
                  size = 3, color='black', segment.color = 'black') +
  theme_minimal()+theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),     
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    #legend.background = element_rect(fill = "white", color = "black", size = 0.5), 
    legend.key = element_rect(color = "black", size = 0.5)  
  )
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
p

ggsave("X_type_C_vs_H_volcano.png", width = 8, height = 8, units = "in", dpi = 300)
ggsave("X_type_C_vs_H_volcano.pdf", width = 8, height = 8, units = "in", dpi = 300)


dataForVolcanoPlot <- as.data.frame(type_C_vs_S_result)

logFcThreshold <- log2(2)
adjPvalThreshold <- 0.05
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange < logFcThreshold | padj> adjPvalThreshold)] <- 'NS'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange >= logFcThreshold & padj <= adjPvalThreshold)] <- 'UP'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange <= -logFcThreshold & padj <= adjPvalThreshold)] <- 'DOWN'


dataForVolcanoPlot
dataForVolcanoPlot$ID = rownames(dataForVolcanoPlot)

#cols <- c('UP'=="red", 'NS'='gray','DOWN'="green")
cols <- c('UP'="red", 'NS'='gray','DOWN'="blue")


dev.new()

p <- ggplot(dataForVolcanoPlot, aes(x = log2FoldChange, y = -log10(padj))) +
  #xlim(-2,2) +
  labs(x=expression(bold('Log'['2']*'(Fold Change)')), 
       y=(expression(bold('-Log'['10']*'(FDR)'))), 
       title=NULL) +
  geom_point(color='gray', alpha=1, size=2) +
  geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
             color='darkgreen', linetype='dashed') +
  geom_hline(yintercept = -log10(adjPvalThreshold), 
             color='darkgreen',linetype='dashed')+
  scale_color_manual(values = cols) +
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'UP'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="red", 
             size=1, color="red",shape=21)+
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'DOWN'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="blue", 
             size=1, color="blue",shape=21) +
  scale_fill_manual(values = cols) +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'UP' & (-log10(padj) > 2 & log2FoldChange > 4)),
                  aes(label = ID), nudge_x = 0.1,
                  segment.alpha = 0.5, #segment.size = 0.5,
                  min.segment.length = 0.2,max.overlaps = 20,
                  size = 3, color='black', segment.color = 'black') +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'DOWN' & log2FoldChange < -2 & -log10(padj) > 2),
                  aes(label = ID),
                  segment.alpha = 0.3, #segment.size = 0.5,
                  min.segment.length = 0.2, max.overlaps = 20,
                  size = 3, color='black', segment.color = 'black') +
  theme_minimal()+theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),  
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),     
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    #legend.background = element_rect(fill = "white", color = "black", size = 0.5), 
    legend.key = element_rect(color = "black", size = 0.5)  
  )
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
p

ggsave("X_type_C_vs_S_volcano.png", width = 8, height = 8, units = "in", dpi = 300)
ggsave("X_type_C_vs_S_volcano.pdf", width = 8, height = 8, units = "in", dpi = 300)



dataForVolcanoPlot <- as.data.frame(type_S_vs_H_result)

logFcThreshold <- log2(2)
adjPvalThreshold <- 0.05
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange < logFcThreshold | padj> adjPvalThreshold)] <- 'NS'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange >= logFcThreshold & padj <= adjPvalThreshold)] <- 'UP'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange <= -logFcThreshold & padj <= adjPvalThreshold)] <- 'DOWN'


dataForVolcanoPlot
dataForVolcanoPlot$ID = rownames(dataForVolcanoPlot)

#cols <- c('UP'=="red", 'NS'='gray','DOWN'="green")
cols <- c('UP'="red", 'NS'='gray','DOWN'="blue")


dev.new()

p <- ggplot(dataForVolcanoPlot, aes(x = log2FoldChange, y = -log10(padj))) +
  #xlim(-2,2) +
  labs(x=expression(bold('Log'['2']*'(Fold Change)')), 
       y=(expression(bold('-Log'['10']*'(FDR)'))), 
       title=NULL) +
  geom_point(color='gray', alpha=1, size=2) +
  geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
             color='darkgreen', linetype='dashed') +
  geom_hline(yintercept = -log10(adjPvalThreshold), 
             color='darkgreen',linetype='dashed')+
  scale_color_manual(values = cols) +
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'UP'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="red", 
             size=1, color="red",shape=21)+
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'DOWN'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="blue", 
             size=1, color="blue",shape=21) +
  scale_fill_manual(values = cols) +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'UP' & (-log10(padj) > 2 & log2FoldChange > 4)),
                  aes(label = ID), nudge_x = 0.1,
                  segment.alpha = 0.5, #segment.size = 0.5,
                  min.segment.length = 0.2,max.overlaps = 20,
                  size = 3, color='black', segment.color = 'black') +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'DOWN' & log2FoldChange < -4 & -log10(padj) > 2),
                  aes(label = ID),
                  segment.alpha = 0.3, #segment.size = 0.5,
                  min.segment.length = 0.2, max.overlaps = 20,
                  size = 3, color='black', segment.color = 'black') +
  theme_minimal()+theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),  
    legend.spacing.y = unit(1, 'cm'),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),     
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
   #legend.background = element_rect(fill = "white", color = "black", size = 0.5), 
    legend.key = element_rect(color = "black", size = 0.5)  
  )
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))
p

ggsave("X_type_S_vs_H_volcano.png", width = 8, height = 8, units = "in", dpi = 300)
ggsave("X_type_S_vs_H_volcano.pdf", width = 8, height = 8, units = "in", dpi = 300)


library(clusterProfiler)
DEGresult<-as.data.frame(res_type_C_vs_H)
result_up_gene<-subset(DEGresult,padj<0.05&log2FoldChange>1)
result_down_gene<-subset(DEGresult,padj<0.05& (-log2FoldChange)>1)

#up
goAll_up <- enrichGO(rownames(result_up_gene),
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "ALL",
                     pAdjustMethod= 'BH',
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.2)   

p2<-dotplot(goAll_up, font.size = 8,
            label_format= 60,
            showCategory=20)+
  ggtitle("GO analysis for upregulated genes in Plasma BC vs HC") +
  theme(plot.title = element_text(size = 8, hjust = 0.5))+
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),    
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),) +
  theme(legend.position = "right")
p2
ggsave("X_type_C_vs_H_up_KEGG.png", width = 8, height = 8, units = "in", dpi = 300)
ggsave("X_type_C_vs_H_up_KEGG.pdf", width = 8, height = 8, units = "in", dpi = 300)


#down
goAll_down <- enrichGO(rownames(result_down_gene),
                       keyType = "SYMBOL",
                       OrgDb = "org.Hs.eg.db",
                       ont = "ALL",
                       pAdjustMethod= 'BH',
                       pvalueCutoff= 0.05,
                       qvalueCutoff= 0.2)   

p2<-dotplot(goAll_down, font.size = 8,
            label_format= 60,
            showCategory=20)+
  ggtitle("GO analysis for downregulated genes in Plasma BC vs HC") +
  theme(plot.title = element_text(size = 8, hjust = 0.5)) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),    
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),) +
  theme(legend.position = "right")
p2
ggsave("X_type_C_vs_H_down_KEGG.png", width = 8, height = 8, units = "in", dpi = 300)
ggsave("X_type_C_vs_H_down_KEGG.pdf", width = 8, height = 8, units = "in", dpi = 300)



DEGresult<-as.data.frame(res_type_C_vs_S)
result_up_gene<-subset(DEGresult,padj<0.05&log2FoldChange>1)
result_down_gene<-subset(DEGresult,padj<0.05& (-log2FoldChange)>1)

#up
goAll_up <- enrichGO(rownames(result_up_gene),
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "ALL",
                     pAdjustMethod= 'BH',
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.2)   

p2<-dotplot(goAll_up, font.size = 8,
            label_format= 60,
            showCategory=20)+
  ggtitle("GO analysis for upregulated genes in Plasma BC vs BS") +
  theme(plot.title = element_text(size = 8, hjust = 0.5))+
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),    
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),) +
  theme(legend.position = "right")
p2
ggsave("X_type_C_vs_S_up_KEGG.png", width = 8, height = 8, units = "in", dpi = 300)
ggsave("X_type_C_vs_S_up_KEGG.pdf", width = 8, height = 8, units = "in", dpi = 300)

#down
goAll_down <- enrichGO(rownames(result_down_gene),
                       keyType = "SYMBOL",
                       OrgDb = "org.Hs.eg.db",
                       ont = "ALL",
                       pAdjustMethod= 'BH',
                       pvalueCutoff= 0.05,
                       qvalueCutoff= 0.2)   

p2<-dotplot(goAll_down, font.size = 8,
            label_format= 60,
            showCategory=20)+
  ggtitle("GO analysis for downregulated genes in Plasma BC vs BS") +
  theme(plot.title = element_text(size = 8, hjust = 0.5))+
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),    
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),) +
  theme(legend.position = "right")
p2
#ggsave("X_type_C_vs_S_down_KEGG.png", width = 8, height = 8, units = "in", dpi = 300)
#ggsave("X_type_C_vs_S_down_KEGG.pdf", width = 8, height = 8, units = "in", dpi = 300)

DEGresult<-as.data.frame(res_type_S_vs_H)
result_up_gene<-subset(DEGresult,padj<0.05&log2FoldChange>1)
result_down_gene<-subset(DEGresult,padj<0.05& (-log2FoldChange)>1)

#up
goAll_up <- enrichGO(rownames(result_up_gene),
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "ALL",
                     pAdjustMethod= 'BH',
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.2)  
p2<-dotplot(goAll_up, font.size = 8,
            label_format= 60,
            showCategory=20)+
  ggtitle("GO analysis for upregulated genes in Plasma BS vs HC") +
  theme(plot.title = element_text(size = 8, hjust = 0.5))+
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),    
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),) +
  theme(legend.position = "right")
p2
ggsave("X_type_S_vs_H_up_KEGG.png", width = 8, height = 8, units = "in", dpi = 300)
ggsave("X_type_S_vs_H_up_KEGG.pdf", width = 8, height = 8, units = "in", dpi = 300)


#down
goAll_down <- enrichGO(rownames(result_down_gene),
                       keyType = "SYMBOL",
                       OrgDb = "org.Hs.eg.db",
                       ont = "ALL",
                       pAdjustMethod= 'BH',
                       pvalueCutoff= 0.05,
                       qvalueCutoff= 0.2) 

p2<-dotplot(goAll_down, font.size = 8,
            label_format= 60,
            showCategory=20)+
  ggtitle("GO analysis for downregulated genes in Plasma BS vs HC") +
  theme(plot.title = element_text(size = 8, hjust = 0.5))+
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),    
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),) +
  theme(legend.position = "right")
p2
ggsave("X_type_S_vs_H_down_KEGG.png", width = 8, height = 8, units = "in", dpi = 300)
ggsave("X_type_S_vs_H_down_KEGG.pdf", width = 8, height = 8, units = "in", dpi = 300)

