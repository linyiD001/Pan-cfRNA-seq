library(Biostrings)
library(BSgenome)
library(eisaR)
library(GenomicFeatures)
library(SummarizedExperiment)
library(tximeta)
library(rjson)
library(reticulate)
library(SingleCellExperiment)
library(scater)
library(BiocFileCache)
library(tximportData)

####TRY USE TXIMPORT

setwd("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_salmonstandard")

#test tximeta on process
name<-read.table("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/cfRNA.txt",header = F,sep = "\t",fill=T)
name<-as.vector(name$V1)

files <- file.path("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.standard", name, "quant.sf") 
file.exists(files)
additional_list<-name[which(file.exists(files)==TRUE)]

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

files <- file.path("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.standard", additional_list, "quant.sf") 



#coldata <- data.frame(files, names=files, condition="A", stringsAsFactors=FALSE)
#coldata
tx2gene<-read.csv("/dssg/home/acct-dahan/share/references/gencode/gencode.v44.tx.to.gene.csv",header=F)

res1 <- str_split(tx2gene[,1],"\\|")
res2 <- as.data.frame(res1)
#head(res2)
res3<-as.data.frame(t(res2)[,c(1,6)])
colnames(res3)<-c("enst","gene_symbol")



library(tximport)

txi.salmon <- tximport(files, type = "salmon", tx2gene = res3)

sampleTable <- data.frame(row.names=additional_list)
idx_S<-grep('S',rownames(sampleTable))
sampleTable$type[idx_S]<-'S'
idx_C<-grep('C',rownames(sampleTable))
sampleTable$type[idx_C]<-'C'
idx_H<-grep('H',rownames(sampleTable))
sampleTable$type[idx_H]<-'H'

dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~type)

dds 

smallestGroupSize <- 10
keep <- rowSums(counts(dds) > 0) >= smallestGroupSize
dds <- dds[keep,]
dds$type <- relevel(dds$type, ref = "H")

dds <- DESeq(dds)
resultsNames(dds)

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

nor_count_N <- assay(vst(dds,blind=F))
write.csv(nor_count_N,"nor_count_N.csv")

plotPCA(vsd, intgroup=c("type"))
# Assuming 'vsd' is your variance-stabilizing transformation object
pca_plot <- plotPCA(vsd, intgroup = c("type"), returnData = TRUE)

# Create the plot and add sample names as labels
ggplot(pca_plot, aes(x = PC1, y = PC2, color = type, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
  theme_minimal() +
  labs(title = "PCA Plot with Sample Names")


res_type_C_vs_H <- results(dds, name="type_C_vs_H")
write.csv(res_type_C_vs_H,"N_type_C_vs_H_result.csv")
type_C_vs_H_result<- lfcShrink(dds, coef ="type_C_vs_H",type = 'apeglm',res=res_type_C_vs_H)
write.csv(type_C_vs_H_result,"lfc_N_type_C_vs_H_result.csv")


res_type_S_vs_H <- results(dds, name="type_S_vs_H")
write.csv(res_type_S_vs_H,"N_type_S_vs_H_result.csv")
type_S_vs_H_result<- lfcShrink(dds, coef ="type_S_vs_H",type = 'apeglm',res=res_type_S_vs_H)
write.csv(type_S_vs_H_result,"lfc_N_type_S_vs_H_result.csv")

dds$type <- relevel(dds$type, ref = "S")

dds <- DESeq(dds)
resultsNames(dds)

res_type_C_vs_S <- results(dds, name="type_C_vs_S")
write.csv(res_type_C_vs_S,"N_type_C_vs_S_result.csv")

type_C_vs_S_result<- lfcShrink(dds, coef ="type_C_vs_S",type = 'apeglm',res=res_type_C_vs_S)
write.csv(type_C_vs_S_result,"lfc_N_type_C_vs_S_result.csv")


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
                                Significance == 'DOWN' & log2FoldChange < -4 & -log10(padj) > 2),
                  aes(label = ID),
                  segment.alpha = 0.3, #segment.size = 0.5,
                  min.segment.length = 0.2, max.overlaps = 20,
                  size = 3, color='black', segment.color = 'black') +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        strip.text = element_text(size=14, face='bold'))
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))

p


DEGresult<-as.data.frame(res_type_C_vs_S )
result_up_gene<-subset(DEGresult,padj<0.05&log2FoldChange>2)
result_down_gene<-subset(DEGresult,padj<0.05& (-log2FoldChange)>2)


goAll_up <- enrichGO(rownames(result_up_gene),
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "ALL",
                     pAdjustMethod= 'BH',
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.2)   

goAll_down <- enrichGO(rownames(result_down_gene),
                       keyType = "SYMBOL",
                       OrgDb = "org.Hs.eg.db",
                       ont = "ALL",
                       pAdjustMethod= 'BH',
                       pvalueCutoff= 0.05,
                       qvalueCutoff= 0.2)  

go_for_plot<-goAll_up

p1<-barplot(go_for_plot, font.size = 10,
            title= "GO Bar plot",
            label_format= 50,
            showCategory=20)
p1

p2<-dotplot(go_for_plot, font.size = 10,
            title= "GO Dot plot",
            label_format= 50,
            showCategory=20)
p2

go_for_plot<-goAll_down

p1<-barplot(go_for_plot, font.size = 10,
            title= "GO Bar plot",
            label_format= 50,
            showCategory=20)
p1

p2<-dotplot(go_for_plot, font.size = 10,
            title= "GO Dot plot",
            label_format= 50,
            showCategory=20)
p2


########### N or X 
idx<-grep('X',additional_list_all)
additional_list<-additional_list_all[idx]

files <- file.path("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.standard", additional_list, "quant.sf") 

#coldata <- data.frame(files, names=files, condition="A", stringsAsFactors=FALSE)
#coldata
tx2gene<-read.csv("/dssg/home/acct-dahan/share/references/gencode/gencode.v44.tx.to.gene.csv",header=F)

res1 <- str_split(tx2gene[,1],"\\|")
res2 <- as.data.frame(res1)
#head(res2)
res3<-as.data.frame(t(res2)[,c(1,6)])
colnames(res3)<-c("enst","gene_symbol")



txi.salmon <- tximport(files, type = "salmon", tx2gene = res3)

sampleTable <- data.frame(row.names=additional_list)
idx_S<-grep('S',rownames(sampleTable))
sampleTable$type[idx_S]<-'S'
idx_C<-grep('C',rownames(sampleTable))
sampleTable$type[idx_C]<-'C'
idx_H<-grep('H',rownames(sampleTable))
sampleTable$type[idx_H]<-'H'

dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~type)

dds 

smallestGroupSize <- 10
keep <- rowSums(counts(dds) > 0) >= smallestGroupSize
dds <- dds[keep,]
dds$type <- relevel(dds$type, ref = "H")

dds <- DESeq(dds)
resultsNames(dds)

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
write.csv(nor_count_X,"nor_count_X.csv")

plotPCA(vsd, intgroup=c("type"))
# Assuming 'vsd' is your variance-stabilizing transformation object
pca_plot <- plotPCA(vsd, intgroup = c("type"), returnData = TRUE)

# Create the plot and add sample names as labels
ggplot(pca_plot, aes(x = PC1, y = PC2, color = type, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
  theme_minimal() +
  labs(title = "PCA Plot with Sample Names")


res_type_C_vs_H <- results(dds, name="type_C_vs_H")
write.csv(res_type_C_vs_H,"X_type_C_vs_H_result.csv")
type_C_vs_H_result<- lfcShrink(dds, coef ="type_C_vs_H",type = 'apeglm',res=res_type_C_vs_H)
write.csv(type_C_vs_H_result,"lfc_X_type_C_vs_H_result.csv")


res_type_S_vs_H <- results(dds, name="type_S_vs_H")
write.csv(res_type_S_vs_H,"X_type_S_vs_H_result.csv")
type_S_vs_H_result<- lfcShrink(dds, coef ="type_S_vs_H",type = 'apeglm',res=res_type_S_vs_H)
write.csv(type_S_vs_H_result,"lfc_X_type_S_vs_H_result.csv")

dds$type <- relevel(dds$type, ref = "S")

dds <- DESeq(dds)
resultsNames(dds)

res_type_C_vs_S <- results(dds, name="type_C_vs_S")
write.csv(res_type_C_vs_S,"X_type_C_vs_S_result.csv")

type_C_vs_S_result<- lfcShrink(dds, coef ="type_C_vs_S",type = 'apeglm',res=res_type_C_vs_S)
write.csv(type_C_vs_S_result,"lfc_X_type_C_vs_S_result.csv")


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
                                Significance == 'DOWN' & log2FoldChange < -4 & -log10(padj) > 2),
                  aes(label = ID),
                  segment.alpha = 0.3, #segment.size = 0.5,
                  min.segment.length = 0.2, max.overlaps = 20,
                  size = 3, color='black', segment.color = 'black') +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        strip.text = element_text(size=14, face='bold'))
#theme(plot.margin =  margin(t = 0.25, r = 0.25, b = 0.25, l = 0.25, unit = "cm"))

p




DEGresult<-as.data.frame(res_type_C_vs_S )
result_up_gene<-subset(DEGresult,padj<0.05&log2FoldChange>2)
result_down_gene<-subset(DEGresult,padj<0.05& (-log2FoldChange)>2)


goAll_up <- enrichGO(rownames(result_up_gene),
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "ALL",
                     pAdjustMethod= 'BH',
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.2)   


goAll_down <- enrichGO(rownames(result_down_gene),
                       keyType = "SYMBOL",
                       OrgDb = "org.Hs.eg.db",
                       ont = "ALL",
                       pAdjustMethod= 'BH',
                       pvalueCutoff= 0.05,
                       qvalueCutoff= 0.2)    

go_for_plot<-goAll_up

p1<-barplot(go_for_plot, font.size = 10,
            title= "GO Bar plot",
            label_format= 50,
            showCategory=20)
p1

p2<-dotplot(go_for_plot, font.size = 10,
            title= "GO Dot plot",
            label_format= 50,
            showCategory=20)
p2

go_for_plot<-goAll_down

p1<-barplot(go_for_plot, font.size = 10,
            title= "GO Bar plot",
            label_format= 50,
            showCategory=20)
p1

p2<-dotplot(go_for_plot, font.size = 10,
            title= "GO Dot plot",
            label_format= 50,
            showCategory=20)
p2