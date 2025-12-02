##########downstream limma-voom in starstandard
#BiocManager::install("org.Hs.eg.db")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#devtools::install_github("GuangchuangYu/clusterProfiler")

#BiocManager::install("Glimma")

library(DESeq2)
library(apeglm)
library(regionReport)
library(ggplot2)
library(biomaRt)
library(ggrepel)
library(clusterProfiler)
library(limma)                                              
library(Glimma)
library(edgeR)
library(Homo.sapiens)

####DEG on standard
setwd("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar")
count_all<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/count.csv",header=T,row.names=1)

idx_ERCC<-grep('ERCC',rownames(count_all))
count_all<-count_all[-idx_ERCC,]

QC_list<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/QC/starstandard/QC_starstandard_list_withoutR.csv",row.names = 1)
count_analysis<-count_all[,which(colnames(count_all) %in% QC_list[,1])]

count_analysis<-as.data.frame.matrix(count_analysis)

#ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=110)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",mirror='www')

geneid=rownames(count_analysis)
mart<-useMart("ensembl")
data=listDatasets(mart)
#listFilters(ensembl)
symbols <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'),
                 filters = 'ensembl_gene_id', values = geneid,
                 mart = ensembl)
symbols

dedu_symbols<-symbols[!duplicated(symbols$ensembl_gene_id),]

dedu_symbols

#idx_NA<-which(dedu_symbols$hgnc_symbol=='')
#dedu_symbols[idx_NA,2]<-dedu_symbols[idx_NA,1]

count_analysis$ensembl_gene_id<-rownames(count_analysis)
countdata_merge<-merge(count_analysis,dedu_symbols,by='ensembl_gene_id',all.x=T)

idx_NA<-which(countdata_merge$hgnc_symbol=='')
countdata_merge[idx_NA,225]<-countdata_merge[idx_NA,1]

Anno_count <- aggregate(countdata_merge[,2:224],by=list(countdata_merge$hgnc_symbol), FUN=sum)
#write.csv(Anno_count,"star_standard_anno_count.csv")

count_analysis<-Anno_count[,-1]
rownames(count_analysis)<-Anno_count[,1]

####误诊(此处进行了排除)
#idx_R<-grep('C89',colnames(count_analysis))
#count_analysis<-count_analysis[,-idx_R]
#idx_R<-grep('C103',colnames(count_analysis))
#count_analysis<-count_analysis[,-idx_R]
#idx_R<-grep('C107',colnames(count_analysis))
#count_analysis<-count_analysis[,-idx_R]

idx_N<-grep('N',colnames(count_analysis))
count_N<-count_analysis[,idx_N]

idx_X<-grep('X',colnames(count_analysis))
count_X<-count_analysis[,idx_X]

####comparing X and N (paired)

count_analysis_all<-count_analysis

sample_name <-colnames(count_analysis)
sample_name <- sub(paste0("*", "_X"), "", sample_name)
sample_name <- sub(paste0("*", "_N"), "", sample_name)
sample_name <-as.factor(sample_name)
sample_name

coldata = data.frame(row.names = colnames(count_analysis_all))
num_paired<-duplicated(sample_name) | duplicated(sample_name, fromLast = TRUE)
paired<-as.vector(num_paired)
idx_paired<- grep(TRUE, paired)
count_analysis_paired<-count_analysis[,idx_paired]

#count_analysis_paired <- count_analysis_paired[,grep("H",colnames(count_analysis_paired))]
coldata_paired = data.frame(row.names = colnames(count_analysis_paired))
sample_name <-colnames(count_analysis_paired)
sample_name <- sub(paste0("*", "_X"), "", sample_name)
sample_name <- sub(paste0("*", "_N"), "", sample_name)
sample_name <-as.factor(sample_name)
sample_name

idx_X<-grep('X',colnames(count_analysis_paired))
coldata_paired$type[idx_X]<-'X'
idx_N<-grep('N',colnames(count_analysis_paired))
coldata_paired$type[idx_N]<-'N'
coldata_paired$sample<-sample_name

dds <- DESeqDataSetFromMatrix(countData = round(count_analysis_paired),
                              colData = coldata_paired,
                              design = ~ sample + type)
dds

smallestGroupSize <- 10
keep <- rowSums(counts(dds) > 5) >= smallestGroupSize
dds <- dds[keep,]
dds$type <- relevel(dds$type, ref = "N")

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
colnames(sampleDistMatrix) <- colnames(count_analysis_paired)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
#nor_count_N<-vst(dds, blind=T)

save.image("DEG_starstandard.RData")
plotPCA(vsd,ntop=1000, intgroup=c("type"))


# 创建 PCA 数据
pcaData <- plotPCA(vsd,ntop=1000,  intgroup = c("type"), returnData = TRUE)

# 自定义颜色
colors <- c("X" = "pink", "N" = "lightblue")

# 绘制 PCA 图
ggplot(pcaData, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 3) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  stat_ellipse(aes(fill = type, color = type), geom = "polygon", alpha = 0.15) +  # 添加半透明圆圈
  theme_bw() +
  labs(title = "PCA Plot",
       x = "PC1 81%",
       y = "PC2 7%") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank())

ggsave("NX_PCA_plot.png", width = 6, height = 6, dpi = 600)
ggsave("NX_PCA_plot.pdf", width = 6, height = 6, dpi = 600)

#BiocManager::install("umap")
library(umap)
# 提取归一化后的表达矩阵
mat <- assay(vsd)

# 运行 UMAP
set.seed(123)
topVarGenes <- head(order(rowVars(mat), decreasing=TRUE), 5000)
map_data <- umap(t(mat[topVarGenes, ]),
                 n_neighbors = 30,
                 min_dist = 0.3,
                 metric = "euclidean",
                 random_state = 123)
# 合并结果和样本信息
umap_df <- data.frame(map_data$layout, colData(vsd))
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")
umap_df$type <- coldata_paired$type

# 画图
ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=type)) +
  geom_point(size=3, alpha=0.8) +
  theme_classic(base_size = 14) +
  labs(title="RNA-seq UMAP", color="Group") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("NX_UMAP_plot.png", width = 8, height = 8, dpi = 600)
ggsave("NX_UMAP_plot.pdf", width = 8, height = 8, dpi = 600)


res_type_X_vs_N <- results(dds, name="type_X_vs_N")
hist(res_type_X_vs_N$pvalue, col = "grey50", border = "white", main = "", xlab = "p-values")
res_type_X_vs_N_res <- as.data.frame(res_type_X_vs_N)
#res_type_X_vs_N_res$ensembl_gene_id <- rownames(res_type_X_vs_N)
#res_type_X_vs_N_res <- merge(res_type_X_vs_N_res,dedu_symbols,by='ensembl_gene_id')
write.csv(res_type_X_vs_N_res,"type_X_vs_N_result.csv")

type_X_vs_N_result<- lfcShrink(dds, coef ="type_X_vs_N",type = 'apeglm',res=res_type_X_vs_N)
DESeq2::plotMA(type_X_vs_N_result)

type_X_vs_N_result_res <- as.data.frame(type_X_vs_N_result)
#type_X_vs_N_result_res$ensembl_gene_id <- rownames(type_X_vs_N_result)
#type_X_vs_N_result_res <- merge(type_X_vs_N_result_res,dedu_symbols,by='ensembl_gene_id')
write.csv(type_X_vs_N_result,"lfc_type_X_vs_N_result.csv")

dataForVolcanoPlot <- as.data.frame(type_X_vs_N_result_res)

logFcThreshold <- log2(2)
adjPvalThreshold <- 0.05
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange < logFcThreshold | padj> adjPvalThreshold)] <- 'NS'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange >= logFcThreshold & padj <= adjPvalThreshold)] <- 'UP'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange <= -logFcThreshold & padj <= adjPvalThreshold)] <- 'DOWN'


dataForVolcanoPlot
dataForVolcanoPlot$ID = rownames(type_X_vs_N_result_res)

#cols <- c('UP'=="red", 'NS'='gray','DOWN'="green")
cols <- c('UP'="pink", 'NS'='gray','DOWN'="lightblue")


dev.new()

p <- ggplot(dataForVolcanoPlot, aes(x = log2FoldChange, y = -log10(padj))) +
  #xlim(-2,2) +
  labs(x=expression(bold('Log'['2']*'(Fold Change)')), 
       y=(expression(bold('-Log'['10']*'(FDR)'))), 
       title=NULL) +
  geom_point(color='gray', alpha=1, size=1) +
  geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
             color='darkgreen', linetype='dashed') +
  geom_hline(yintercept = -log10(adjPvalThreshold), 
             color='darkgreen',linetype='dashed')+
  scale_color_manual(values = cols) +
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'UP'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="pink", 
             size=1, color="pink",shape=21)+
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'DOWN'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="lightblue", 
             size=1, color="lightblue",shape=21) +
  scale_fill_manual(values = cols) +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'UP' & (-log10(padj) > 10 & log2FoldChange > 4)),
                  aes(label = ID), nudge_x = 0.1,
                  segment.alpha = 0.5, #segment.size = 0.5,
                  min.segment.length = 0.2,
                  size = 3, color='black', segment.color = 'black') +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'DOWN' & log2FoldChange < -4 & -log10(padj) > 30),
                  aes(label = ID),
                  segment.alpha = 0.5, #segment.size = 0.5,
                  min.segment.length = 0.2, max.overlaps = 10,
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

ggsave("N_X_volcano.png", width = 6, height = 6, dpi = 600)
ggsave("N_X_volcano.pdf", width = 6, height = 6, dpi = 600)

analysisresult<-as.data.frame(res_type_X_vs_N_res)
logFCThreshold=1
adj.P.ValThreshold=0.05
result_up_gene_list<-subset(analysisresult,padj<adj.P.ValThreshold&log2FoldChange>logFCThreshold)
result_down_gene_list<-subset(analysisresult,padj<adj.P.ValThreshold&log2FoldChange<(-logFCThreshold))
#result_DEG<-subset(analysisresult,padj<adj.P.ValThreshold&abs(log2FoldChange)>logFCThreshold)
load("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/cell_type_marker/database/humanGeneSets.RData")
result_up_gene<-rownames(result_up_gene_list)
result_down_gene<-rownames(result_down_gene_list)
#up
goAll_up <- enrichGO(result_up_gene,
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "ALL",
                     pAdjustMethod= 'BH',
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.2)  # 以上三种  


#down
goAll_down <- enrichGO(result_down_gene,
                       keyType = "SYMBOL",
                       OrgDb = "org.Hs.eg.db",
                       ont = "ALL",
                       pAdjustMethod= 'BH',
                       pvalueCutoff= 0.05,
                       qvalueCutoff= 0.2)  # 以上三种  
go_for_plot<-goAll_up
p1<-dotplot(go_for_plot, font.size = 10,
            title= "Up in Urine cfRNA",
            label_format= 50,
            showCategory=10)
p1
ggsave("N_up_GO_All.png", width = 8, height = 6, dpi = 600)
ggsave("N_up_GO_All.pdf", width = 8, height = 6, dpi = 600)

go_for_plot<-goAll_down
p2<-dotplot(go_for_plot, font.size = 10,
            title= "Up in Plasma cfRNA",
            label_format= 50,
            showCategory=10)
p2
ggsave("X_up_GO_All.png", width = 8, height = 6, dpi = 600)
ggsave("X_up_GO_All.pdf", width = 8, height = 6, dpi = 600)


tissue_significance <- read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/cell_type_marker/database/PanglaoDB_markers_27_Mar_2020.csv")
tissue_significance_human <- tissue_significance %>%
  dplyr::filter(species == "Hs") %>%
  dplyr::group_by(cell.type) %>%
  dplyr::summarise(genes = list(unique(official.gene.symbol))) %>%
  dplyr::filter(lengths(genes) >= 5)  # 

library(clusterProfiler)
type_X_vs_N_result_genelist<-type_X_vs_N_result_res#[res_type_X_vs_N_res$padj<0.05,]
geneList_pre<-cbind(rownames(type_X_vs_N_result_genelist),type_X_vs_N_result_genelist$log2FoldChange)
geneList = as.numeric(geneList_pre[,2])
names(geneList) = as.character(geneList_pre[,1])
geneList = sort(geneList, decreasing = TRUE)
geneList <- geneList[names(geneList) != ""]  
geneList<- geneList[!is.na(names(geneList))]

TERM2GENE <- tissue_significance_human[,] %>%
  tidyr::unnest(genes) %>%
  dplyr::select(term = cell.type, gene = genes)

gsea_result <- GSEA(
  geneList = geneList,
  TERM2GENE = TERM2GENE,  # 把 list 转成两列：term, gene
  pvalueCutoff = 0.05,
  minGSSize = 5,
  maxGSSize = 500,
  verbose = FALSE
)
gsea_result

# 保存为 PNG
png("GSEA_DC.png", width = 6, height = 6,res=300,unit='in')
gseaplot(
  gsea_result,           
  geneSetID = 1,         
  title = gsea_result@result$Description[1],  
  by = "all"             
)
dev.off()

# 保存为 PDF
png("GSEA_DC.pdf", width = 6, height = 6,res=300,unit='in')
gseaplot(
  gsea_result,           
  geneSetID = 1,         
  title = gsea_result@result$Description[1],  
  by = "all"             
)
dev.off()


# 保存为 PNG
png("GSEA_2.png", width = 8, height = 8,res=300,unit='in')
gseaplot(
  gsea_result,           
  geneSetID = 2,         
  title = gsea_result@result$Description[2],  
  by = "all"             
)
dev.off()

# 保存为 PDF
png("GSEA_2.pdf", width = 8, height = 8,res=300,unit='in')
gseaplot(
  gsea_result,           
  geneSetID = 2,         
  title = gsea_result@result$Description[2],  
  by = "all"             
)
dev.off()

head(humanGeneSets)[1:3]

TERM2GENE <- do.call(
  rbind,
  lapply(colnames(humanGeneSets), function(term) {
    genes <- rownames(humanGeneSets)[humanGeneSets[, term] != 0]
    if (length(genes) > 0) {
      data.frame(term = term, gene = genes, stringsAsFactors = FALSE)
    }
  })
)

gsea_result <- GSEA(
  geneList = geneList,
  TERM2GENE = TERM2GENE,  # 把 list 转成两列：term, gene
  pvalueCutoff = 1,
  minGSSize = 5,
  maxGSSize = 500,
  verbose = FALSE,
  eps =1e-200
)
load("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/cell_type_marker/database/humanGeneSetsInfo.RData")
head(humanGeneSets)
head(gsea_result@result, 10)

res <- gsea_result@result
colnames(humanGeneSetsInfo)[colnames(humanGeneSetsInfo) == "gset"] <- "ID"

res_annotated <- merge(res, humanGeneSetsInfo, by = "ID", all.x = TRUE)

plot_df <- res_annotated %>%
  mutate(logP = -log10(p.adjust)) %>%
  arrange(gsetBioName, desc(logP)) %>%
  mutate(x = row_number())

# 每个类别的边界
category_boundaries <- plot_df %>%
  group_by(gsetBioName) %>%
  summarise(start = min(x), end = max(x)) %>%
  mutate(mid = (start + end) / 2)

# 在每个类别区间内随机生成点
set.seed(1234)
plot_df <- plot_df %>%
  left_join(category_boundaries, by = "gsetBioName") %>%
  rowwise() %>%
  mutate(x_jitter = runif(1, start, end)) %>%
  ungroup()

# 标注每组前3个显著集
plot_df <- plot_df %>%
  group_by(gsetBioName) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  mutate(significance = if_else(row_number() <= 3, "anno", "NS")) %>%
  ungroup()

ymax= max(plot_df$logP)+5
# 绘图

ggplot(plot_df, aes(x = x_jitter, y = logP, color = gsetBioName, size = setSize)) +
  ylim(0, ymax) +
  geom_rect(
    data = category_boundaries,
    inherit.aes = FALSE,
    aes(xmin = start, xmax = end, ymin = 0, ymax = ymax, fill = gsetBioName),
    alpha = 0.08, color = NA
  ) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(
    data = category_boundaries[-1, ],
    aes(xintercept = start - 1),
    color = "gray50", linetype = "dotted", linewidth = 0.4
  ) +
  geom_text(
    data = category_boundaries,
    aes(x = mid, y = -0.5, label = gsetBioName, color = gsetBioName),
    inherit.aes = FALSE,
    size = 3.5, angle = 30, hjust = 1
  ) +
  scale_size_continuous(range = c(0.5, 2)) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  labs(
    y = expression(-log[10](p.adjust)),
    size = "Gene set size",
    color = "Category",
    fill = "Category"
  ) +
  geom_text_repel(
    data = subset(plot_df, significance == "anno"),
    aes(label = ID),             # 不在这里改颜色
    color = ifelse(subset(plot_df, significance=="anno")$NES > 0,"#DB498E","darkblue"),  # 文字颜色单独指定
    segment.alpha = 1,
    text.alpha=1,
    min.segment.length = 0.2,
    max.overlaps = 5,
    size = 3.5,
    face="bold",
    segment.color = "black"
  )


ggsave("N_X_GSEA_summary_plot.png", width = 15, height = 6, dpi = 600)
ggsave("N_X_GSEA_summary_plot.pdf", width = 15, height = 6, dpi = 600)


library(pheatmap)
library(RColorBrewer)
type_X_vs_N_result <- read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/lfc_type_X_vs_N_result.csv")
# 提取结果表
res_df <- as.data.frame(type_X_vs_N_result)

# 筛选显著基因（可根据 padj 或 logFC）
sig_res <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
sig_genes <- sig_res$X

# 若显著基因太多，可取前50或前100个
top_genes <- head(sig_genes[order(sig_res$padj)],12)

# 获取这些基因的标准化表达量

mat <- as.data.frame(assay(vsd))
mat <-mat[top_genes, ]
head(mat)
# Z-score 标准化每个基因（行标准化）
mat_z <- t(scale(t(mat)))

# 构建分组注释
annotation_col <- data.frame(
  Type = colData(vsd)$type
)
rownames(annotation_col) <- colnames(vsd)
ann_colors <- list(Type = c("N" = "lightblue", "X" = "pink"))

# 绘制热图
mat_z_t <- t(mat_z)

#png("NX_top12_heatmap.png", width =6 , height =8, unit='in', res = 300)
pdf("NX_top12_heatmap.pdf", width =6 , height =8)

pheatmap(mat_z_t,
         annotation_row = annotation_col,   # 注意这里改成 row 注释
         annotation_colors = ann_colors,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(255),
         show_rownames = FALSE,
         show_colnames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 10,
         fontsize_col = 8,
         scale='column'
         #main = "Top Differentially Expressed Genes (Rotated)"
         )
dev.off()




# 绘制条形图
ggplot(gsea_df, aes(x = NES, y = ID, fill = NES)) +
  geom_col(width = 0.6) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_bw(base_size = 12) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "",
    title = "GSEA of Cell-type specific genes"
  ) +
  theme(
    axis.text.y = element_text(, size = 10),
    plot.title = element_text(, hjust = 0.5)
  )
ggsave("N_X_GSEA.png", width = 10, height = 6, dpi = 600)
ggsave("N_X_GSEA.pdf", width = 10, height = 6, dpi = 600)





####N only

count_analysis_N<-count_N

coldata_N = data.frame(row.names = colnames(count_analysis_N))
idx_S<-grep('S',colnames(count_analysis_N))
coldata_N$type[idx_S]<-'BS'
idx_C<-grep('C',colnames(count_analysis_N))
coldata_N$type[idx_C]<-'BC'
idx_H<-grep('H',colnames(count_analysis_N))
coldata_N$type[idx_H]<-'HC'


data <-count_N

coldata<-matrix(ncol = 2,nrow = 106)
rownames(coldata)<-colnames(data)

idx_N_counts<-grep("N",colnames(data))
idx_S_counts<-grep("S",colnames(data))
idx_H_counts<-grep("H",colnames(data))
idx_C_counts<-grep("C",colnames(data))

coldata[idx_N_counts,1]<-"Urine"
coldata[idx_H_counts,2]<-"HC"
coldata[idx_S_counts,2]<-"BS"
coldata[idx_C_counts,2]<-"BC"


colnames(coldata)<-c("category","Group")

library(readxl)

BC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BC")
BC_cd<-as.data.frame(cbind(BC_clinicial_data$Patient,BC_clinicial_data$MIBC,BC_clinicial_data$AGE,BC_clinicial_data$SEX,
                           BC_clinicial_data$Tumor_size,BC_clinicial_data$`T`,BC_clinicial_data$Tumor_number,BC_clinicial_data$Grade))
colnames(BC_cd)<-c("Patients","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade')

BS_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BS")
BS_cd<-as.data.frame(cbind(BS_clinicial_data$Patient,NA,BS_clinicial_data$AGE,BS_clinicial_data$SEX,NA,NA,NA,NA))
colnames(BS_cd)<-c("Patients","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade')

HC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "HC")
HC_cd<-as.data.frame(cbind(HC_clinicial_data$Patient,NA,HC_clinicial_data$AGE,HC_clinicial_data$SEX,NA,NA,NA,NA))
colnames(HC_cd)<-c("Patients","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade')


clinicial_design<-rbind(BC_cd,BS_cd,HC_cd)

coldata<-as.data.frame(coldata)
coldata$sample_num<-rownames(coldata)
coldata$Patients<-rownames(coldata)
coldata$Patients <- sub(paste0("*", "_X"), "", coldata$Patients)
coldata$Patients <- sub(paste0("*", "_N"), "", coldata$Patients)


coldata_all<-merge(coldata,clinicial_design,by="Patients",all.x=T)
rownames(coldata_all)<-coldata_all$sample_num
coldata_all<-coldata_all[,-1]
coldata_all<-coldata_all[,-3]

coldata_all <- coldata_all[colnames(data), , drop = FALSE]
dds <- DESeqDataSetFromMatrix(countData = round(data),
                              colData = coldata_all,
                              design = ~ Sex)
dds

smallestGroupSize <- 10
keep <- rowSums(counts(dds) > 5) >= smallestGroupSize
dds <- dds[keep,]
#dds$type <- relevel(dds$type, ref = "HC")

dds <- DESeq(dds)
resultsNames(dds)

ntd <- normTransform(dds)
library("pheatmap")
df <- as.data.frame(colData(dds)[,c("Sex")])
vsd <- vst(dds, blind=FALSE)
#rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Sex
colnames(sampleDistMatrix) <- colnames(count_analysis_N)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
#nor_count_N<-vst(dds, blind=T)


plotPCA(vsd,ntop=1000, intgroup=c("Sex"))

# 创建 PCA 数据
pcaData <- plotPCA(vsd,ntop=1000,  intgroup = c("Sex"), returnData = TRUE)

# 自定义颜色
colors <- c("Female" = "pink", "Male" = "lightblue")

# 绘制 PCA 图
ggplot(pcaData, aes(x = PC1, y = PC2, color = Sex)) +
  geom_point(size = 3) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  stat_ellipse(aes(fill = Sex, color = Sex), geom = "polygon", alpha = 0.15) +  # 添加半透明圆圈
  theme_bw() +
  labs(title = "PCA Plot",
       x = "PC1 49%",
       y = "PC2 9%") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank())
ggsave("PCA_for_gender.png",height=6,width=6,unit='in',res=300)
ggsave("PCA_for_gender.pdf",height=6,width=6)


##备用编组
#idx_R<-grep('C89',colnames(count_analysis_N))
#coldata$type[idx_R]<-'R'
#idx_R<-grep('C103',colnames(count_analysis_N))
#coldata$type[idx_R]<-'R'
#idx_R<-grep('C107',colnames(count_analysis_N))
#coldata$type[idx_R]<-'R'


dds <- DESeqDataSetFromMatrix(countData = round(count_analysis_N),
                              colData = coldata_N,
                              design = ~ type)
dds

smallestGroupSize <- 10
keep <- rowSums(counts(dds) > 5) >= smallestGroupSize
dds <- dds[keep,]
dds$type <- relevel(dds$type, ref = "HC")

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
colnames(sampleDistMatrix) <- colnames(count_analysis_N)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
#nor_count_N<-vst(dds, blind=T)


plotPCA(vsd,ntop=1000, intgroup=c("type"))

# 提取绘图所用的矩阵
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing=TRUE)[seq_len(min(1000, length(rv)))]

# 对选出的基因做 PCA
pca <- prcomp(t(assay(vsd)[select, ],scale=F))

# 提取 loadings (每个基因在每个PC上的权重)
loadings <- pca$rotation

# 找出最影响 PC1 的基因
pc1_load <- sort(abs(loadings[,1]), decreasing = TRUE)

# 查看前 20 个对 PC1 影响最大的基因
head(pc1_load, 20)

# 创建 PCA 数据
pcaData <- plotPCA(vsd, ntop=1000 ,intgroup = c("type"), returnData = TRUE)

cpm <- read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/N_cpm_matrix_from_star.csv",header=T,row.names=1)
MT_genes <- cpm[grep("MT-",rownames(cpm)),]
log2_MT_genes <- as.data.frame(log2(MT_genes+1))
mean_log2_MT_genes <- colSums(log2_MT_genes)
mean_log2_MT_genes

data_combine<-cbind(pcaData$PC1,as.data.frame(mean_log2_MT_genes))
cor_value <- cor(pcaData$PC1, mean_log2_MT_genes, method = "pearson")
cat("相关系数 r =", round(cor_value, 3), "\n")
data_combine$ID[grepl("C", rownames(data_combine))] <- "BC"
data_combine$ID[grepl("S", rownames(data_combine))] <- "BS"
data_combine$ID[grepl("H", rownames(data_combine))] <- "HC"



ggplot(data_combine, aes(x = pcaData$PC1, y = mean_log2_MT_genes)) +
  stat_cor(method = "pearson", 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           size = 5) +
  geom_point(size = 3, alpha = 0.7,aes(color=ID)) +       # 散点
  geom_smooth(method = "lm", se = TRUE,     # 拟合线（线性回归）
              color = "blue", linetype = "solid") +
  labs(x = "PC1", y = "mean_log2_MT_genes")+ 
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13)
  )


# 自定义颜色
colors <- c("HC" = "#379FB4", "BS" = "#FCAE59", "BC" = "#DB498E")

# 绘制 PCA 图
p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 3) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  stat_ellipse(aes(fill = type, color = type), geom = "polygon", alpha = 0.15) +  # 添加半透明圆圈
  theme_bw() +
  labs(
    title = "PCA Plot of Urine cfRNA",
    x = "PC1 51%",
    y = "PC2 8%"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )

p
ggsave("N_PCA.png", width = 8, height = 8, dpi = 600)
ggsave("N_PCA.pdf", width = 8, height = 8, dpi = 600)
#BiocManager::install("umap")
library(umap)
# 提取归一化后的表达矩阵
mat <- assay(vsd)

# 运行 UMAP
set.seed(123)
topVarGenes <- head(order(rowVars(mat), decreasing=TRUE), 1000)
map_data <- umap(t(mat[topVarGenes, ]),
                  n_neighbors = 30,
                  min_dist = 0.3,
                  metric = "euclidean",
                  random_state = 123)
# 合并结果和样本信息
umap_df <- data.frame(map_data$layout, colData(vsd))
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")
umap_df$type <- coldata_N$type

# 画图
ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=type)) +
  geom_point(size=3, alpha=0.8) +
  theme_classic(base_size = 14) +
  labs(title="", color="Group") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("N_UMAP.png", width = 8, height = 8, dpi = 600)
ggsave("N_UMAP.pdf", width = 8, height = 8, dpi = 600)


nor_count_N <- assay(vst(dds))
write.csv(nor_count_N ,"nor_N_count.csv")

resultsNames(dds)
res_type_C_vs_H <- results(dds, name="type_BC_vs_HC")
hist(res_type_C_vs_H$pvalue, col = "grey50", border = "white", main = "", xlab = "p-values")
res_type_C_vs_H_res <- as.data.frame(res_type_C_vs_H)
#res_type_C_vs_H_res$ensembl_gene_id <- rownames(res_type_C_vs_H)
#res_type_C_vs_H_res <- merge(res_type_C_vs_H_res,dedu_symbols,by='ensembl_gene_id')
write.csv(res_type_C_vs_H_res,"N_type_C_vs_H_result.csv")

type_C_vs_H_result<- lfcShrink(dds, coef ="type_BC_vs_HC",type = 'apeglm',res=res_type_C_vs_H)
type_C_vs_H_result_res <- as.data.frame(type_C_vs_H_result)
DESeq2::plotMA(type_C_vs_H_result)
#type_C_vs_H_result_res$ensembl_gene_id <- rownames(type_C_vs_H_result)
#type_C_vs_H_result_res <- merge(type_C_vs_H_result_res,dedu_symbols,by='ensembl_gene_id')
write.csv(type_C_vs_H_result,"lfc_N_type_C_vs_H_result.csv")

res_type_S_vs_H <- results(dds, name="type_BS_vs_HC")
hist(res_type_S_vs_H$pvalue, col = "grey50", border = "white", main = "", xlab = "p-values")
res_type_S_vs_H_res <- as.data.frame(res_type_S_vs_H)
#res_type_S_vs_H_res$ensembl_gene_id <- rownames(res_type_S_vs_H)
#res_type_S_vs_H_res <- merge(res_type_S_vs_H_res,dedu_symbols,by='ensembl_gene_id')
write.csv(res_type_S_vs_H_res,"N_type_S_vs_H_result.csv")

type_S_vs_H_result<- lfcShrink(dds, coef ="type_BS_vs_HC",type = 'apeglm',res=res_type_S_vs_H)
type_S_vs_H_result_res <- as.data.frame(type_S_vs_H_result)
#type_S_vs_H_result_res$ensembl_gene_id <- rownames(type_S_vs_H_result)
#type_S_vs_H_result_res <- merge(type_S_vs_H_result_res,dedu_symbols,by='ensembl_gene_id')
write.csv(type_S_vs_H_result,"lfc_N_type_S_vs_H_result.csv")


dds$type <- relevel(dds$type, ref = "BS")

dds <- DESeq(dds)
resultsNames(dds)

res_type_C_vs_S <- results(dds, name="type_BC_vs_BS")
hist(res_type_C_vs_S$pvalue, col = "grey50", border = "white", main = "", xlab = "p-values")
res_type_C_vs_S_res <- as.data.frame(res_type_C_vs_S)
#res_type_C_vs_S_res$ensembl_gene_id <- rownames(res_type_C_vs_S)
#res_type_C_vs_S_res <- merge(res_type_C_vs_S_res,dedu_symbols,by='ensembl_gene_id')
write.csv(res_type_C_vs_S_res,"N_type_C_vs_S_result.csv")

type_C_vs_S_result<- lfcShrink(dds, coef ="type_BC_vs_BS",type = 'apeglm',res=res_type_C_vs_S)
type_C_vs_S_result_res <- as.data.frame(type_C_vs_S_result)
#type_C_vs_S_result_res$ensembl_gene_id <- rownames(type_C_vs_S_result)
#type_C_vs_S_result_res <- merge(type_C_vs_S_result_res,dedu_symbols,by='ensembl_gene_id')
write.csv(type_C_vs_S_result,"lfc_N_type_C_vs_S_result.csv")



dataForVolcanoPlot <- as.data.frame(type_S_vs_H_result_res)

logFcThreshold <- log2(2)
adjPvalThreshold <- 0.05
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange < logFcThreshold | padj> adjPvalThreshold)] <- 'NS'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange >= logFcThreshold & padj <= adjPvalThreshold)] <- 'UP'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange <= -logFcThreshold & padj <= adjPvalThreshold)] <- 'DOWN'


dataForVolcanoPlot
dataForVolcanoPlot$ID = rownames(dataForVolcanoPlot)

#cols <- c('UP'=="red", 'NS'='gray','DOWN'="green")
cols <- c('UP'="pink", 'NS'='gray','DOWN'="lightblue")


dev.new()

p <- ggplot(dataForVolcanoPlot, aes(x = log2FoldChange, y = -log10(padj))) +
  #xlim(-2,2) +
  labs(x=expression(bold('Log'['2']*'(Fold Change)')), 
       y=(expression(bold('-Log'['10']*'(FDR)'))), 
       title=NULL) +
  geom_point(color='gray', alpha=1, size=1) +
  geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
             color='darkgreen', linetype='dashed') +
  geom_hline(yintercept = -log10(adjPvalThreshold), 
             color='darkgreen',linetype='dashed')+
  scale_color_manual(values = cols) +
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'UP'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="pink", 
             size=1, color="pink",shape=21)+
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'DOWN'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="lightblue", 
             size=1, color="lightblue",shape=21) +
  scale_fill_manual(values = cols) +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'UP' & (-log10(padj) > 2 & log2FoldChange > 1)),
                  aes(label = ID), nudge_x = 0.1,
                  segment.alpha = 0.5, #segment.size = 0.5,
                  min.segment.length = 0.2,
                  size = 3, color='black', segment.color = 'black') +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'DOWN' & log2FoldChange < -1 & -log10(padj) > 2),
                  aes(label = ID),
                  segment.alpha = 0.5, #segment.size = 0.5,
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

ggsave("N_S_H_volcano.png", width = 12, height = 8, dpi = 600)
ggsave("N_S_H_volcano.pdf", width = 12, height = 8, dpi = 600)


DEGresult<-as.data.frame(type_C_vs_H_result)
#logFCThreshold=1
#adj.P.ValThreshold=0.05
#result_up_gene<-subset(DEGresult,padj<adj.P.ValThreshold&log2FoldChange>logFCThreshold)
#result_down_gene<-subset(DEGresult,padj<adj.P.ValThreshold&log2FoldChange<(-logFCThreshold))
#result_DEG<-subset(DEGresult,padj<adj.P.ValThreshold&abs(log2FoldChange)>logFCThreshold)


result_DEG<-subset(DEGresult,padj<0.05&abs(log2FoldChange)>2)

select<-which(rownames(dds) %in% rownames(result_DEG))
df <- as.data.frame(coldata_N)


analysisresult<-as.data.frame(type_S_vs_H_result_res)
logFCThreshold=1
adj.P.ValThreshold=0.05
result_up_gene_list<-subset(analysisresult,padj<adj.P.ValThreshold&log2FoldChange>logFCThreshold)
result_down_gene_list<-subset(analysisresult,padj<adj.P.ValThreshold&log2FoldChange<(-logFCThreshold))
#result_DEG<-subset(analysisresult,padj<adj.P.ValThreshold&abs(log2FoldChange)>logFCThreshold)

result_up_gene<-rownames(result_up_gene_list)
result_down_gene<-rownames(result_down_gene_list)
#up
#goMF_up <- enrichGO(result_up_gene,
#                    keyType = "SYMBOL",    # 类型
#                    OrgDb = "org.Hs.eg.db",  # 
#                    ont = "MF",
#                    pAdjustMethod= 'BH',
#                    pvalueCutoff= 0.05,
#                    qvalueCutoff= 0.2)    # 分子功能
#goBP_up <- enrichGO(result_up_gene,
#                    keyType = "SYMBOL",
#                    OrgDb = "org.Hs.eg.db",
#                    ont = "BP",
#                    pAdjustMethod= 'BH',
#                    pvalueCutoff= 0.05,
#                    qvalueCutoff= 0.2)  # 生物学过程
#goCC_up <- enrichGO(result_up_gene,
#                    keyType = "SYMBOL",
#                    OrgDb = "org.Hs.eg.db",
#                    ont = "CC",
#                    pAdjustMethod= 'BH',
#                    pvalueCutoff= 0.05,
#                    qvalueCutoff= 0.2) # 细胞组分
goAll_up <- enrichGO(result_up_gene,
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "ALL",
                     pAdjustMethod= 'BH',
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.2)  # 以上三种  


#down
#goMF_down <- enrichGO(result_down_gene,
#                      keyType = "SYMBOL",    # degEntrezIdl类型
#                      OrgDb = "org.Hs.eg.db",  # 人的数据库
#                      ont = "MF",
#                      pAdjustMethod= 'BH',
#                      pvalueCutoff= 0.05,
#                      qvalueCutoff= 0.2)   # 分子功能
#goBP_down <- enrichGO(result_down_gene,
#                      keyType = "SYMBOL",
#                      OrgDb = "org.Hs.eg.db",
#                      ont = "BP",
#                      pAdjustMethod= 'BH',
#                      pvalueCutoff= 0.05,
#                      qvalueCutoff= 0.2)  # 生物学过程
#goCC_down <- enrichGO(result_down_gene,
#                      keyType = "SYMBOL",
#                      OrgDb = "org.Hs.eg.db",
#                      ont = "CC",
#                      pAdjustMethod= 'BH',
#                      pvalueCutoff= 0.05,
#                      qvalueCutoff= 0.2) # 细胞组分
goAll_down <- enrichGO(result_down_gene,
                       keyType = "SYMBOL",
                       OrgDb = "org.Hs.eg.db",
                       ont = "ALL",
                       pAdjustMethod= 'BH',
                       pvalueCutoff= 0.05,
                       qvalueCutoff= 0.2)  # 以上三种  

go_for_plot <- goAll_up
library(enrichplot)
p2 <- dotplot(go_for_plot, font.size = 10,
            title= "GO Dot plot",
            label_format= 50,
            showCategory=10)
p2

ggsave("N_S_H_up_GO.png", width = 6, height = 4, dpi = 600)
ggsave("N_S_H_up_GO.pdf", width = 6, height = 4, dpi = 600)

go_for_plot <- goAll_down
library(enrichplot)
p3 <- dotplot(go_for_plot, font.size = 10,
            title= "GO Dot plot",
            label_format= 50,
            showCategory=10)
p3

ggsave("N_S_H_down_GO.png", width = 6, height = 4, dpi = 600)
ggsave("N_S_H_down_GO.pdf", width = 6, height = 4, dpi = 600)




geneList_pre<-cbind(rownames(type_C_vs_S_result_res),type_C_vs_S_result_res$log2FoldChange)
geneList = as.numeric(geneList_pre[,2])
names(geneList) = as.character(geneList_pre[,1])
geneList = sort(geneList, decreasing = TRUE)
geneList<- geneList[!is.na(names(geneList))]

ego3 <- gseGO(geneList     = geneList,
              OrgDb        = "org.Hs.eg.db",
              keyType      = "SYMBOL",
              ont          = "ALL",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

ridgeplot(ego3,
          showCategory = 10,      # 显示前多少条通路
          fill = "p.adjust",      # 填充颜色，常用 p.adjust / NES
          label_format = 50       # 通路名称换行宽度
)
ggsave("N_C_H_gseGO.png", width = 8, height = 4, dpi = 600)
ggsave("N_C_H_gseGO.pdf", width = 8, height = 4, dpi = 600)




gsea_result <- GSEA(
  geneList = geneList,
  TERM2GENE = TERM2GENE,  # 把 list 转成两列：term, gene
  pvalueCutoff = 0.05,
  minGSSize = 5,
  maxGSSize = 500,
  verbose = FALSE
)
gsea_result
head(gsea_result@result, 10)
gsea_df <- gsea_result@result %>%
  as.data.frame() %>%
  dplyr::mutate(ID = factor(ID, levels = ID[order(NES)]))  # 按 NES 排序


# 绘制条形图
ggplot(gsea_df, aes(x = NES, y = ID, fill = NES)) +
  geom_col(width = 0.6) +
  #scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  #scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_bw(base_size = 12) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "",
    title = "GSEA of Cell-type specific genes"
  ) +
  theme(
    axis.text.y = element_text(, size = 10),
    plot.title = element_text(, hjust = 0.5)
  )
ggsave("N_CH_GSEA.png", width = 10, height = 6, dpi = 600)
ggsave("N_CH_GSEA.pdf", width = 10, height = 6, dpi = 600)



####X only

count_analysis_X<-count_X

coldata_X = data.frame(row.names = colnames(count_analysis_X))
idx_S<-grep('S',colnames(count_analysis_X))
coldata_X$type[idx_S]<-'BS'
idx_C<-grep('C',colnames(count_analysis_X))
coldata_X$type[idx_C]<-'BC'
idx_H<-grep('H',colnames(count_analysis_X))
coldata_X$type[idx_H]<-'HC'


##备用编组
#idx_R<-grep('C89',colnames(count_analysis_N))
#coldata$type[idx_R]<-'R'
#idx_R<-grep('C103',colnames(count_analysis_N))
#coldata$type[idx_R]<-'R'
#idx_R<-grep('C107',colnames(count_analysis_N))
#coldata$type[idx_R]<-'R'


dds <- DESeqDataSetFromMatrix(countData = round(count_analysis_X),
                              colData = coldata_X,
                              design = ~ type)
dds

smallestGroupSize <- 10
keep <- rowSums(counts(dds) > 5) >= smallestGroupSize
dds <- dds[keep,]
dds$type <- relevel(dds$type, ref = "HC")

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
colnames(sampleDistMatrix) <- colnames(count_analysis_X)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


plotPCA(vsd,ntop=1000, intgroup=c("type"))
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing=TRUE)[seq_len(min(1000, length(rv)))]

# 对选出的基因做 PCA
pca <- prcomp(t(assay(vsd)[select, ]))

# 提取 loadings (每个基因在每个PC上的权重)
loadings <- pca$rotation

# 找出最影响 PC1 的基因
pc1_load <- sort(abs(loadings[,1]), decreasing = TRUE)

# 查看前 20 个对 PC1 影响最大的基因
head(pc1_load, 20)

# 创建 PCA 数据
pcaData <- plotPCA(vsd, ntop=1000 ,intgroup = c("type"), returnData = TRUE)


# 自定义颜色
colors <- c("HC" = "#379FB4", "BS" = "#FCAE59", "BC" = "#DB498E")

# 绘制 PCA 图
ggplot(pcaData, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 3) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  stat_ellipse(aes(fill = type, color = type), geom = "polygon", alpha = 0.15) +  # 添加半透明圆圈
  theme_bw() +
  labs(title = "PCA Plot of Urine cfRNA",
       x = "PC1 76%",
       y = "PC2 3%") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank())
ggsave("X_PCA.png", width = 8, height = 8, dpi = 600)
ggsave("X_PCA.pdf", width = 8, height = 8, dpi = 600)

#BiocManager::install("umap")
library(umap)
# 提取归一化后的表达矩阵
mat <- assay(vsd)

# 运行 UMAP
set.seed(123)
topVarGenes <- head(order(rowVars(mat), decreasing=TRUE), 1000)
map_data <- umap(t(mat[topVarGenes, ]),
                 n_neighbors = 30,
                 min_dist = 0.3,
                 metric = "euclidean",
                 random_state = 123)
# 合并结果和样本信息
umap_df <- data.frame(map_data$layout, colData(vsd))
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")
umap_df$type <- coldata_N$type

# 画图
ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=type)) +
  geom_point(size=3, alpha=0.8) +
  theme_classic(base_size = 14) +
  labs(title="", color="Group") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("X_UMAP.png", width = 8, height = 8, dpi = 600)
ggsave("X_UMAP.pdf", width = 8, height = 8, dpi = 600)


nor_count_X <- assay(vst(dds))
write.csv(nor_count_X ,"nor_X_count.csv")

resultsNames(dds)
res_type_C_vs_H <- results(dds, name="type_BC_vs_HC")
hist(res_type_C_vs_H$pvalue, col = "grey50", border = "white", main = "", xlab = "p-values")
res_type_C_vs_H_res <- as.data.frame(res_type_C_vs_H)
#res_type_C_vs_H_res$ensembl_gene_id <- rownames(res_type_C_vs_H)
#res_type_C_vs_H_res <- merge(res_type_C_vs_H_res,dedu_symbols,by='ensembl_gene_id')
write.csv(res_type_C_vs_H_res,"X_type_C_vs_H_result.csv")

type_C_vs_H_result<- lfcShrink(dds, coef ="type_BC_vs_HC",type = 'apeglm',res=res_type_C_vs_H)
type_C_vs_H_result_res <- as.data.frame(type_C_vs_H_result)
DESeq2::plotMA(type_C_vs_H_result)
#type_C_vs_H_result_res$ensembl_gene_id <- rownames(type_C_vs_H_result)
#type_C_vs_H_result_res <- merge(type_C_vs_H_result_res,dedu_symbols,by='ensembl_gene_id')
write.csv(type_C_vs_H_result,"lfc_X_type_C_vs_H_result.csv")

res_type_S_vs_H <- results(dds, name="type_BS_vs_HC")
hist(res_type_S_vs_H$pvalue, col = "grey50", border = "white", main = "", xlab = "p-values")
res_type_S_vs_H_res <- as.data.frame(res_type_S_vs_H)
#res_type_S_vs_H_res$ensembl_gene_id <- rownames(res_type_S_vs_H)
#res_type_S_vs_H_res <- merge(res_type_S_vs_H_res,dedu_symbols,by='ensembl_gene_id')
write.csv(res_type_S_vs_H_res,"X_type_S_vs_H_result.csv")

type_S_vs_H_result<- lfcShrink(dds, coef ="type_BS_vs_HC",type = 'apeglm',res=res_type_S_vs_H)
type_S_vs_H_result_res <- as.data.frame(type_S_vs_H_result)
#type_S_vs_H_result_res$ensembl_gene_id <- rownames(type_S_vs_H_result)
#type_S_vs_H_result_res <- merge(type_S_vs_H_result_res,dedu_symbols,by='ensembl_gene_id')
write.csv(type_S_vs_H_result,"lfc_X_type_S_vs_H_result.csv")


dds$type <- relevel(dds$type, ref = "BS")

dds <- DESeq(dds)
resultsNames(dds)

res_type_C_vs_S <- results(dds, name="type_BC_vs_BS")
hist(res_type_C_vs_S$pvalue, col = "grey50", border = "white", main = "", xlab = "p-values")
res_type_C_vs_S_res <- as.data.frame(res_type_C_vs_S)
#res_type_C_vs_S_res$ensembl_gene_id <- rownames(res_type_C_vs_S)
#res_type_C_vs_S_res <- merge(res_type_C_vs_S_res,dedu_symbols,by='ensembl_gene_id')
write.csv(res_type_C_vs_S_res,"X_type_C_vs_S_result.csv")

type_C_vs_S_result<- lfcShrink(dds, coef ="type_BC_vs_BS",type = 'apeglm',res=res_type_C_vs_S)
type_C_vs_S_result_res <- as.data.frame(type_C_vs_S_result)
#type_C_vs_S_result_res$ensembl_gene_id <- rownames(type_C_vs_S_result)
#type_C_vs_S_result_res <- merge(type_C_vs_S_result_res,dedu_symbols,by='ensembl_gene_id')
write.csv(type_C_vs_S_result,"lfc_X_type_C_vs_S_result.csv")



dataForVolcanoPlot <- as.data.frame(type_S_vs_H_result_res)

logFcThreshold <- log2(2)
adjPvalThreshold <- 0.05
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange < logFcThreshold | padj> adjPvalThreshold)] <- 'NS'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange >= logFcThreshold & padj <= adjPvalThreshold)] <- 'UP'
dataForVolcanoPlot$Significance[with(dataForVolcanoPlot,log2FoldChange <= -logFcThreshold & padj <= adjPvalThreshold)] <- 'DOWN'


dataForVolcanoPlot
dataForVolcanoPlot$ID = rownames(dataForVolcanoPlot)

#cols <- c('UP'=="red", 'NS'='gray','DOWN'="green")
cols <- c('UP'="pink", 'NS'='gray','DOWN'="lightblue")


dev.new()

p <- ggplot(dataForVolcanoPlot, aes(x = log2FoldChange, y = -log10(padj))) +
  #xlim(-2,2) +
  labs(x=expression(bold('Log'['2']*'(Fold Change)')), 
       y=(expression(bold('-Log'['10']*'(FDR)'))), 
       title=NULL) +
  geom_point(color='gray', alpha=1, size=1) +
  geom_vline(xintercept = c(-logFcThreshold, logFcThreshold),
             color='darkgreen', linetype='dashed') +
  geom_hline(yintercept = -log10(adjPvalThreshold), 
             color='darkgreen',linetype='dashed')+
  scale_color_manual(values = cols) +
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'UP'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="pink", 
             size=1, color="pink",shape=21)+
  geom_point(data = subset(dataForVolcanoPlot,
                           Significance == 'DOWN'),
             aes(x=log2FoldChange, y=-log10(padj)), fill="lightblue", 
             size=1, color="lightblue",shape=21) +
  scale_fill_manual(values = cols) +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'UP' & (-log10(padj) > 4 & log2FoldChange > 1)),
                  aes(label = ID), nudge_x = 0.1,
                  segment.alpha = 0.5, #segment.size = 0.5,
                  min.segment.length = 0.2,
                  size = 3, color='black', segment.color = 'black') +
  geom_text_repel(data = subset(dataForVolcanoPlot,
                                Significance == 'DOWN' & log2FoldChange < -1 & -log10(padj) > 4),
                  aes(label = ID),
                  segment.alpha = 0.5, #segment.size = 0.5,
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

ggsave("X_S_H_volcano.png", width = 12, height = 8, dpi = 600)
ggsave("X_S_H_volcano.pdf", width = 12, height = 8, dpi = 600)


DEGresult<-as.data.frame(type_C_vs_H_result)
#logFCThreshold=1
#adj.P.ValThreshold=0.05
#result_up_gene<-subset(DEGresult,padj<adj.P.ValThreshold&log2FoldChange>logFCThreshold)
#result_down_gene<-subset(DEGresult,padj<adj.P.ValThreshold&log2FoldChange<(-logFCThreshold))
#result_DEG<-subset(DEGresult,padj<adj.P.ValThreshold&abs(log2FoldChange)>logFCThreshold)


result_DEG<-subset(DEGresult,padj<0.05&abs(log2FoldChange)>2)

select<-which(rownames(dds) %in% rownames(result_DEG))
df <- as.data.frame(coldata_N)


analysisresult<-as.data.frame(type_S_vs_H_result_res)
logFCThreshold=1
adj.P.ValThreshold=0.05
result_up_gene_list<-subset(analysisresult,padj<adj.P.ValThreshold&log2FoldChange>logFCThreshold)
result_down_gene_list<-subset(analysisresult,padj<adj.P.ValThreshold&log2FoldChange<(-logFCThreshold))
#result_DEG<-subset(analysisresult,padj<adj.P.ValThreshold&abs(log2FoldChange)>logFCThreshold)

result_up_gene<-rownames(result_up_gene_list)
result_down_gene<-rownames(result_down_gene_list)
#up
#goMF_up <- enrichGO(result_up_gene,
#                    keyType = "SYMBOL",    # 类型
#                    OrgDb = "org.Hs.eg.db",  # 
#                    ont = "MF",
#                    pAdjustMethod= 'BH',
#                    pvalueCutoff= 0.05,
#                    qvalueCutoff= 0.2)    # 分子功能
#goBP_up <- enrichGO(result_up_gene,
#                    keyType = "SYMBOL",
#                    OrgDb = "org.Hs.eg.db",
#                    ont = "BP",
#                    pAdjustMethod= 'BH',
#                    pvalueCutoff= 0.05,
#                    qvalueCutoff= 0.2)  # 生物学过程
#goCC_up <- enrichGO(result_up_gene,
#                    keyType = "SYMBOL",
#                    OrgDb = "org.Hs.eg.db",
#                    ont = "CC",
#                    pAdjustMethod= 'BH',
#                    pvalueCutoff= 0.05,
#                    qvalueCutoff= 0.2) # 细胞组分
goAll_up <- enrichGO(result_up_gene,
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "ALL",
                     pAdjustMethod= 'BH',
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.2)  # 以上三种  


#down
#goMF_down <- enrichGO(result_down_gene,
#                      keyType = "SYMBOL",    # degEntrezIdl类型
#                      OrgDb = "org.Hs.eg.db",  # 人的数据库
#                      ont = "MF",
#                      pAdjustMethod= 'BH',
#                      pvalueCutoff= 0.05,
#                      qvalueCutoff= 0.2)   # 分子功能
#goBP_down <- enrichGO(result_down_gene,
#                      keyType = "SYMBOL",
#                      OrgDb = "org.Hs.eg.db",
#                      ont = "BP",
#                      pAdjustMethod= 'BH',
#                      pvalueCutoff= 0.05,
#                      qvalueCutoff= 0.2)  # 生物学过程
#goCC_down <- enrichGO(result_down_gene,
#                      keyType = "SYMBOL",
#                      OrgDb = "org.Hs.eg.db",
#                      ont = "CC",
#                      pAdjustMethod= 'BH',
#                      pvalueCutoff= 0.05,
#                      qvalueCutoff= 0.2) # 细胞组分
goAll_down <- enrichGO(result_down_gene,
                       keyType = "SYMBOL",
                       OrgDb = "org.Hs.eg.db",
                       ont = "ALL",
                       pAdjustMethod= 'BH',
                       pvalueCutoff= 0.05,
                       qvalueCutoff= 0.2)  # 以上三种  

go_for_plot <- goAll_up
library(enrichplot)
p2 <- dotplot(go_for_plot, font.size = 10,
              title= "GO Dot plot",
              label_format= 50,
              showCategory=10)
p2

ggsave("X_S_H_up_GO.png", width = 6, height = 4, dpi = 600)
ggsave("X_S_H_up_GO.pdf", width = 6, height = 4, dpi = 600)

go_for_plot <- goAll_down
library(enrichplot)
p3 <- dotplot(go_for_plot, font.size = 10,
              title= "GO Dot plot",
              label_format= 50,
              showCategory=10)
p3

ggsave("X_S_H_down_GO.png", width = 6, height = 4, dpi = 600)
ggsave("X_S_H_down_GO.pdf", width = 6, height = 4, dpi = 600)



geneList_pre<-cbind(rownames(type_C_vs_H_result_res),type_C_vs_H_result_res$log2FoldChange)
geneList = as.numeric(geneList_pre[,2])
names(geneList) = as.character(geneList_pre[,1])
geneList = sort(geneList, decreasing = TRUE)
geneList<- geneList[!is.na(names(geneList))]

ego3 <- gseGO(geneList     = geneList,
              OrgDb        = "org.Hs.eg.db",
              keyType      = "SYMBOL",
              ont          = "ALL",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

ridgeplot(ego3,
          showCategory = 10,      # 显示前多少条通路
          fill = "p.adjust",      # 填充颜色，常用 p.adjust / NES
          label_format = 50       # 通路名称换行宽度
)
ggsave("X_S_H_gseGO.png", width = 8, height = 4, dpi = 600)
ggsave("X_S_H_gseGO.pdf", width = 8, height = 4, dpi = 600)





