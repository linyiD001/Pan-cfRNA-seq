###intersection between all groups
#BiocManager::install("ComplexHeatmap")
library(devtools)
library(org.Hs.eg.db)
library(stringr)
#devtools::install_github("hms-dbmi/UpSetR")
#install_github("jokergoo/ComplexHeatmap",force = TRUE)
#BiocManager::install("Matrix",force = TRUE)
setwd("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar")
library(UpSetR)
library(ComplexHeatmap)
N_C_H_DEG<-read.csv("N_type_C_vs_H_result.csv",header=T,row.names = 1)
X_C_H_DEG<-read.csv("X_type_C_vs_H_result.csv",header=T,row.names = 1)
N_S_H_DEG<-read.csv("N_type_S_vs_H_result.csv",header=T,row.names = 1)
X_S_H_DEG<-read.csv("X_type_S_vs_H_result.csv",header=T,row.names = 1)
N_C_S_DEG<-read.csv("N_type_C_vs_S_result.csv",header=T,row.names = 1)
X_C_S_DEG<-read.csv("X_type_C_vs_S_result.csv",header=T,row.names = 1)
#T_P_paired_DEG<-read.csv("/dssg/home/acct-dahan/share/BC/tissue/RESULT/NORMALIZATION/paired_type_T_vs_P_result.csv",header=T,row.names = 1)
T_P_unpaired_DEG<-read.csv("/dssg/home/acct-dahan/share/BC/tissue/RESULT/NORMALIZATION/unpaired_type_T_vs_P_result.csv",header=T,row.names = 1)

padjthreshold<-0.05
log2FoldChangethreshold<-2
basemeanthreshold<-10
tissuelog2FoldChangethreshold<-log2(2)
up_N_C_H_DEG<-rownames(N_C_H_DEG[which(N_C_H_DEG$padj<padjthreshold&N_C_H_DEG$log2FoldChange>log2FoldChangethreshold&N_C_H_DEG$baseMean>basemeanthreshold),])
up_N_C_S_DEG<-rownames(N_C_S_DEG[which(N_C_S_DEG$padj<padjthreshold&N_C_S_DEG$log2FoldChange>log2FoldChangethreshold&N_C_S_DEG$baseMean>basemeanthreshold),])
up_N_S_H_DEG<-rownames(N_S_H_DEG[which(N_S_H_DEG$padj<padjthreshold&N_S_H_DEG$log2FoldChange>log2FoldChangethreshold&N_S_H_DEG$baseMean>basemeanthreshold),])

up_X_C_H_DEG<-rownames(X_C_H_DEG[which(X_C_H_DEG$padj<padjthreshold&X_C_H_DEG$log2FoldChange>log2FoldChangethreshold&X_C_H_DEG$baseMean>basemeanthreshold),])
up_X_C_S_DEG<-rownames(X_C_S_DEG[which(X_C_S_DEG$padj<padjthreshold&X_C_S_DEG$log2FoldChange>log2FoldChangethreshold&X_C_S_DEG$baseMean>basemeanthreshold),])
up_X_S_H_DEG<-rownames(X_S_H_DEG[which(X_S_H_DEG$padj<padjthreshold&X_S_H_DEG$log2FoldChange>log2FoldChangethreshold&X_S_H_DEG$baseMean>basemeanthreshold),])

#up_T_P_paired_DEG<-rownames(T_P_paired_DEG[which(T_P_paired_DEG$padj<padjthreshold&T_P_paired_DEG$log2FoldChange>tissuelog2FoldChangethreshold&T_P_paired_DEG$baseMean>0),])
up_T_P_unpaired_DEG<-rownames(T_P_unpaired_DEG[which(T_P_unpaired_DEG$padj<padjthreshold&T_P_unpaired_DEG$log2FoldChange>tissuelog2FoldChangethreshold&T_P_unpaired_DEG$baseMean>0),])

down_N_C_H_DEG<-rownames(N_C_H_DEG[which(N_C_H_DEG$padj<padjthreshold&N_C_H_DEG$log2FoldChange<(-log2FoldChangethreshold)&N_C_H_DEG$baseMean>basemeanthreshold),])
down_N_C_S_DEG<-rownames(N_C_S_DEG[which(N_C_S_DEG$padj<padjthreshold&N_C_S_DEG$log2FoldChange<(-log2FoldChangethreshold)&N_C_S_DEG$baseMean>basemeanthreshold),])
down_N_S_H_DEG<-rownames(N_S_H_DEG[which(N_S_H_DEG$padj<padjthreshold&N_S_H_DEG$log2FoldChange<(-log2FoldChangethreshold)&N_S_H_DEG$baseMean>basemeanthreshold),])

down_X_C_H_DEG<-rownames(X_C_H_DEG[which(X_C_H_DEG$padj<padjthreshold&X_C_H_DEG$log2FoldChange<(-log2FoldChangethreshold)&X_C_H_DEG$baseMean>basemeanthreshold),])
down_X_C_S_DEG<-rownames(X_C_S_DEG[which(X_C_S_DEG$padj<padjthreshold&X_C_S_DEG$log2FoldChange<(-log2FoldChangethreshold)&X_C_S_DEG$baseMean>basemeanthreshold),])
down_X_S_H_DEG<-rownames(X_S_H_DEG[which(X_S_H_DEG$padj<padjthreshold&X_S_H_DEG$log2FoldChange<(-log2FoldChangethreshold)&X_S_H_DEG$baseMean>basemeanthreshold),])

#lt = list(Upregulated_in_BCvsHC_Urine_cfRNA = up_N_C_H_DEG,
#          Upregulated_in_BCvsBS_Urine_cfRNA = up_N_C_S_DEG,
#          Upregulated_in_BSvsHC_Urine_cfRNA = up_N_S_H_DEG,
#          Upregulated_in_BCvsHC_Plasma_cfRNA = up_X_C_H_DEG,
#          Upregulated_in_BCvsBS_Plasma_cfRNA = up_X_C_S_DEG,
#          Upregulated_in_BSvsHC_Plasma_cfRNA = up_X_S_H_DEG,
#          Upregulated_in_Tumor_tissue = up_T_P_unpaired_DEG)
#          #Urine_down_BC_HC = down_N_C_H_DEG,
#          #Urine_down_BC_BS = down_N_C_S_DEG,
#          #Plasma_down_BC_HC = down_X_C_H_DEG,
#          #Plasma_down_BC_BS = down_X_C_S_DEG,
#          Tissue_down = down_T_P_unpaired_DEG)

lt = list(Up_in_BC_Urine_cfRNA = up_N_C_H_DEG,
          Up_in_BC_Plasma_cfRNA = up_X_C_H_DEG,
          Up_in_tumor = up_T_P_unpaired_DEG
          #Up_in_NAT_tissue= up_P_N_DEG
          #Only_Up_in_BC_Urine_cfRNA = up_N_C_S_DEG,
          #Only_Up_in_BC_Plasma_cfRNA = up_X_C_S_DEG,
          #Urine_down_BC_HC = down_N_C_H_DEG,
          #Plasma_down_BC_HC = down_X_C_H_DEG,
          #Only_Up_in_Tumor_tissue= up_T_P_DEG
          )

lt <- lapply(lt, unique)

# Create UpSet data

m = make_comb_mat(lt,mode = "intersect")

comb_name(m)

m_filtered = m[comb_degree(m) >= 1]
#m_filtered <- m[comb_degree(m) >= 2 & comb_name(m) != "1100"]
# 绘制 UpSet 图并注释每个类别的数量

library(Cairo)
png("UpSet_plot.png",width=8,height=5,res=300, units = "in")  # 设置分辨率
#CairoPDF("UpSet_plot.pdf",width=8,height=5,family='Arial')  # 设置分辨率

UpSet(
  m_filtered,
  comb_col = "darkblue",
  set_order = c( "Up_in_BC_Urine_cfRNA", 
                 "Up_in_BC_Plasma_cfRNA",
                 #"Only_Up_in_BC_Urine_cfRNA",
                 #"Only_Up_in_BC_Plasma_cfRNA",
                 "Up_in_tumor"
                 #"Up_in_NAT_tissue"
                 #"Only_Up_in_Tumor_tissue"
                 ),
  top_annotation = HeatmapAnnotation(
    "Intersection Size" = anno_barplot(
      comb_size(m_filtered),
      height = unit(10, "cm"),
      add_numbers = TRUE,
      numbers_gp = gpar(fontsize = 10, col = "black")
    ),
    annotation_name_gp = gpar(fontsize = 12),       # 控制字体大小
    annotation_name_rot = 90,                       # 设置标题文字为竖着显示
    gp = gpar(fill = "steelblue")),right_annotation = NULL )
dev.off()

genelist<-Reduce(intersect,list(up_N_C_H_DEG,up_N_C_S_DEG,up_N_S_H_DEG))#,up_T_P_unpaired_DEG))
genelist
#write.csv(genelist,"N_genelist_up.csv")

genelist<-Reduce(intersect,list(down_N_C_H_DEG,down_N_C_S_DEG,down_N_S_H_DEG))#,down_T_P_unpaired_DEG))
genelist
#write.csv(genelist,"N_genelist_down.csv")

intersection_gene_N_tumor<-Reduce(intersect,list(up_N_C_H_DEG,up_T_P_unpaired_DEG))
#write.csv(intersection_gene_N_up,"N_tissue_intersection_up_gene.csv")

intersection_gene_X_tumor<-Reduce(intersect,list(up_X_C_H_DEG,up_T_P_unpaired_DEG))
#write.csv(intersection_gene_X_up,"X_tissue_intersection_up_gene.csv")

######N_up_Tumor_up
select_N_CH_gene<-N_C_H_DEG[which(rownames(N_C_H_DEG)%in%intersection_gene_N_tumor),]
select_T_N_gene<-T_P_unpaired_DEG[which(rownames(T_P_unpaired_DEG)%in%intersection_gene_N_tumor),]

select_N_CH_gene$ID<-rownames(select_N_CH_gene)
select_T_N_gene$ID<-rownames(select_T_N_gene)

dataplot<-merge(select_N_CH_gene,select_T_N_gene,by='ID')
library(ggplot2)

result <- cor.test(dataplot$log2FoldChange.x,dataplot$log2FoldChange.y, method = "pearson")
result$estimate  # 输出相关系数
result$p.value   # 输出P值
p_value <- result$p.value  # 提取P值
ymax <- max(dataplot$log2FoldChange.y, na.rm = TRUE)
ymin <- min(dataplot$log2FoldChange.x, na.rm = TRUE)

xmax <- max(dataplot$log2FoldChange.x, na.rm = TRUE)
xmin <- min(dataplot$log2FoldChange.x, na.rm = TRUE)

p <- ggplot(dataplot, aes(x = log2FoldChange.x, y = log2FoldChange.y)) +
  geom_point(color = "black", alpha = 0.6, size = 2) +
  
  # 拟合线（线性回归）
  geom_smooth(method = "lm", se = FALSE, color = "#E64B35FF",
              linewidth = 1, linetype = "solid") +
  
  # 每个点加上注释（如基因名或ID）
  geom_text_repel(
    aes(label = dataplot$ID),  # 假设行名是基因名
    size = 2.5, 
    family = "Arial",
    color = "black",
    max.overlaps = 50,  # 控制重叠数量（可调大或调小）
    box.padding = 0.3,
    point.padding = 0.2
  ) +
  
  # 标注 P 值
  geom_text(
    aes(
      x = (max(log2FoldChange.x) - min(log2FoldChange.x)) * 0.8 + min(log2FoldChange.x),
      y = (max(log2FoldChange.y) - min(log2FoldChange.y)) * 0.9 + min(log2FoldChange.y)
    ),
    label = paste("P =", format(result$p.value, digits = 3)),
    color = "black", size = 3.5, family = "Arial"
  ) +
  
  # 标注相关系数
  geom_text(
    aes(
      x = (max(log2FoldChange.x) - min(log2FoldChange.x)) * 0.8 + min(log2FoldChange.x),
      y = (max(log2FoldChange.y) - min(log2FoldChange.y)) * 0.8 + min(log2FoldChange.y)
    ),
    label = paste("Correlation =", format(result$estimate, digits = 3)),
    color = "black", size = 3.5, family = "Arial"
  ) +
  
  labs(
    x = "Log2FC in Urine cfRNA",
    y = "Log2FC in Tumor Tissue"
  ) +
  
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

p
#-------------------------------------------
# 使用 Cairo 保存
#-------------------------------------------
CairoPDF("N_Tumor_DEG_up_cor.pdf",height=50,width=50,dpi=1200)
print(p)
dev.off()

CairoPNG("N_Tumor_DEG_up_cor.png",height=4,width=4,dpi=1200,units = "in",res = 300)
print(p)
dev.off()

library(clusterProfiler)

goAll_up <- enrichGO(intersection_gene_N_tumor,
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod= 'BH',
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.2)  # 以上三种  

top10_go <- goAll_up@result %>%
  dplyr::arrange(p.adjust) %>%      
  dplyr::slice_head(n = 5) %>%      
  dplyr::mutate(GeneRatio = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))) %>%  # 转为数值
  dplyr::arrange(desc(GeneRatio))  
go_for_plot <- top10_go
go_for_plot$Description <- str_wrap(go_for_plot$Description, width = 30)
p2_ggplot <- ggplot(go_for_plot, aes(x = reorder(Description, GeneRatio), y = GeneRatio, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.7) +  # 使用点绘制
  coord_flip() +  # 翻转坐标轴
  scale_size_continuous(range = c(3, 10)) +  # 设置点的大小范围
  scale_color_gradient(low = "blue", high = "red") +  # p-value 渐变色
  theme_minimal() +
  labs(x = "Gene Ontology enrichment analysis", y = "Gene Ratio", size = "Count", color = "Adjusted p-value") +
  theme(
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size = 20, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   # 设置图例项的水平间距
    legend.spacing.y = unit(1, 'cm'),   # 设置图例项的垂直间距
    legend.key.spacing.y = unit(0.1, "cm"),    # 设置图例色块的高度
    legend.key.height = unit(0.5, "cm"),    # 设置图例色块的高度
    legend.key.width = unit(0.5, "cm"),     # 设置图例色块的宽度
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    # 图例背景和边框
    #legend.background = element_rect(fill = "white", color = "black", size = 0.5), # 图例背景和边框
    #legend.key = element_rect(color = "black", size = 0.5)  # 每个图例项的边框
  ) +
  theme(legend.position = "right")


# 显示图形
print(p2_ggplot)

ggsave("N_Tumor_DEG_overlapping_up_KEGG.png", width = 8, height = 5,dpi=600,units = "in")
ggsave("N_Tumor_DEG_overlapping_up_KEGG.pdf", width = 8, height = 5,dpi=600,units = "in", device = cairo_pdf)

