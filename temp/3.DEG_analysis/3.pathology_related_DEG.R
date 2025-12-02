###intersection 
library(devtools)
library(org.Hs.eg.db)
library(Cairo)

repeataware_genelist<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome/all_genelist_annotation.csv",header=T,row.names = 1)

non_re_idx<-grep('ENST',repeataware_genelist$enst)
noRE_repeataware_genelist<-repeataware_genelist[non_re_idx,]
genename<-noRE_repeataware_genelist[!duplicated(noRE_repeataware_genelist[,4]),4]

library(UpSetR)
library(ComplexHeatmap)
N_C_H_DEG<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_repeataware/salmon_repeataware_stranded_allgenome/N_type_C_vs_H_result.csv",header=T,row.names = 1)
X_C_H_DEG<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_repeataware/salmon_repeataware_stranded_allgenome/X_type_C_vs_H_result.csv",header=T,row.names = 1)
N_S_H_DEG<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_repeataware/salmon_repeataware_stranded_allgenome/N_type_S_vs_H_result.csv",header=T,row.names = 1)
X_S_H_DEG<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_repeataware/salmon_repeataware_stranded_allgenome/X_type_S_vs_H_result.csv",header=T,row.names = 1)
N_C_S_DEG<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_repeataware/salmon_repeataware_stranded_allgenome/N_type_C_vs_S_result.csv",header=T,row.names = 1)
X_C_S_DEG<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_repeataware/salmon_repeataware_stranded_allgenome/X_type_C_vs_S_result.csv",header=T,row.names = 1)
T_P_paired_DEG<-read.csv("/dssg/home/acct-dahan/share/tissue/RESULT/NORMALIZATION/paired_type_T_vs_P_result.csv",header=T,row.names = 1)
T_P_unpaired_DEG<-read.csv("/dssg/home/acct-dahan/share/tissue/RESULT/NORMALIZATION/unpaired_type_T_vs_P_result.csv",header=T,row.names = 1)

idx_noname<-which(T_P_unpaired_DEG$hgnc_symbol=='')
T_P_unpaired_DEG[idx_noname,7]<-T_P_unpaired_DEG[idx_noname,1]


idx_noname<-which(T_P_paired_DEG$hgnc_symbol=='')
T_P_paired_DEG[idx_noname,7]<-T_P_paired_DEG[idx_noname,1]


padjthreshold<-0.05
log2FoldChangethreshold<-1
basemeanthreshold<-0
tissuelog2FoldChangethreshold<-log2(2)
up_N_C_H_DEG<-rownames(N_C_H_DEG[which(N_C_H_DEG$padj<padjthreshold&N_C_H_DEG$log2FoldChange>log2FoldChangethreshold&N_C_H_DEG$baseMean>basemeanthreshold&rownames(N_C_H_DEG)%in%genename),])
up_N_C_S_DEG<-rownames(N_C_S_DEG[which(N_C_S_DEG$padj<padjthreshold&N_C_S_DEG$log2FoldChange>log2FoldChangethreshold&N_C_S_DEG$baseMean>basemeanthreshold&rownames(N_C_S_DEG)%in%genename),])
up_X_C_H_DEG<-rownames(X_C_H_DEG[which(X_C_H_DEG$padj<padjthreshold&X_C_H_DEG$log2FoldChange>log2FoldChangethreshold&X_C_H_DEG$baseMean>basemeanthreshold&rownames(X_C_H_DEG)%in%genename),])
up_X_C_S_DEG<-rownames(X_C_S_DEG[which(X_C_S_DEG$padj<padjthreshold&X_C_S_DEG$log2FoldChange>log2FoldChangethreshold&X_C_S_DEG$baseMean>basemeanthreshold&rownames(X_C_S_DEG)%in%genename),])
up_T_P_unpaired_DEG<-T_P_unpaired_DEG[which(T_P_unpaired_DEG$padj<padjthreshold&T_P_unpaired_DEG$log2FoldChange>tissuelog2FoldChangethreshold&T_P_unpaired_DEG$baseMean>0&T_P_unpaired_DEG$hgnc_symbol!= ''),7]
down_N_C_H_DEG<-rownames(N_C_H_DEG[which(N_C_H_DEG$padj<padjthreshold&N_C_H_DEG$log2FoldChange<(-log2FoldChangethreshold)&N_C_H_DEG$baseMean>basemeanthreshold&rownames(N_C_H_DEG)%in%genename),])
down_N_C_S_DEG<-rownames(N_C_S_DEG[which(N_C_S_DEG$padj<padjthreshold&N_C_S_DEG$log2FoldChange<(-log2FoldChangethreshold)&N_C_S_DEG$baseMean>basemeanthreshold&rownames(N_C_S_DEG)%in%genename),])
down_X_C_H_DEG<-rownames(X_C_H_DEG[which(X_C_H_DEG$padj<padjthreshold&X_C_H_DEG$log2FoldChange<(-log2FoldChangethreshold)&X_C_H_DEG$baseMean>basemeanthreshold&rownames(X_C_H_DEG)%in%genename),])
down_X_C_S_DEG<-rownames(X_C_S_DEG[which(X_C_S_DEG$padj<padjthreshold&X_C_S_DEG$log2FoldChange<(-log2FoldChangethreshold)&X_C_S_DEG$baseMean>basemeanthreshold&rownames(X_C_S_DEG)%in%genename),])
down_T_P_unpaired_DEG<-T_P_unpaired_DEG[which(T_P_unpaired_DEG$padj<padjthreshold&T_P_unpaired_DEG$log2FoldChange<(-tissuelog2FoldChangethreshold)&T_P_unpaired_DEG$baseMean>0&T_P_unpaired_DEG$hgnc_symbol!= ''),7]


lt = list(Upregulated_in_Urine_cfRNA = up_N_C_H_DEG,
          Upregulated_in_Plasma_cfRNA = up_X_C_H_DEG,
          Upregulated_in_Tumor_tissue= up_T_P_unpaired_DEG,
          #Urine_down_BC_HC = down_N_C_H_DEG,
          #Plasma_down_BC_HC = down_X_C_H_DEG,
          Upregulated_in_NAT_tissue= down_T_P_unpaired_DEG)


lt <- lapply(lt, unique)

# Create UpSet data

m = make_comb_mat(lt,mode = "intersect")
m_filtered = m[comb_degree(m) >= 2]

png("UpSet_plot.png",width=8,height=5,res=300, units = "in")  
CairoPDF("UpSet_plot.pdf",width=8,height=5,family='Arial')  

UpSet(
  m_filtered,
  comb_col = "darkblue",
  set_order = c( "Upregulated_in_Urine_cfRNA", "Upregulated_in_Plasma_cfRNA","Upregulated_in_Tumor_tissue", "Upregulated_in_NAT_tissue"),
  top_annotation = HeatmapAnnotation(
    "Intersection Size" = anno_barplot(
      comb_size(m_filtered),
      height = unit(10, "cm"),
      add_numbers = TRUE,
      numbers_gp = gpar(fontsize = 10, col = "black")
    ),
    annotation_name_gp = gpar(fontsize = 12),       
    annotation_name_rot = 90,                       
    gp = gpar(fill = "steelblue")),right_annotation = NULL )
dev.off()

intersection_gene_N_up<-Reduce(intersect,list(up_N_C_H_DEG,up_T_P_unpaired_DEG))
#write.csv(intersection_gene_N_up,"N_tissue_intersection_up_gene.csv")

intersection_gene_N_down<-Reduce(intersect,list(up_N_C_H_DEG,down_T_P_unpaired_DEG))
#write.csv(intersection_gene_N_down,"N_tissue_intersection_down_gene.csv")

intersection_gene_X_up<-Reduce(intersect,list(up_X_C_H_DEG,up_T_P_unpaired_DEG))
#write.csv(intersection_gene_X_up,"X_tissue_intersection_up_gene.csv")

intersection_gene_X_down<-Reduce(intersect,list(up_X_C_H_DEG,down_T_P_unpaired_DEG))
#write.csv(intersection_gene_X_down,"X_tissue_intersection_down_gene.csv")

######N_up_tissue_up
select_N_CH_gene<-N_C_H_DEG[which(rownames(N_C_H_DEG)%in%intersection_gene_N_up),]
select_T_P_gene<-T_P_unpaired_DEG[which(T_P_unpaired_DEG$hgnc_symbol%in%intersection_gene_N_up),]

select_N_CH_gene$ID<-rownames(select_N_CH_gene)
select_T_P_gene$ID<-select_T_P_gene$hgnc_symbol

dataplot<-merge(select_N_CH_gene,select_T_P_gene,by='ID')
library(ggplot2)

result <- cor.test(dataplot$log2FoldChange.x,dataplot$log2FoldChange.y, method = "pearson")
result$estimate  
result$p.value   
p_value <- result$p.value  
ymax <- max(dataplot$log2FoldChange.y, na.rm = TRUE)
ymin <- min(dataplot$log2FoldChange.x, na.rm = TRUE)

xmax <- max(dataplot$log2FoldChange.x, na.rm = TRUE)
xmin <- min(dataplot$log2FoldChange.x, na.rm = TRUE)

ggplot(dataplot, aes(x = log2FoldChange.x, y = log2FoldChange.y)) +
  geom_point(color = "black", alpha = 0.6, size = 2) +
  geom_text(aes(x = (xmax-xmin)*0.8+xmin, y =(ymax-ymin)*0.9+xmin),
            label = paste("P =", format(result$p.value, digits = 3)),
            color = "black", size = 3.5,family = "sans") +
  geom_text(aes(x = (xmax-xmin)*0.8+xmin,y =(ymax-ymin)*0.8+xmin, label = paste("correlation =", format(result$estimate, digits = 3))),
            color = "black", size = 3.5,family = "sans") +  
  labs(title = "", x = "Log2FC_in_Urine_CfRNA", y = "Log2FC_in_tissue") +
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
ggsave("N_Tissue_DEG_up_cor.png",height=4,width=4,dpi=1200,units = "in")
ggsave("N_Tissue_DEG_up_cor.pdf",height=4,width=4,dpi=1200,units = "in", device = cairo_pdf)

library(clusterProfiler)

goAll_up <- enrichGO(intersection_gene_N_up,
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod= 'BH',
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.2)  

top10_go <- goAll_up@result %>%
  dplyr::arrange(p.adjust) %>%      
  dplyr::slice_head(n = 8) %>%      
  dplyr::mutate(GeneRatio = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))) %>%  # 转为数值
  dplyr::arrange(desc(GeneRatio))  
go_for_plot <- top10_go

p2_ggplot <- ggplot(go_for_plot, aes(x = reorder(Description, GeneRatio), y = GeneRatio, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.7) +  
  coord_flip() +  
  scale_size_continuous(range = c(3, 10)) +  
  scale_color_gradient(low = "blue", high = "red") +  
  theme_minimal() +
  labs(title = "Both Upregulated in Urine and Tumor", x = "Gene Ontology enrichment analysis", y = "Gene Ratio", size = "Count", color = "Adjusted p-value") +
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
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    #legend.background = element_rect(fill = "white", color = "black", size = 0.5), 
    #legend.key = element_rect(color = "black", size = 0.5)  
  ) +
  theme(legend.position = "right")

print(p2_ggplot)

ggsave("N_Tissue_DEG_overlapping_up_KEGG.png", width = 7, height = 5,dpi=600,units = "in")
ggsave("N_Tissue_DEG_overlapping_up_KEGG.pdf", width = 7, height = 5,dpi=600,units = "in", device = cairo_pdf)


go_for_plot_graph<-as.data.frame(go_for_plot)
go_for_plot_graph_selected <- go_for_plot_graph[order(go_for_plot_graph$p.adjust, decreasing = F), ]
go_for_plot_graph_selected <-go_for_plot_graph_selected[1:5,]
colnames(go_for_plot_graph_selected)[2]<-'pathway'

gene_list_for_graph<-str_split(go_for_plot_graph_selected$geneID,"\\/")
pathway_list_for_graph<-go_for_plot_graph_selected$pathway

edge_list <- tibble(
  pathway = rep(pathway_list_for_graph, sapply(gene_list_for_graph, length)),
  gene = unlist(gene_list_for_graph)
)

logFC_vector<-N_C_H_DEG[which(rownames(N_C_H_DEG) %in% unique(unlist(gene_list_for_graph))),]
logFC_vector$gene <- rownames(logFC_vector)

edge_list <- left_join(edge_list, go_for_plot_graph_selected, by = "pathway")
edge_list <- left_join(edge_list, logFC_vector, by = "gene")

logFC_vector_2<-T_P_unpaired_DEG[which(T_P_unpaired_DEG$hgnc_symbol %in% unique(unlist(gene_list_for_graph))),-c(1:2,4:6)]
logFC_vector_2$gene <- logFC_vector_2$hgnc_symbol

edge_list <- left_join(edge_list, logFC_vector_2, by = "gene")


sector_colors <- c("mitotic nuclear division" = "#F5F5F5", 
                   "nuclear division" = "#FFB6C1", 
                   "sister chromatid segregation" = "#B0E0E6", 
                   "mitotic sister chromatid segregation" = "#FFFACD", 
                   "chromosome segregation" = "#E6E6FA")
circos.clear()

edge_list_plot<-edge_list[,c(1,2)]

logfc_col_fun <- colorRamp2(
  breaks = c(-8,0, 8),                      
  colors = c("blue","white", "red")        
)

sector_colors_grid <- setNames(
  logfc_col_fun(edge_list$log2FoldChange.x),
  edge_list$gene
)

sector_colors_combined <- c(sector_colors_grid, sector_colors)


#CairoPNG("N_Tissue_DEG_up_chord.png",width = 7, height = 5,res=600,units = "in")  
#CairoPDF("N_Tissue_DEG_up_chord.pdf",width = 7, height = 5,family = "Arial")  
link_colors <- sector_colors[edge_list_plot$pathway]

chordDiagram(edge_list_plot, 
             col = link_colors,  
             annotationTrack = c("grid"),  
             transparency = 0.1,  
             grid.col = sector_colors_combined
)
                                
circos.trackPlotRegion(
  track.index = 1,
  bg.border = NA,
  panel.fun = function(x, y) {
    sector_name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    circos.text(x = mean(xlim),
                y = ylim[2]+0.5,                     
                labels = sector_name,
                facing = "bending.inside",
                niceFacing = F,
                col = 'black',  
                cex = 0.4,
                family = "Arial")
  }
)

col_fun <- colorRamp2(c(-8,0,8), c("blue","white", "red"))
sector_colors <- setNames(col_fun(edge_list$log2FoldChange.y), edge_list_plot$gene)

circos.track(
  track.index = 2,  
  ylim = c(0, 1),
  bg.border = NA,
  track.height = 0.08,
  panel.fun = function(x, y) {
    sector_name <- get.cell.meta.data("sector.index")
    if (!(sector_name %in% names(sector_colors))) {
      return(NULL)
    }
    color <- sector_colors[sector_name]
    circos.rect(xleft = CELL_META$xlim[1],
                xright = CELL_META$xlim[2],
                ybottom = 0,
                ytop = 1,
                col = color,
                border = NA)
  }
)
lgd_logfc <- Legend(col_fun = logfc_col_fun,
                    title = "logFC",
                    at = c(-8, 0, 8),
                    title_gp = gpar(fontsize = 10, fontface = "plain", col = "black"), 
                    labels_gp = gpar(fontsize = 10, fontface = "plain", col = "black"),
                    labels = c("-8", "0", "8"))
draw(lgd_logfc, 
     x = unit(1, "npc") - unit(10, "mm"), 
     y = unit(0.5, "npc")) 
dev.off()


######N_up_tissue_down
select_N_CH_gene<-N_C_H_DEG[which(rownames(N_C_H_DEG)%in%intersection_gene_N_down),]
select_T_P_gene<-T_P_unpaired_DEG[which(T_P_unpaired_DEG$hgnc_symbol%in%intersection_gene_N_down),]

select_N_CH_gene$ID<-rownames(select_N_CH_gene)
select_T_P_gene$ID<-select_T_P_gene$hgnc_symbol

dataplot<-merge(select_N_CH_gene,select_T_P_gene,by='ID')
library(ggplot2)

result <- cor.test(dataplot$log2FoldChange.x,dataplot$log2FoldChange.y, method = "pearson")
result$estimate  
result$p.value   
p_value <- result$p.value  

ymax <- max(dataplot$log2FoldChange.y, na.rm = TRUE)
ymin <- min(dataplot$log2FoldChange.x, na.rm = TRUE)

xmax <- max(dataplot$log2FoldChange.x, na.rm = TRUE)
xmin <- min(dataplot$log2FoldChange.x, na.rm = TRUE)
ggplot(dataplot, aes(x = log2FoldChange.x, y = log2FoldChange.y)) +
  geom_point(color = "black", alpha = 0.6, size = 2) +
  geom_text(aes(x = (xmax-xmin)*0.8+xmin, y =(ymax-ymin)*0.9+xmin),
            label = paste("P =", format(result$p.value, digits = 3)),
            color = "black", size = 3.5,family = "sans") +
  geom_text(aes(x = (xmax-xmin)*0.8+xmin,y =(ymax-ymin)*0.8+xmin, label = paste("correlation =", format(result$estimate, digits = 3))),
            color = "black", size = 3.5,family = "sans") +  
  labs(title = "", x = "Log2FC_in_Urine_CfRNA", y = "Log2FC_in_tissue") +
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
#ggsave("N_Tissue_DEG_down_cor.png", height=4,width=4,dpi=1200,units = "in")
#ggsave("N_Tissue_DEG_down_cor.pdf", height=4,width=4,dpi=1200,units = "in",device = cairo_pdf)


goAll_down <- enrichGO(intersection_gene_X_down,
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod= 'BH',
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.2)    
top10_go <- goAll_down@result %>%
  dplyr::arrange(p.adjust) %>%      
  dplyr::slice_head(n = 8) %>%      
  dplyr::mutate(GeneRatio = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))) %>%  # 转为数值
  dplyr::arrange(desc(GeneRatio))  
go_for_plot <- top10_go

p2_ggplot <- ggplot(go_for_plot, aes(x = reorder(Description, GeneRatio), y = GeneRatio, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.7) +  
  coord_flip() +  
  scale_size_continuous(range = c(3, 10)) +  
  scale_color_gradient(low = "blue", high = "red") +  
  theme_minimal() +
  labs(title = "Both Upregulated in Urine and Tumor", x = "Gene Ontology enrichment analysis", y = "Gene Ratio", size = "Count", color = "Adjusted p-value") +
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
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    #legend.background = element_rect(fill = "white", color = "black", size = 0.5), 
    #legend.key = element_rect(color = "black", size = 0.5)  
  ) +
  theme(legend.position = "right")
print(p2_ggplot)
ggsave("N_Tissue_DEG_overlapping_down_KEGG.png", width = 7, height = 5,dpi=600,units = "in")
ggsave("N_Tissue_DEG_overlapping_down_KEGG.pdf", width = 7, height = 5,dpi=600,units = "in", device = cairo_pdf)


go_for_plot<-goAll_down
go_for_plot_graph<-as.data.frame(go_for_plot)
go_for_plot_graph_selected <- go_for_plot_graph[order(go_for_plot_graph$p.adjust, decreasing = F), ]
go_for_plot_graph_selected <-go_for_plot_graph_selected[1:5,]
colnames(go_for_plot_graph_selected)[2]<-'pathway'

gene_list_for_graph<-str_split(go_for_plot_graph_selected$geneID,"\\/")
pathway_list_for_graph<-go_for_plot_graph_selected$pathway

edge_list <- tibble(
  pathway = rep(pathway_list_for_graph, sapply(gene_list_for_graph, length)),
  gene = unlist(gene_list_for_graph)
)

logFC_vector<-N_C_H_DEG[which(rownames(N_C_H_DEG) %in% unique(unlist(gene_list_for_graph))),]
logFC_vector$gene <- rownames(logFC_vector)

edge_list <- left_join(edge_list, go_for_plot_graph_selected, by = "pathway")
edge_list <- left_join(edge_list, logFC_vector, by = "gene")

logFC_vector_2<-T_P_unpaired_DEG[which(T_P_unpaired_DEG$hgnc_symbol %in% unique(unlist(gene_list_for_graph))),-c(1:2,4:6)]
logFC_vector_2$gene <- logFC_vector_2$hgnc_symbol

edge_list <- left_join(edge_list, logFC_vector_2, by = "gene")



sector_colors <- c("adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains" = "#F5F5F5", 
                   "regulation of T cell activation" = "#FFB6C1", 
                   "positive regulation of leukocyte activation" = "#B0E0E6", 
                   "leukocyte mediated immunity" = "#FFFACD", 
                   "positive regulation of lymphocyte activation" = "#E6E6FA")
circos.clear()

edge_list_plot<-edge_list[,c(1,2)]

logfc_col_fun <- colorRamp2(
  breaks = c(-8,0, 8),                      
  colors = c("blue","white", "red")        
)

sector_colors_grid <- setNames(
  logfc_col_fun(edge_list$log2FoldChange.x),
  edge_list$gene
)

sector_colors_combined <- c(sector_colors_grid, sector_colors)

#CairoPNG("N_Tissue_DEG_down_chord.png",width = 7, height = 5,res=600,units = "in")  

CairoPDF("N_Tissue_DEG_down_chord.pdf",width = 7, height = 5,family = "Arial")  
link_colors <- sector_colors[edge_list_plot$pathway]

chordDiagram(edge_list_plot, 
             col = link_colors,  
             annotationTrack = c("grid"),  
             transparency = 0.1,  
             grid.col = sector_colors_combined
)

circos.trackPlotRegion(
  track.index = 1,
  bg.border = NA,
  panel.fun = function(x, y) {
    sector_name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    
    label_lines <- strwrap(sector_name, width = 40)
    n_lines <- length(label_lines)

    line_spacing <- 0.4
 
    total_height <- (n_lines - 1) * line_spacing
    start_y <- ylim[2] + 0.5 + total_height / 2 
    
    for (i in seq_along(label_lines)) {
      circos.text(
        x = mean(xlim),
        y = start_y - (i - 1) * line_spacing,
        labels = label_lines[i],
        facing = "bending.inside",
        niceFacing = FALSE,
        col = "black",
        cex = 0.3,
        family = "Arial"
      )
    }
  }
)

col_fun <- colorRamp2(c(-8,0,8), c("blue","white", "red"))
sector_colors <- setNames(col_fun(edge_list$log2FoldChange.y), edge_list_plot$gene)

circos.track(
  track.index = 2, 
  ylim = c(0, 1),
  bg.border = NA,
  track.height = 0.08,
  panel.fun = function(x, y) {
    sector_name <- get.cell.meta.data("sector.index")
    if (!(sector_name %in% names(sector_colors))) {
      return(NULL)
    }
    color <- sector_colors[sector_name]
    circos.rect(xleft = CELL_META$xlim[1],
                xright = CELL_META$xlim[2],
                ybottom = 0,
                ytop = 1,
                col = color,
                border = NA)
  }
)
lgd_logfc <- Legend(col_fun = logfc_col_fun,
                    title = "logFC",
                    at = c(-8, 0, 8),
                    title_gp = gpar(fontsize = 10, fontface = "plain", col = "black"), 
                    labels_gp = gpar(fontsize = 10, fontface = "plain", col = "black"),
                    labels = c("-8", "0", "8"))
draw(lgd_logfc, 
     x = unit(1, "npc") - unit(10, "mm"), 
     y = unit(0.5, "npc")) 
dev.off()


######X_up_tissue_up
select_X_CH_gene<-X_C_H_DEG[which(rownames(X_C_H_DEG)%in%intersection_gene_X_up),]
select_T_P_gene<-T_P_unpaired_DEG[which(T_P_unpaired_DEG$hgnc_symbol%in%intersection_gene_X_up),]

select_X_CH_gene$ID<-rownames(select_X_CH_gene)
select_T_P_gene$ID<-select_T_P_gene$hgnc_symbol

dataplot<-merge(select_X_CH_gene,select_T_P_gene,by='ID')
library(ggplot2)

result <- cor.test(dataplot$log2FoldChange.x,dataplot$log2FoldChange.y, method = "pearson")
result$estimate 
result$p.value  
p_value <- result$p.value  

ymax <- max(dataplot$log2FoldChange.y, na.rm = TRUE)
ymin <- min(dataplot$log2FoldChange.x, na.rm = TRUE)

xmax <- max(dataplot$log2FoldChange.x, na.rm = TRUE)
xmin <- min(dataplot$log2FoldChange.x, na.rm = TRUE)

ggplot(dataplot, aes(x = log2FoldChange.x, y = log2FoldChange.y)) +
  geom_point(color = "black", alpha = 0.6, size = 2) +
  geom_text(aes(x = (xmax-xmin)*0.8+xmin, y =(ymax-ymin)*0.9+xmin),
            label = paste("P =", format(result$p.value, digits = 3)),
            color = "black", size = 3.5,family = "sans") +
  geom_text(aes(x = (xmax-xmin)*0.8+xmin,y =(ymax-ymin)*0.8+xmin, label = paste("correlation =", format(result$estimate, digits = 3))),
            color = "black", size = 3.5,family = "sans") +  
  labs(title = "", x = "Log2FC_in_Plasma_CfRNA", y = "Log2FC_in_tissue") +
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


ggsave("X_Tissue_DEG_up_cor.png",height=4,width=4,dpi=1200,units = "in")
ggsave("X_Tissue_DEG_up_cor.pdf",height=4,width=4,dpi=1200,units = "in", device = cairo_pdf)


goAll_up <- enrichGO(intersection_gene_X_up,
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     ont = "ALL",
                     pAdjustMethod= 'BH',
                     pvalueCutoff= 0.05,
                     qvalueCutoff= 0.2)  

top10_go <- goAll_up@result %>%
  dplyr::arrange(p.adjust) %>%      
  dplyr::slice_head(n = 8) %>%      
  dplyr::mutate(GeneRatio = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))) %>%  # 转为数值
  dplyr::arrange(desc(GeneRatio))  
go_for_plot <- top10_go

p2_ggplot <- ggplot(go_for_plot, aes(x = reorder(Description, GeneRatio), y = GeneRatio, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.7) +  
  coord_flip() +  
  scale_size_continuous(range = c(3, 10)) +  
  scale_color_gradient(low = "blue", high = "red") +  
  theme_minimal() +
  labs(title = "Both Upregulated in Urine and Tumor", x = "Gene Ontology enrichment analysis", y = "Gene Ratio", size = "Count", color = "Adjusted p-value") +
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
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    #legend.background = element_rect(fill = "white", color = "black", size = 0.5), 
    #legend.key = element_rect(color = "black", size = 0.5)  
  ) +
  theme(legend.position = "right")
print(p2_ggplot)


ggsave("X_Tissue_DEG_overlapping_up_KEGG.png", width = 7, height = 5,dpi=600,units = "in")
ggsave("X_Tissue_DEG_overlapping_up_KEGG.pdf", width = 7, height = 5,dpi=600,units = "in", device = cairo_pdf)


go_for_plot<-goAll_up
#png("N_Tissue_DEG_overlapping_KEGG.png", width = 300, height = 300, res = 150)  
p2<-clusterProfiler::dotplot(go_for_plot, font.size = 10,
                             title= "GO Dot plot",
                             label_format= 50,
                             showCategory=10)
p2

#ggsave("X_Tissue_DEG_overlapping_up_KEGG.png", width = 8, height = 8,dpi=600)


######X_up_tissue_down
library(ggplot2)
select_X_CH_gene<-X_C_H_DEG[which(rownames(X_C_H_DEG)%in%intersection_gene_X_down),]
select_T_P_gene<-T_P_unpaired_DEG[which(T_P_unpaired_DEG$hgnc_symbol%in%intersection_gene_X_down),]

select_X_CH_gene$ID<-rownames(select_X_CH_gene)
select_T_P_gene$ID<-select_T_P_gene$hgnc_symbol

dataplot<-merge(select_X_CH_gene,select_T_P_gene,by='ID')
library(ggplot2)

result <- cor.test(dataplot$log2FoldChange.x,dataplot$log2FoldChange.y, method = "pearson")
result$estimate  
result$p.value   
p_value <- result$p.value  

ymax <- max(dataplot$log2FoldChange.y, na.rm = TRUE)
ymin <- min(dataplot$log2FoldChange.x, na.rm = TRUE)

xmax <- max(dataplot$log2FoldChange.x, na.rm = TRUE)
xmin <- min(dataplot$log2FoldChange.x, na.rm = TRUE)

ggplot(dataplot, aes(x = log2FoldChange.x, y = log2FoldChange.y)) +
  geom_point(color = "black", alpha = 0.6, size = 2) +
  geom_text(aes(x = (xmax-xmin)*0.8+xmin, y =(ymax-ymin)*0.9+xmin),
            label = paste("P =", format(result$p.value, digits = 3)),
            color = "black", size = 3.5,family = "sans") +
  geom_text(aes(x = (xmax-xmin)*0.8+xmin,y =(ymax-ymin)*0.8+xmin, label = paste("correlation =", format(result$estimate, digits = 3))),
            color = "black", size = 3.5,family = "sans") +  # 显示P值
  labs(title = "", x = "Log2FC_in_Plasma_CfRNA", y = "Log2FC_in_tissue") +
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


ggsave("X_Tissue_DEG_down_cor.png",height=4,width=4,dpi=1200,units = "in")
ggsave("X_Tissue_DEG_down_cor.pdf",height=4,width=4,dpi=1200,units = "in", device = cairo_pdf)
library(clusterProfiler)
goAll_down <- enrichGO(intersection_gene_X_down,
                       keyType = "SYMBOL",
                       OrgDb = "org.Hs.eg.db",
                       ont = "ALL",
                       pAdjustMethod= 'BH',
                       pvalueCutoff= 0.05,
                       qvalueCutoff= 0.2)  

top10_go <- goAll_down@result %>%
  dplyr::arrange(p.adjust) %>%      
  dplyr::slice_head(n = 8) %>%      
  dplyr::mutate(GeneRatio = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))) %>%  # 转为数值
  dplyr::arrange(desc(GeneRatio))  
go_for_plot <- top10_go

p2_ggplot <- ggplot(go_for_plot, aes(x = reorder(Description, GeneRatio), y = GeneRatio, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.7) +  
  coord_flip() +  
  scale_size_continuous(range = c(3, 10)) + 
  scale_color_gradient(low = "blue", high = "red") +  
  theme_minimal() +
  labs(title = "Both Upregulated in Urine and Tumor", x = "Gene Ontology enrichment analysis", y = "Gene Ratio", size = "Count", color = "Adjusted p-value") +
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
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    #legend.background = element_rect(fill = "white", color = "black", size = 0.5), 
    #legend.key = element_rect(color = "black", size = 0.5)  
  ) +
  theme(legend.position = "right")
print(p2_ggplot)
ggsave("X_Tissue_DEG_overlapping_down_KEGG.png", width = 7, height = 5,dpi=600,units = "in")
ggsave("X_Tissue_DEG_overlapping_down_KEGG.pdf", width = 7, height = 5,dpi=600,units = "in", device = cairo_pdf)




go_for_plot<-goAll_down
go_for_plot_graph<-as.data.frame(go_for_plot)
go_for_plot_graph_selected <- go_for_plot_graph[order(go_for_plot_graph$p.adjust, decreasing = F), ]
go_for_plot_graph_selected <-go_for_plot_graph_selected[1:5,]
colnames(go_for_plot_graph_selected)[2]<-'pathway'

gene_list_for_graph<-str_split(go_for_plot_graph_selected$geneID,"\\/")
pathway_list_for_graph<-go_for_plot_graph_selected$pathway

edge_list <- tibble(
  pathway = rep(pathway_list_for_graph, sapply(gene_list_for_graph, length)),
  gene = unlist(gene_list_for_graph)
)

logFC_vector<-X_C_H_DEG[which(rownames(X_C_H_DEG) %in% unique(unlist(gene_list_for_graph))),]
logFC_vector$gene <- rownames(logFC_vector)

edge_list <- left_join(edge_list, go_for_plot_graph_selected, by = "pathway")
edge_list <- left_join(edge_list, logFC_vector, by = "gene")

logFC_vector_2<-T_P_unpaired_DEG[which(T_P_unpaired_DEG$hgnc_symbol %in% unique(unlist(gene_list_for_graph))),-c(1:2,4:6)]
logFC_vector_2$gene <- logFC_vector_2$hgnc_symbol

edge_list <- left_join(edge_list, logFC_vector_2, by = "gene")



sector_colors <- c("myeloid leukocyte activation" = "#F5F5F5", 
                   "regulation of cardiocyte differentiation" = "#FFB6C1", 
                   "positive regulation of T cell migration" = "#B0E0E6", 
                   "T cell migration" = "#FFFACD", 
                   "cellular extravasation" = "#E6E6FA")
circos.clear()

edge_list_plot<-edge_list[,c(1,2)]

logfc_col_fun <- colorRamp2(
  breaks = c(-8,0, 8),                      
  colors = c("blue","white", "red")        
)

sector_colors_grid <- setNames(
  logfc_col_fun(edge_list$log2FoldChange.x),
  edge_list$gene
)

sector_colors_combined <- c(sector_colors_grid, sector_colors)


#CairoPNG("X_Tissue_DEG_down_chord.png",width = 7, height = 5,res=600,units = "in")  
library(Cairo)

CairoPDF("X_Tissue_DEG_down_chord.pdf",width = 7, height = 5,family = "Arial")  
link_colors <- sector_colors[edge_list_plot$pathway]
chordDiagram(edge_list_plot, 
             col = link_colors,  
             annotationTrack = c("grid"),  
             transparency = 0.1,  
             grid.col = sector_colors_combined
)

circos.trackPlotRegion(
  track.index = 1,
  bg.border = NA,
  panel.fun = function(x, y) {
    sector_name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    
    label_lines <- strwrap(sector_name, width = 40)
    n_lines <- length(label_lines)
    
    line_spacing <- 0.4
    
    total_height <- (n_lines - 1) * line_spacing
    start_y <- ylim[2] + 0.5 + total_height / 2  
    
    for (i in seq_along(label_lines)) {
      circos.text(
        x = mean(xlim),
        y = start_y - (i - 1) * line_spacing,
        labels = label_lines[i],
        facing = "bending.inside",
        niceFacing = FALSE,
        col = "black",
        cex = 0.3,
        family = "Arial"
      )
    }
  }
)

col_fun <- colorRamp2(c(-8,0,8), c("blue","white", "red"))
sector_colors <- setNames(col_fun(edge_list$log2FoldChange.y), edge_list_plot$gene)

circos.track(
  track.index = 2,  
  ylim = c(0, 1),
  bg.border = NA,
  track.height = 0.08,
  panel.fun = function(x, y) {
    sector_name <- get.cell.meta.data("sector.index")
    if (!(sector_name %in% names(sector_colors))) {
      return(NULL)
    }
    color <- sector_colors[sector_name]
    circos.rect(xleft = CELL_META$xlim[1],
                xright = CELL_META$xlim[2],
                ybottom = 0,
                ytop = 1,
                col = color,
                border = NA)
  }
)
lgd_logfc <- Legend(col_fun = logfc_col_fun,
                    title = "logFC",
                    at = c(-8, 0, 8),
                    title_gp = gpar(fontsize = 10, fontface = "plain", col = "black"), 
                    labels_gp = gpar(fontsize = 10, fontface = "plain", col = "black"),
                    labels = c("-8", "0", "8"))

draw(lgd_logfc, 
     x = unit(1, "npc") - unit(10, "mm"), 
     y = unit(0.5, "npc")) 
dev.off()