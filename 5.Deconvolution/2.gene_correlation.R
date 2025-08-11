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
all.equal(rownames(paired_T_TPM_tissue),rownames(paired_P_TPM_tissue),
          rownames(paired_X_TPM_cfRNA),rownames(paired_N_TPM_cfRNA))

all.equal(colnames(paired_T_TPM_tissue),colnames(paired_P_TPM_tissue),
          colnames(paired_X_TPM_cfRNA),colnames(paired_N_TPM_cfRNA))

log2_paired_T_TPM_tissue <- log2(paired_T_TPM_tissue + 0.1)
log2_paired_P_TPM_tissue <- log2(paired_P_TPM_tissue + 0.1)
log2_paired_X_TPM_cfRNA <- log2(paired_X_TPM_cfRNA + 0.1)
log2_paired_N_TPM_cfRNA <- log2(paired_N_TPM_cfRNA + 0.1)


cor_list<-matrix(ncol = 6,nrow = 61228)
rownames(cor_list)<-rownames(paired_T_TPM_tissue)
colnames(cor_list)<-c("Tumor_NAT","Tumor_Plasma","Tumor_Urine","NAT_Plasma","NAT_Urine","Plasma_Urine")
p_list<-matrix(ncol = 6,nrow = 61228)
rownames(p_list)<-rownames(paired_T_TPM_tissue)
colnames(p_list)<-c("Tumor_NAT","Tumor_Plasma","Tumor_Urine","NAT_Plasma","NAT_Urine","Plasma_Urine")
for (k in 1:61228) {
  
  # T vs P
  if (!(all(paired_T_TPM_tissue[k,] == 0) | all(paired_P_TPM_tissue[k,] == 0))) {
    result_TP<-cor.test(t(log2_paired_T_TPM_tissue[k,]), t(log2_paired_P_TPM_tissue[k,]), method = 'pearson')
    cor_list[k,1] <- result_TP$estimate
    p_list[k,1] <- result_TP$p.value
    
  }
  
  # T vs X
  if (!(all(paired_T_TPM_tissue[k,] == 0) | all(paired_X_TPM_cfRNA[k,] == 0))) {
    result_TX<-cor.test(t(log2_paired_T_TPM_tissue[k,]), t(log2_paired_X_TPM_cfRNA[k,]), method = 'pearson')
    cor_list[k,2] <- result_TX$estimate
    p_list[k,2] <- result_TX$p.value
  }
  
  # T vs N
  if (!(all(paired_T_TPM_tissue[k,] == 0) | all(paired_N_TPM_cfRNA[k,] == 0))) {
    result_TN<-cor.test(t(log2_paired_T_TPM_tissue[k,]), t(log2_paired_N_TPM_cfRNA[k,]), method = 'pearson')
    cor_list[k,3] <- result_TN$estimate
    p_list[k,3] <- result_TN$p.value
  }
  
  # P vs X
  if (!(all(paired_P_TPM_tissue[k,] == 0) | all(paired_X_TPM_cfRNA[k,] == 0))) {
    result_PX<-cor.test(t(log2_paired_P_TPM_tissue[k,]), t(log2_paired_X_TPM_cfRNA[k,]), method = 'pearson')
    cor_list[k,4] <- result_PX$estimate
    p_list[k,4] <- result_PX$p.value
  }
  
  # P vs N
  if (!(all(paired_P_TPM_tissue[k,] == 0) | all(paired_N_TPM_cfRNA[k,] == 0))) {
    result_PN<-cor.test(t(log2_paired_P_TPM_tissue[k,]), t(log2_paired_N_TPM_cfRNA[k,]), method = 'pearson')
    cor_list[k,5] <- result_PN$estimate
    p_list[k,5] <- result_PN$p.value
  }
  
  # X vs N
  if (!(all(paired_X_TPM_cfRNA[k,] == 0) | all(paired_N_TPM_cfRNA[k,] == 0))) {
    result_XN<-cor.test(t(log2_paired_X_TPM_cfRNA[k,]), t(log2_paired_N_TPM_cfRNA[k,]), method = 'pearson')
    cor_list[k,6] <- result_XN$estimate
    p_list[k,6] <- result_XN$p.value
  }
  
}

write.csv(cor_list,"gene_cor_list.csv")
write.csv(p_list,"gene_p_list.csv")

colSums(!is.na(cor_list))


cor_list<-read.csv("gene_cor_list.csv",header=T,row.names=1)
p_list<-read.csv("gene_p_list.csv",header=T,row.names=1)


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
basemeanthreshold<-1
tissuelog2FoldChangethreshold<-log2(2)
N_C_H_DEG_rownames<-rownames(N_C_H_DEG[which(N_C_H_DEG$padj<padjthreshold&abs(N_C_H_DEG$log2FoldChange)>log2FoldChangethreshold&N_C_H_DEG$baseMean>basemeanthreshold),])
N_C_S_DEG_rownames<-rownames(N_C_S_DEG[which(N_C_S_DEG$padj<padjthreshold&abs(N_C_S_DEG$log2FoldChange)>log2FoldChangethreshold&N_C_S_DEG$baseMean>basemeanthreshold),])
X_C_H_DEG_rownames<-rownames(X_C_H_DEG[which(X_C_H_DEG$padj<padjthreshold&abs(X_C_H_DEG$log2FoldChange)>log2FoldChangethreshold&X_C_H_DEG$baseMean>basemeanthreshold),])
X_C_S_DEG_rownames<-rownames(X_C_S_DEG[which(X_C_S_DEG$padj<padjthreshold&abs(X_C_S_DEG$log2FoldChange)>log2FoldChangethreshold&X_C_S_DEG$baseMean>basemeanthreshold),])
T_P_unpaired_DEG<-T_P_unpaired_DEG[which(T_P_unpaired_DEG$padj<padjthreshold&T_P_unpaired_DEG$log2FoldChange>tissuelog2FoldChangethreshold&T_P_unpaired_DEG$baseMean>0&T_P_unpaired_DEG$hgnc_symbol!= ''),7]

genelist<-rownames(N_C_H_DEG[which(N_C_H_DEG$baseMean>basemeanthreshold),])


data_for_plot<-cbind(cor_list[,c(5)],p_list[,c(5)])
rownames(data_for_plot)<-rownames(cor_list)
data_for_plot<-data_for_plot[which(rownames(data_for_plot) %in% genelist),]
data_for_plot<-as.data.frame(data_for_plot)
colSums(!is.na(data_for_plot))

colnames(data_for_plot)<-c("correlation","p")
data_for_plot$ID<-rownames(data_for_plot)

data_for_plot$Significance <- 'NS'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_H_DEG_rownames &!(ID %in% N_C_S_DEG_rownames))] <- 'N_C_H_only_DEG'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_S_DEG_rownames &!(ID %in% N_C_H_DEG_rownames))] <- 'N_C_S_only_DEG'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_H_DEG_rownames & ID %in% N_C_S_DEG_rownames)] <- 'both_DEG'

dim(data_for_plot[which(data_for_plot$p<0.01&data_for_plot$correlation>0),])
dim(data_for_plot[which(data_for_plot$p<0.01&data_for_plot$correlation<0),])
dim(data_for_plot[which(data_for_plot$p<0.01&data_for_plot$correlation>0&data_for_plot$Significance=='N_C_H_only_DEG'),])
dim(data_for_plot[which(data_for_plot$p<0.01&data_for_plot$correlation>0&data_for_plot$Significance=='N_C_S_only_DEG'),])
dim(data_for_plot[which(data_for_plot$p<0.01&data_for_plot$correlation>0&data_for_plot$Significance=='both_DEG'),])

data_for_plot[which(data_for_plot$p<0.01&data_for_plot$correlation>0&data_for_plot$Significance=='N_C_H_only_DEG'),]
data_for_plot[which(data_for_plot$p<0.01&data_for_plot$correlation>0&data_for_plot$Significance=='N_C_S_only_DEG'),]
data_for_plot[which(data_for_plot$p<0.01&data_for_plot$correlation>0&data_for_plot$Significance=='both_DEG'),]



cols <- c('N_C_H_only_DEG'="darkgreen", 
          'N_C_S_only_DEG'='orange',
          'both_DEG'="darkblue",
          'NS'='gray')
p1 <- ggplot(data_for_plot, aes(x = correlation, y = -log10(p), color = Significance)) +
  labs(x = expression('correlation'), 
       y = expression('-Log'['10']*'(pvalue)'), 
       title = NULL) +
  geom_point(size = 2, shape = 16, alpha = 0.8) +   # <<<<<< 实心点，shape=16
  scale_color_manual(values = cols) +
  #geom_vline(xintercept = c(-1, 1), color = 'darkgreen', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), color = 'darkgreen', linetype = 'dashed') +
  geom_text_repel(data = subset(data_for_plot, (-log10(p) > 2.8 & correlation > 0.5 & Significance != "NS")),
                  aes(label = ID), nudge_x = 0.1,
                  segment.alpha = 0.5, min.segment.length = 0.3,max.overlaps = 50,
                  size = 3, color = 'black', segment.color = 'black') +
  geom_text_repel(data = subset(data_for_plot, correlation < (-0.5) & -log10(p) > 2.8 & Significance != "NS"),
                  aes(label = ID),
                  segment.alpha = 0.5, min.segment.length = 0.3, max.overlaps = 50,
                  size = 3, color = 'black', segment.color = 'black') +
 theme_minimal()+theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10,family='Arial'),
    legend.text = element_text(size = 10,family='Arial'),
    axis.title = element_text(size = 10, margin = margin(t = 5),family='Arial'),
    axis.title.y = element_text(margin = margin(r = 0),family='Arial'),
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
  )
p1

ggsave("NAT_N.png",height=5,width=6,units='in',dpi=1200)
ggsave("NAT_N.pdf",height=5,width=6,units='in',device = cairo_pdf,family='Arial')


data_for_plot<-cbind(cor_list[,c(3)],p_list[,c(3)])
rownames(data_for_plot)<-rownames(cor_list)
data_for_plot<-data_for_plot[which(rownames(data_for_plot) %in% genelist),]
colSums(!is.na(data_for_plot))

data_for_plot<-as.data.frame(data_for_plot)
colnames(data_for_plot)<-c("correlation","p")
data_for_plot$ID<-rownames(data_for_plot)
data_for_plot$Significance <- 'NS'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_H_DEG_rownames &!(ID %in% N_C_S_DEG_rownames))] <- 'N_C_H_only_DEG'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_S_DEG_rownames &!(ID %in% N_C_H_DEG_rownames))] <- 'N_C_S_only_DEG'
data_for_plot$Significance[with(data_for_plot,ID %in% N_C_H_DEG_rownames & ID %in% N_C_S_DEG_rownames)] <- 'both_DEG'

dim(data_for_plot[which(data_for_plot$p<0.01&data_for_plot$correlation>0),])
dim(data_for_plot[which(data_for_plot$p<0.01&data_for_plot$correlation<0),])
dim(data_for_plot[which(data_for_plot$p<0.01&data_for_plot$correlation>0&data_for_plot$Significance=='N_C_H_only_DEG'),])
dim(data_for_plot[which(data_for_plot$p<0.01&data_for_plot$correlation>0&data_for_plot$Significance=='N_C_S_only_DEG'),])
dim(data_for_plot[which(data_for_plot$p<0.01&data_for_plot$correlation>0&data_for_plot$Significance=='both_DEG'),])

data_for_plot[which(data_for_plot$p<0.01&data_for_plot$correlation>0&data_for_plot$Significance=='N_C_H_only_DEG'),]
data_for_plot[which(data_for_plot$p<0.01&data_for_plot$correlation>0&data_for_plot$Significance=='N_C_S_only_DEG'),]
data_for_plot[which(data_for_plot$p<0.01&data_for_plot$correlation>0&data_for_plot$Significance=='both_DEG'),]


cols <- c('N_C_H_only_DEG'="darkgreen", 
          'N_C_S_only_DEG'='orange',
          'both_DEG'="darkblue",
          'NS'='gray')
p2 <- ggplot(data_for_plot, aes(x = correlation, y = -log10(p), color = Significance)) +
  labs(x = expression('correlation'), 
       y = expression('-Log'['10']*'(pvalue)'), 
       title = NULL) +
  geom_point(size = 2, shape = 16, alpha = 0.8) +   # <<<<<< 实心点，shape=16
  scale_color_manual(values = cols) +
  #geom_vline(xintercept = c(-1, 1), color = 'darkgreen', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), color = 'darkgreen', linetype = 'dashed') +
  geom_text_repel(data = subset(data_for_plot, (-log10(p) > 2.8 & correlation > 0.5 & Significance != "NS")),
                  aes(label = ID), nudge_x = 0.1,
                  segment.alpha = 0.5, min.segment.length = 0.3,max.overlaps = 50,
                  size = 3, color = 'black', segment.color = 'black') +
  geom_text_repel(data = subset(data_for_plot, correlation < (-0.5) & -log10(p) > 2.8 & Significance != "NS"),
                  aes(label = ID),
                  segment.alpha = 0.5, min.segment.length = 0.3, max.overlaps = 50,
                  size = 3, color = 'black', segment.color = 'black') +
  theme_minimal()+theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10,family='Arial'),
    legend.text = element_text(size = 10,family='Arial'),
    axis.title = element_text(size = 10, margin = margin(t = 5),family='Arial'),
    axis.title.y = element_text(margin = margin(r = 0),family='Arial'),
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
  )
p2
ggsave("Tumor_N.png",height=5,width=6,units='in',dpi=1200)
ggsave("Tumor_N.pdf",height=5,width=6,units='in',device = cairo_pdf,family='Arial')
