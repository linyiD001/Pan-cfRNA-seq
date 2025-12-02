genelist<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome/all_genelist_annotation.csv",header=T,row.names = 1)
genelist<-genelist[,4:5]

genelist[grep("DNA",genelist$genetype),2]<-"Other"
genelist[grep("IG",genelist$genetype),2]<-"Other"
genelist[grep("TR",genelist$genetype),2]<-"Other"
genelist[grep("LTR",genelist$genetype),2]<-"LTR"
genelist[grep("rRNA",genelist$genetype),2]<-"rRNA"
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
  arrange(genetype == "Other") %>%  # 把 "other" 排在后面
  ungroup()

genelist_clean<-genelist_clean[!duplicated(genelist_clean$gene_symbol),]


DEG<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/DEG_repeataware/salmon_repeataware_stranded_allgenome/lfc_X_type_C_vs_H_result.csv",header=T)

colnames(DEG)[1]<-c("gene_symbol")
data_merge<-merge(DEG,genelist_clean,by='gene_symbol',all.x=T)

DEG<-data_merge[which(data_merge$pvalue<0.05&(abs(data_merge$log2FoldChange)>2)),]

length(unique(DEG[which(DEG$genetype=='Simple_repeat'),1]))
unique(DEG[which(DEG$genetype=='Satellite'),1])
unique(DEG[which(DEG$genetype=='Low_complexity'),1])
unique(DEG[which(DEG$genetype=='SINE'),1])
unique(DEG[which(DEG$genetype=='LINE'),1])
unique(DEG[which(DEG$genetype=='RC'),1])
unique(DEG[which(DEG$genetype=='LTR'),1])

RE_list<-c("Simple_repeat","Satellite","Low_complexity","SINE","LINE","RC","LTR")

length(unique(DEG[which(DEG$genetype%in%RE_list),1]))



RE_list<-c("Simple_repeat","Satellite","Low_complexity","SINE","LINE","RC","LTR")

library(ggrepel)
library(ggplot2)

p <- ggplot(data_merge, aes(x = log2FoldChange, y = -log10(padj))) +
  labs(x = expression('Log'['2']*'(Fold Change)'), 
       y = expression('-Log'['10']*'(FDR)'), 
       color = "Biotype" ) +
  # Color all points gray by default
  geom_point(alpha = 0.6, size = 2, shape = 16, color = "gray") +
  # Color only the points where log2FoldChange > 2 based on their genetype
  geom_point(data = subset(data_merge,abs(log2FoldChange)>2&(-log10(padj)>2)), 
             aes(color = genetype), 
             alpha = 0.6, size = 2, shape = 16) +
  geom_vline(xintercept = c(-2, 2), color = 'darkgreen', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), color = 'darkgreen', linetype = 'dashed') +
  geom_text_repel(data = subset(data_merge, (-log10(padj) > 2 & log2FoldChange > 6)),
                  aes(label = gene_symbol), nudge_x = 0.1,
                  segment.alpha = 0.5, min.segment.length = 0.2,
                  size = 3, color = 'black', segment.color = 'black') +
  geom_text_repel(data = subset(data_merge, log2FoldChange < -6 & -log10(padj) > 2),
                  aes(label = gene_symbol),
                  segment.alpha = 0.5, min.segment.length = 0.2, max.overlaps = 100,
                  size = 3, color = 'black', segment.color = 'black') +
  scale_color_manual(values = c('LINE' = 'darkred', 
                                'SINE' = 'blue', 
                                'protein_coding' = 'darkgreen', 
                                'Simple_repeat' = 'pink', 
                                'LTR' = 'orange', 
                                'Low_complexity' = 'purple',
                                'lncRNA'='lightblue')) +
  theme_bw() + 
  theme(legend.position = "right") +
  theme(axis.text = element_text(size = 10),
        strip.text = element_text(size = 10, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 10, margin = margin(t = 5)),
        axis.title.y = element_text(margin = margin(r = 5)),
        axis.text.x = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.x = unit(0.5, 'cm'),   
        legend.spacing.y = unit(0.5, 'cm'),   
        legend.key.height = unit(0.5, "cm"),    
        legend.key.spacing.y = unit(0.1, "cm"),    
            legend.key.width = unit(0.5, "cm"),     
            legend.box.spacing = unit(0.5, 'cm'))
p


ggsave("X_type_C_vs_H_result_volcano_genetype.png",height=5,width=7, dpi=1200)
ggsave("X_type_C_vs_H_result_volcano_genetype.pdf",height=5,width=7,dpi=1200, units = "in", device = cairo_pdf)

