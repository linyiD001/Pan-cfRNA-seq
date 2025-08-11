###for repeataware
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
library(stringr)
library(viridis)
library(ggplot2)
library(tidyverse)
library(colorspace)    
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(patchwork)

setwd("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome")
tximeta::loadLinkedTxome("/dssg/home/acct-dahan/share/references/gencode/RE_allgenome/gencode.v44.index.repeat.aware.json")
name<-read.table("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/cfRNA.txt",header = F,sep = "\t",fill=T)
name<-as.vector(name$V1)

files <- file.path("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome", name, "quant.sf") 
#file.exists(files)
additional_list<-name[which(file.exists(files)==TRUE)]

#repeataware_all_count <- read.csv("salmon.repeataware/repeataware_all_count.csv",row.names = 1)
files <- file.path("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome", additional_list, "quant.sf") 
coldata <- data.frame(files, names=name, condition="ISR", stringsAsFactors=FALSE)

se <- tximeta(coldata)
#se <- tximeta(coldata,skipMeta=TRUE)
data<-as.data.frame(assay(se))

res3<-colnames(data)
data_1<-data
colnames(data_1)<-res3
write.csv(data_1,'cfRNA_repeataware_count.csv')

countdata<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome/cfRNA_repeataware_count.csv")
rmsk_genelist<-read.csv("/dssg/home/acct-dahan/share/references/gencode/RE_allgenome/gencode.v44.ucsc.rmsk.tx.to.gene.csv",header=F)

res1 <- str_split(rmsk_genelist[,1],"\\|")
res2 <- as.data.frame(res1)
head(res2)
res3<-as.data.frame(t(res2)[,c(1,2,6,8)])
res4<-cbind(rmsk_genelist[,1],res3)
colnames(res4)<-c("rowname","enst","ensg","gene_symbol","genetype")

all_genelist_test<-read.csv("/dssg/home/acct-dahan/share/references/gencode/RE_allgenome/ucsc.rm/rmsk.csv",row.names = 1,header=F)
colnames(all_genelist_test)<-c('rowname','gene_symbol','V3','genetype')

ensg<-res4[grep('ENSG',res4$ensg),]
RE<-res4[-grep('ENSG',res4$ensg),1:3]

RE<-merge(RE,all_genelist_test[,c(1,2,4)],all.x=T,by='rowname')
all_genelist_combine<-rbind(ensg,RE)

write.csv(all_genelist_combine,"all_genelist_annotation.csv")

#all_genelist_combine<-read.csv("all_genelist_annotation.csv",header=T,row.names=1)
countdata$enst<-countdata$X
countdata<-countdata[,-1]
countdata_merge<-merge(countdata,all_genelist_combine,by='enst',all.x=T)

genename_count<- aggregate(countdata_merge[,2:277],by=list(countdata_merge$gene_symbol), FUN=sum)
write.csv(genename_count,'salmon.repeataware_simple_genename.csv')
genetype_count<- aggregate(countdata_merge[,2:277],by=list(countdata_merge$genetype), FUN=sum)
write.csv(genetype_count,"genetype_cfRNA.csv")

#genetype_count<-read.csv("genetype_cfRNA.csv",header = T,row.names=1)

percent_count<-as.data.frame(apply(genetype_count[,2:277], 2, function (x) x/sum(x)*100 ))
percent_count<-cbind(genetype_count$Group.1,percent_count)

rearranged_genetype<-as.data.frame(percent_count$`genetype_count$Group.1`)
colnames(rearranged_genetype)[1]<-'genetype'
rearranged_genetype$new_type<-rearranged_genetype$genetype

rearranged_genetype[grep("DNA",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("IG",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("TR",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("LTR",rearranged_genetype$genetype),2]<-"LTR"
rearranged_genetype[grep("rRNA",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("tRNA",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("decay",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("lncRNA",rearranged_genetype$genetype),2]<-"lncRNA"
rearranged_genetype[grep("ribozyme",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[which(rearranged_genetype$genetype=='RNA'),2]<-"Other"
rearranged_genetype[grep("protein_coding",rearranged_genetype$genetype),2]<-"protein_coding"
rearranged_genetype[grep("RC",rearranged_genetype$genetype),2]<-"RC"
rearranged_genetype[grep("pseudogene",rearranged_genetype$genetype),2]<-"pseudogene"
rearranged_genetype[grep("SINE",rearranged_genetype$genetype),2]<-"SINE"
rearranged_genetype[grep("artifact",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("scaRNA",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("scRNA",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("snoRNA",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("snRNA",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("sRNA",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("srpRNA",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("Unknown",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("vault_RNA",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("TEC",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("retained_intron",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("processed_transcript",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("misc_RNA",rearranged_genetype$genetype),2]<-"Other"
rearranged_genetype[grep("miRNA",rearranged_genetype$genetype),2]<-"Other"

colnames(percent_count)[1]<-'genetype'
percent_count_rearranged<-merge(percent_count,rearranged_genetype,by='genetype',all=T)
percent_count_rearranged<- aggregate(percent_count_rearranged[,2:277],by=list(percent_count_rearranged$new_type), FUN=sum)

rownames(percent_count_rearranged)<-percent_count_rearranged$Group.1
percent_count_rearranged_1<-percent_count_rearranged[,-1]

write.csv(percent_count,"genetype_percent_cfRNA.csv")
write.csv(percent_count_rearranged_1,"genetype_percent_cfRNA_rearranged.csv")

#percent_count_rearranged_1<-read.csv("genetype_percent_cfRNA_rearranged.csv",header=T,row.names = 1)
QC_list<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/QC_repeataware_list_withoutR.csv",row.names = 1)
percent_count_rearranged_1<-percent_count_rearranged_1[,which(colnames(percent_count_rearranged_1) %in% QC_list[,1])]
#write.csv(percent_count_rearranged_1,"genetype_percent_cfRNA_rearranged_filtered.csv")

###Urine
#############
N_data<-percent_count_rearranged_1[,grep("N",colnames(percent_count_rearranged_1))]
df_long <-N_data %>%
  rownames_to_column(var = "Biotype") %>% 
  pivot_longer(cols = -Biotype,           
               names_to = "sample",
               values_to = "percentage")

df_long$group<-NA
df_long[grep("^C", df_long$sample),4]<-"BC"
df_long[grep("^S", df_long$sample),4]<-"BS"
df_long[grep("^H", df_long$sample),4]<-"HC"

df_line <- df_long %>% 
  filter(Biotype == "LINE")
df_line_summary <- df_line %>%
  group_by(group) %>%
  dplyr::summarise(
    mean_pct = mean(percentage),
    sd_pct = sd(percentage)
  )
custom_colors <- c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4")

y_max <- max(df_line$percentage)

urine_LINE_plot<-ggplot()+geom_bar(
    data = df_line_summary,
    aes(x = group, y = mean_pct, fill = group),
    stat = "identity", width = 0.6, color = "black",
    alpha = 0.6
  ) +
  geom_errorbar(
    data = df_line_summary,
    aes(x = group, ymin = mean_pct,ymax = mean_pct + sd_pct),
    width = 0.2
  ) +）
  geom_jitter(
    data = df_line,
    aes(x = group, y = percentage, color = group),
    width = 0.15, size = 2, alpha = 0.8,
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(
    x = "Group",
    y = "Percentage of LINE in Urine CfRNA",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 8, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.line = element_line(color = "black"),  
    axis.text.x = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(0.5, 'cm'),   
    legend.key.height = unit(0.5, "cm"),    
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.width = unit(0.5, "cm"),     
    legend.box.spacing = unit(0.5, 'cm'),
    legend.key = element_rect(color = "black", size = 0.1),  
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  scale_y_continuous(expand = c(0, 0))+
  annotate("text", x = 2.5, y = y_max*1.1, label = "****",size=5) +
  expand_limits(y = y_max*1.3) + 
  annotate("segment", x = 2, xend = 3, y = y_max*1.08, yend = y_max*1.08,linewidth=0.3)  +
  annotate("text", x = 2, y = y_max*1.2, label = "***",size=5) +
  annotate("segment", x = 1, xend = 3, y = y_max*1.18, yend = y_max*1.18,linewidth=0.3)  

#urine_LINE_plot


df_line <- df_long %>% 
  filter(Biotype == "SINE")
df_line_summary <- df_line %>%
  group_by(group) %>%
  dplyr::summarise(
    mean_pct = mean(percentage),
    sd_pct = sd(percentage)
  )
custom_colors <- c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4")

y_max <- max(df_line$percentage)

urine_SINE_plot<-ggplot()+geom_bar(
  data = df_line_summary,
  aes(x = group, y = mean_pct, fill = group),
  stat = "identity", width = 0.6, color = "black",
  alpha = 0.6
) +
  geom_errorbar(
    data = df_line_summary,
    aes(x = group, ymin = mean_pct,ymax = mean_pct + sd_pct),
    width = 0.2
  ) +
  geom_jitter(
    data = df_line,
    aes(x = group, y = percentage, color = group),
    width = 0.15, size = 2, alpha = 0.8,
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(
    x = "Group",
    y = "Percentage of SINE in Urine CfRNA",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 8, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.line = element_line(color = "black"),  
    axis.text.x = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(0.5, 'cm'),   
    legend.key.height = unit(0.5, "cm"),    
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.width = unit(0.5, "cm"),    
    legend.box.spacing = unit(0.5, 'cm'),
    legend.key = element_rect(color = "black", size = 0.1),  
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  scale_y_continuous(expand = c(0, 0))+
  annotate("text", x = 2.5, y = y_max*1.1, label = "****",size=5) +
  expand_limits(y = y_max*1.3) + 
  annotate("segment", x = 2, xend = 3, y = y_max*1.08, yend = y_max*1.08,linewidth=0.3)  +
  annotate("text", x = 2, y = y_max*1.2, label = "***",size=5) +
  annotate("segment", x = 1, xend = 3, y = y_max*1.18, yend = y_max*1.18,linewidth=0.3)  

urine_SINE_plot


df_line <- df_long %>% 
  filter(Biotype == "Simple_repeat")
df_line_summary <- df_line %>%
  group_by(group) %>%
  dplyr::summarise(
    mean_pct = mean(percentage),
    sd_pct = sd(percentage)
  )
custom_colors <- c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4")

y_max <- max(df_line$percentage)

urine_Simple_repeat_plot<-ggplot()+geom_bar(
  data = df_line_summary,
  aes(x = group, y = mean_pct, fill = group),
  stat = "identity", width = 0.6, color = "black",
  alpha = 0.6
) +
  geom_errorbar(
    data = df_line_summary,
    aes(x = group, ymin = mean_pct,ymax = mean_pct + sd_pct),
    width = 0.2
  ) +
  geom_jitter(
    data = df_line,
    aes(x = group, y = percentage, color = group),
    width = 0.15, size = 2, alpha = 0.8,
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(
    x = "Group",
    y = "Percentage of Simple repeat in Urine CfRNA",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 8, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.line = element_line(color = "black"),  
    axis.text.x = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(0.5, 'cm'),   
    legend.key.height = unit(0.5, "cm"),    
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.width = unit(0.5, "cm"),     #
    legend.box.spacing = unit(0.5, 'cm'),
    legend.key = element_rect(color = "black", size = 0.1),  
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  scale_y_continuous(expand = c(0, 0))+
  expand_limits(y = y_max*1.3) + 
  annotate("text", x = 2, y = y_max*1.2, label = "****",size=5) +
  annotate("segment", x = 1, xend = 3, y = y_max*1.18, yend = y_max*1.18,linewidth=0.3)  

#urine_Simple_repeat_plot


df_line <- df_long %>% 
  filter(Biotype == "LTR")
df_line_summary <- df_line %>%
  group_by(group) %>%
  dplyr::summarise(
    mean_pct = mean(percentage),
    sd_pct = sd(percentage)
  )
custom_colors <- c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4")

y_max <- max(df_line$percentage)

urine_LTR_plot<-ggplot()+geom_bar(
  data = df_line_summary,
  aes(x = group, y = mean_pct, fill = group),
  stat = "identity", width = 0.6, color = "black",
  alpha = 0.6
) +
  geom_errorbar(
    data = df_line_summary,
    aes(x = group, ymin = mean_pct,ymax = mean_pct + sd_pct),
    width = 0.2
  ) +
  geom_jitter(
    data = df_line,
    aes(x = group, y = percentage, color = group),
    width = 0.15, size = 2, alpha = 0.8,
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(
    x = "Group",
    y = "Percentage of LTR in Urine CfRNA",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 8, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.line = element_line(color = "black"),  
    axis.text.x = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(0.5, 'cm'),   
    legend.key.height = unit(0.5, "cm"),    
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.width = unit(0.5, "cm"),    
    legend.box.spacing = unit(0.5, 'cm'),
    legend.key = element_rect(color = "black", size = 0.1),  
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  scale_y_continuous(expand = c(0, 0))+
  annotate("text", x = 2.5, y = y_max*1.1, label = "****",size=5) +
  expand_limits(y = y_max*1.3) + 
  annotate("segment", x = 2, xend = 3, y = y_max*1.08, yend = y_max*1.08,linewidth=0.3)  +
  annotate("text", x = 2, y = y_max*1.2, label = "***",size=5) +
  annotate("segment", x = 1, xend = 3, y = y_max*1.18, yend = y_max*1.18,linewidth=0.3)  

#urine_LTR_plot


###Plasma
#############
X_data<-percent_count_rearranged_1[,grep("X",colnames(percent_count_rearranged_1))]
df_long <-X_data %>%
  rownames_to_column(var = "Biotype") %>% 
  pivot_longer(cols = -Biotype,           
               names_to = "sample",
               values_to = "percentage")

df_long$group<-NA
df_long[grep("^C", df_long$sample),4]<-"BC"
df_long[grep("^S", df_long$sample),4]<-"BS"
df_long[grep("^H", df_long$sample),4]<-"HC"

df_line <- df_long %>% 
  filter(Biotype == "LINE")
df_line_summary <- df_line %>%
  group_by(group) %>%
  dplyr::summarise(
    mean_pct = mean(percentage),
    sd_pct = sd(percentage)
  )
custom_colors <- c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4")
y_max <- max(df_line$percentage)


Plasma_LINE_plot<-ggplot()+geom_bar(
  data = df_line_summary,
  aes(x = group, y = mean_pct, fill = group),
  stat = "identity", width = 0.6, color = "black",
  alpha = 0.6
) +
  geom_errorbar(
    data = df_line_summary,
    aes(x = group, ymin = mean_pct,ymax = mean_pct + sd_pct),
    width = 0.2
  ) +
  geom_jitter(
    data = df_line,
    aes(x = group, y = percentage, color = group),
    width = 0.15, size = 2, alpha = 0.8,
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(
    x = "Group",
    y = "Percentage of LINE in Plasma CfRNA",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 8, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.line = element_line(color = "black"), 
    axis.text.x = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(0.5, 'cm'),   
    legend.key.height = unit(0.5, "cm"),    
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.width = unit(0.5, "cm"),     
    legend.box.spacing = unit(0.5, 'cm'),
    legend.key = element_rect(color = "black", size = 0.1),  
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  scale_y_continuous(expand = c(0, 0))+
  annotate("text", x = 2.5, y = y_max*1.1, label = "**",size=5) +
  expand_limits(y = y_max*1.3) + 
  annotate("segment", x = 2, xend = 3, y = y_max*1.08, yend = y_max*1.08,linewidth=0.3)  +
  annotate("text", x = 2, y = y_max*1.2, label = "**",size=5) +
  annotate("segment", x = 1, xend = 3, y = y_max*1.18, yend = y_max*1.18,linewidth=0.3)  
#Plasma_LINE_plot



df_line <- df_long %>% 
  filter(Biotype == "SINE")
df_line_summary <- df_line %>%
  group_by(group) %>%
  dplyr::summarise(
    mean_pct = mean(percentage),
    sd_pct = sd(percentage)
  )
custom_colors <- c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4")

y_max <- max(df_line$percentage)
Plasma_SINE_plot<-ggplot()+geom_bar(
  data = df_line_summary,
  aes(x = group, y = mean_pct, fill = group),
  stat = "identity", width = 0.6, color = "black",
  alpha = 0.6
) +
  geom_errorbar(
    data = df_line_summary,
    aes(x = group, ymin = mean_pct,ymax = mean_pct + sd_pct),
    width = 0.2
  ) +
  geom_jitter(
    data = df_line,
    aes(x = group, y = percentage, color = group),
    width = 0.15, size = 2, alpha = 0.8,
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(
    x = "Group",
    y = "Percentage of SINE in Plasma CfRNA",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 8, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.line = element_line(color = "black"),  
    axis.text.x = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(0.5, 'cm'),   
    legend.key.height = unit(0.5, "cm"),    
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.width = unit(0.5, "cm"),     
    legend.box.spacing = unit(0.5, 'cm'),
    legend.key = element_rect(color = "black", size = 0.1),  
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  scale_y_continuous(expand = c(0, 0))+
  annotate("text", x = 2.5, y = y_max*1.1, label = "***",size=5) +
  expand_limits(y = y_max*1.3) +
  annotate("segment", x = 2, xend = 3, y = y_max*1.08, yend = y_max*1.08,linewidth=0.3)  +
  annotate("text", x = 2, y = y_max*1.2, label = "****",size=5) +
  annotate("segment", x = 1, xend = 3, y = y_max*1.18, yend = y_max*1.18,linewidth=0.3)  
#Plasma_SINE_plot


df_line <- df_long %>% 
  filter(Biotype == "Simple_repeat")
df_line_summary <- df_line %>%
  group_by(group) %>%
  dplyr::summarise(
    mean_pct = mean(percentage),
    sd_pct = sd(percentage)
  )
custom_colors <- c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4")

y_max <- max(df_line$percentage)

Plasma_Simple_repeat_plot<-ggplot()+geom_bar(
  data = df_line_summary,
  aes(x = group, y = mean_pct, fill = group),
  stat = "identity", width = 0.6, color = "black",
  alpha = 0.6
) +
  geom_errorbar(
    data = df_line_summary,
    aes(x = group, ymin = mean_pct,ymax = mean_pct + sd_pct),
    width = 0.2
  ) +
  geom_jitter(
    data = df_line,
    aes(x = group, y = percentage, color = group),
    width = 0.15, size = 2, alpha = 0.8,
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(
    x = "Group",
    y = "Percentage of Simple repeat in Plasma CfRNA",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 8, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.line = element_line(color = "black"),  
    axis.text.x = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(0.5, 'cm'),   
    legend.key.height = unit(0.5, "cm"),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.width = unit(0.5, "cm"),    
    legend.box.spacing = unit(0.5, 'cm'),
    legend.key = element_rect(color = "black", size = 0.1),  
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  scale_y_continuous(expand = c(0, 0))+
   annotate("text", x = 2.5, y = y_max*1.1, label = "*",size=5) +
  expand_limits(y = y_max*1.3) + 
  annotate("segment", x = 2, xend = 3, y = y_max*1.08, yend = y_max*1.08,linewidth=0.3)  +
  annotate("text", x = 2, y = y_max*1.2, label = "***",size=5) +
  annotate("segment", x = 1, xend = 3, y = y_max*1.18, yend = y_max*1.18,linewidth=0.3)  

#Plasma_Simple_repeat_plot


df_line <- df_long %>% 
  filter(Biotype == "LTR")
df_line_summary <- df_line %>%
  group_by(group) %>%
  dplyr::summarise(
    mean_pct = mean(percentage),
    sd_pct = sd(percentage)
  )
custom_colors <- c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4")

y_max <- max(df_line$percentage)

Plasma_LTR_plot<-ggplot()+geom_bar(
  data = df_line_summary,
  aes(x = group, y = mean_pct, fill = group),
  stat = "identity", width = 0.6, color = "black",
  alpha = 0.6
) +
  geom_errorbar(
    data = df_line_summary,
    aes(x = group, ymin = mean_pct,ymax = mean_pct + sd_pct),
    width = 0.2
  ) +
  geom_jitter(
    data = df_line,
    aes(x = group, y = percentage, color = group),
    width = 0.15, size = 2, alpha = 0.8,
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(
    x = "Group",
    y = "Percentage of LTR in Plasma CfRNA",
    fill = "Group",
    color = "Group"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 8, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.line = element_line(color = "black"),  
    axis.text.x = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(0.5, 'cm'),   
    legend.key.height = unit(0.5, "cm"),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.width = unit(0.5, "cm"),     
    legend.box.spacing = unit(0.5, 'cm'),
    legend.key = element_rect(color = "black", size = 0.1),  
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  scale_y_continuous(expand = c(0, 0))+
  annotate("text", x = 2.5, y = y_max*1.1, label = "****",size=5) +
  expand_limits(y = y_max*1.3) + 
  annotate("segment", x = 2, xend = 3, y = y_max*1.08, yend = y_max*1.08,linewidth=0.3)  +
  annotate("text", x = 2, y = y_max*1.2, label = "***",size=5) +
  annotate("segment", x = 1, xend = 3, y = y_max*1.18, yend = y_max*1.18,linewidth=0.3)  
#Plasma_LTR_plot


combined_plot <- (
  urine_LINE_plot + urine_SINE_plot + urine_Simple_repeat_plot + urine_LTR_plot +
    Plasma_LINE_plot + Plasma_SINE_plot + Plasma_Simple_repeat_plot + Plasma_LTR_plot
) +
  plot_layout(ncol = 4, guides = "collect") & 
  theme(legend.position = "right")

combined_plot

ggsave("genetype_compare.png",width=10,height=6,dpi=1200, units = "in")
ggsave("genetype_compare.pdf",width=10,height=6,dpi=1200, units = "in", device = cairo_pdf)
