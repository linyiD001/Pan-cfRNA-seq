library(tidyr)
library(dplyr)
library(ggplot2)

countdata<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome/cfRNA_repeataware_count.csv")
all_genelist_combine<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome/all_genelist_annotation.csv",header=T,row.names=1)

countdata$enst<-countdata$X
countdata<-countdata[,-1]
countdata_merge<-merge(countdata,all_genelist_combine,by='enst',all.x=T)

LINE_count<-countdata_merge[which(countdata_merge$genetype=='LINE'),-grep("X",colnames(countdata_merge))]

LINE_count_symbol<- aggregate(LINE_count[,2:143],by=list(LINE_count$gene_symbol), FUN=sum)
LINE_count_symbol$new_family<-NA
LINE_count_symbol$new_family[grep("CR1",LINE_count_symbol$Group.1)] <- "CR1"
LINE_count_symbol$new_family[grep("L1",LINE_count_symbol$Group.1)] <- "L1"
LINE_count_symbol$new_family[grep("L2",LINE_count_symbol$Group.1)] <- "L2"
LINE_count_symbol$new_family[grep("L3",LINE_count_symbol$Group.1)] <- "L3"
LINE_count_symbol_new_family<- aggregate(LINE_count_symbol[,2:143],by=list(LINE_count_symbol$new_family), FUN=sum)

LINE_count_symbol_new_family_t<-as.data.frame(t(LINE_count_symbol_new_family))
  
colnames(LINE_count_symbol_new_family_t)<-LINE_count_symbol_new_family_t[1,]
LINE_count_symbol_new_family_t<-LINE_count_symbol_new_family_t[-1,]
rownames_backup <- rownames(LINE_count_symbol_new_family_t)
LINE_count_symbol_new_family_t[] <- lapply(LINE_count_symbol_new_family_t, as.numeric)
rownames(LINE_count_symbol_new_family_t) <- rownames_backup

LINE_count_symbol_new_family_t$total_count <-colSums(countdata_merge[,grep("N",colnames(countdata_merge))])
LINE_ratio_symbol <- LINE_count_symbol_new_family_t
rownames(LINE_ratio_symbol)<-rownames(LINE_count_symbol_new_family_t)
LINE_ratio_symbol[, 1:4] <- LINE_count_symbol_new_family_t[, 1:4]/ LINE_count_symbol_new_family_t$total_count
LINE_ratio_symbol$group<-NA
LINE_ratio_symbol$group[grep("C",rownames(LINE_ratio_symbol))] <- "BC"
LINE_ratio_symbol$group[grep("H",rownames(LINE_ratio_symbol))] <- "HC"
LINE_ratio_symbol$group[grep("S",rownames(LINE_ratio_symbol))] <- "BS"

QC_list<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/QC_repeataware_list.csv",row.names = 1)

LINE_ratio_symbol<-LINE_ratio_symbol[rownames(LINE_ratio_symbol) %in% QC_list$x,]
idx_R <- which(rownames(LINE_ratio_symbol) %in% c("C107_N"))
LINE_ratio_symbol<-LINE_ratio_symbol[-idx_R,]
  

ratio_data <- LINE_ratio_symbol[, c("CR1", "L1", "L2", "L3", "group")]
ratio_long <- pivot_longer(ratio_data, cols = c("CR1", "L1", "L2", "L3"),
                           names_to = "Gene", values_to = "Ratio")

group_avg <- ratio_long %>%
  group_by(group, Gene) %>%
  summarise(AvgRatio = mean(Ratio, na.rm = TRUE), .groups = "drop")

ggplot(group_avg, aes(x = reorder(Gene, -AvgRatio), y = AvgRatio, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "LINE Family",
       y = "Average Ratio") +
  theme_light(base_size = 10)+
  scale_fill_manual(values = c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave("LINE_proportion_in_N.png",width = 5, height = 3, units = "in", dpi = 300)
ggsave("LINE_proportion_in_N.pdf",width = 5, height = 3, units = "in", dpi = 300)
