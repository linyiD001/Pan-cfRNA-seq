###OUR DATA
repeataware_mapping_rate<-read.table("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome/multiqc_data/multiqc_salmon.txt",sep='\t',header=T)
standard_mapping_rate<-read.table("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.standard/multiqc_data/multiqc_salmon.txt",sep='\t',header=T)

data_merge<-merge(repeataware_mapping_rate[,c(1,35)],standard_mapping_rate[,c(1,35)],by='Sample')
data_merge<-data_merge[grep("N",data_merge$Sample),]

QC_list<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/QC_repeataware_list.csv",row.names = 1)
data_merge<-data_merge[which(data_merge$Sample %in% QC_list$x),]

colnames(data_merge) <- c("Sample", "RepeatAware", "Standard")
library(reshape2)
library(ggplot2)
library(ggpubr)
library(dplyr)

data_long <- melt(data_merge, id.vars = "Sample", variable.name = "Method", value.name = "MappingRate")

data_long$MappingRate<-as.numeric(data_long$MappingRate)


table_counts <- data_long %>%
  group_by(Sample) %>%
  summarise(n = n())

print(table_counts)

samples_complete <- table_counts %>%
  filter(n == 2) %>%
  pull(Sample)

data_long_complete <- data_long %>%
  filter(Sample %in% samples_complete)


ggplot(data_long_complete, aes(x = Method, y = MappingRate)) +
  geom_boxplot(aes(fill = Method), alpha = 0.3, outlier.shape = NA, width = 0.6) +  
  geom_point(aes(group = Sample), color = 'black', size = 2) + 
  geom_line(aes(group = Sample), color = 'grey', alpha = 0.5) +  
  ylab("Mapping Rate (%)") +
  theme_minimal() + stat_compare_means(method = "t.test",paired = T,label = "p.signif",
                                       label.y = max(data_long_complete$MappingRate) + 5,
                                       label.x = 1.5,size=5)

ggsave("standard_vs_repeataware.png",height=6,width=6,dpi=300)
ggsave("standard_vs_repeataware.pdf",height=6,width=6,dpi=300)


###PUBLIC DATA
###

sample_info<-read.csv("/dssg/home/acct-dahan/share/cfRNA/published/quake/SraRunTable.csv",header=T)
repeataware_mapping_rate<-read.table("/dssg/home/acct-dahan/share/cfRNA/published/quake/salmon_repeataware/multiqc_data/multiqc_salmon.txt",sep='\t',header=T)
standard_mapping_rate<-read.table("/dssg/home/acct-dahan/share/cfRNA/published/quake/salmon_standard/multiqc_data/multiqc_salmon.txt",sep='\t',header=T)

data_merge<-merge(repeataware_mapping_rate[,c(1,35)],standard_mapping_rate[,c(1,35)],by='Sample')
colnames(sample_info)[1]<-'Sample'
data_merge<-merge(data_merge,sample_info[,c(1,34)],by='Sample')
data_merge<-data_merge[which(data_merge$sample_type=='supernatant'),]
data_merge<-data_merge[,1:3]
colnames(data_merge) <- c("Sample", "RepeatAware", "Standard")
library(reshape2)
library(ggplot2)
library(ggpubr)
library(dplyr)

data_long <- melt(data_merge, id.vars = "Sample", variable.name = "Method", value.name = "MappingRate")

data_long$MappingRate<-as.numeric(data_long$MappingRate)

table_counts <- data_long %>%
  group_by(Sample) %>%
  summarise(n = n())

print(table_counts)

samples_complete <- table_counts %>%
  filter(n == 2) %>%
  pull(Sample)

data_long_complete <- data_long %>%
  filter(Sample %in% samples_complete)


ggplot(data_long_complete, aes(x = Method, y = MappingRate)) +
  geom_boxplot(aes(fill = Method), alpha = 0.3, outlier.shape = NA, width = 0.6) +  
  geom_point(aes(group = Sample), color = 'black', size = 2) +  
  geom_line(aes(group = Sample), color = 'grey', alpha = 0.5) +  
  ylab("Mapping Rate (%)") +
  theme_minimal() + stat_compare_means(method = "t.test",paired = T,label = "p.signif",
                                       label.y = max(data_long_complete$MappingRate) + 5,
                                       label.x = 1.5,size=5)


ggsave("quake_standard_vs_repeataware.png",height=6,width=6,dpi=300)
ggsave("quake_standard_vs_repeataware.pdf",height=6,width=6,dpi=300)
