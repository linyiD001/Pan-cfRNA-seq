####correlation between TEcount & repeataware
library(dplyr) 
library(tidyr)
library(stats)
library(ggplot2)
library(ggpubr)
library(scales)
library(readr)

TEcount_data<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/TEcount/cf_Transcript_genename.csv",header=T,row.names = 1)
TEcount_data<-aggregate(TEcount_data[,2:277],by=list(TEcount_data$name), FUN=sum)

QC_list<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/QC_repeataware_list.csv",row.names = 1)

repeataware_data<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome/salmon.repeataware_simple_genename.csv",header=T,row.names = 1)

list<-intersect(TEcount_data$Group.1,repeataware_data$Group.1)

common_genes <- intersect(TEcount_data$Group.1, repeataware_data$Group.1)
TEcount_data_common <- TEcount_data[which(TEcount_data$Group.1 %in% common_genes), colnames(TEcount_data) %in% QC_list[,1]]
repeataware_data_common <- repeataware_data[which(repeataware_data$Group.1 %in% common_genes), colnames(repeataware_data) %in% QC_list[,1]]

rownames(TEcount_data_common) <- TEcount_data$Group.1[which(TEcount_data$Group.1 %in% common_genes)]
rownames(repeataware_data_common) <- repeataware_data$Group.1[which(repeataware_data$Group.1 %in% common_genes)]


total_reads_TEcount <- colSums(TEcount_data_common)  
RPM_TEcount <- sweep(TEcount_data_common, 2, total_reads_TEcount, FUN = "/") * 1e6  

total_reads_repeataware <- colSums(repeataware_data_common)  
RPM_repeataware <- sweep(repeataware_data_common, 2, total_reads_repeataware, FUN = "/") * 1e6  

log2_RPM_TEcount<-log2(RPM_TEcount+1)
log2_RPM_repeataware<-log2(RPM_repeataware+1)
correlation_results_new <- data.frame(Column = colnames(TEcount_data_common), Correlation = NA, P_Value = NA)

for (i in 1:ncol(TEcount_data_common)) {
  test_result <- cor.test(log2_RPM_TEcount[, i], log2_RPM_repeataware[, i], use = "complete.obs", method='spearman',exact = FALSE)
  correlation_results_new$Correlation[i] <- test_result$estimate  
  correlation_results_new$P_Value[i] <- test_result$p.value       
}

write.csv(correlation_results_new,"correlation_results_new_log2RPM.csv")

###HC_urine
RPM_data <- data.frame(
  TEcount_RPM = as.vector(as.matrix(RPM_TEcount$H26_N)),
  Repeataware_RPM = as.vector(as.matrix(RPM_repeataware$H26_N))
)
cor_test <- cor.test(log2(RPM_data$TEcount_RPM + 1), log2(RPM_data$Repeataware_RPM + 1), method = "pearson")
cor_value <- round(cor_test$estimate, 3)
p_value <- cor_test$p.value

ggplot(RPM_data, aes(x = log2(TEcount_RPM+1), y = log2(Repeataware_RPM+1))) +
  geom_point(alpha = 0.3) +                       
  labs(
    title = "Comparison of TEcount and Repeataware RPM Values",
    x = "log2(TEcount RPM+1)",
    y = "log2(Repeataware RPM+1)"
  ) +
  theme_minimal()+
  annotate("text", 
           x = 1, 
           y = max(log2(RPM_data$Repeataware_RPM + 1)) * 0.9,
           label = paste("Pearson r =", cor_value, "\nP-value =", format.pval(p_value, digits = 3)),
           size = 3,
           hjust = 0)
ggsave("TE_H26_N.png",height=4,width=4,dpi=600,units='in')
ggsave("TE_H26_N.pdf",height=4,width=4,dpi=600,units='in')



ggplot(correlation_results_new, aes(x =Correlation)) + 
  geom_density(alpha = 0.5) +
  labs(title = "",
       x = "Correlation between two methods",
       y = "Density") +
  theme_minimal()
ggsave("cor_density_plot.png",height=12,width=12,dpi=600)
