
setwd("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/MAPPING/coverage")

rm(list = ls())

cfRNA_name<- read.table("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/MAPPING/cfRNA.txt",header = F)
cfRNA_name<-cfRNA_name[-1,]
cfRNA_name<-as.data.frame(cfRNA_name)


i="C100_N"
des<-paste("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/MAPPING/coverage/",i,"_coverage.geneBodyCoverage.txt",sep="")
data_tem<-read.table(des,header=T)
rownames(data_tem)<-i
coverage<-data_tem

for (i in cfRNA_name[1:275,]){
  des<-paste("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/MAPPING/coverage/",i,"_coverage.geneBodyCoverage.txt",sep="")
  data_tem<-read.table(des,header=T)
  rownames(data_tem)<-i
  coverage<-rbind(coverage,data_tem)
  print(i)
}

coverage<-coverage[,-1]

write.csv(coverage,"coverage.csv")


coverage<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/MAPPING/coverage/coverage.csv",row.names = 1)

##filtered
filtered_namelist<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/QC/starstandard/QC_starstandard_list_withoutR.csv",header=T,row.names = 1)
idx_filtered_namelist<-which(rownames(coverage) %in% filtered_namelist[,1])
coverage<-coverage[idx_filtered_namelist,]


idx_X<-grep("X",rownames(coverage))
idx_N<-grep("N",rownames(coverage))

data_for_plot<-coverage[idx_N,]
a<-as.data.frame(t(data_for_plot))
rownames(a)<-1:100

processed_data <-as.data.frame(apply(a, 2, function(x) ((x - min(x))/(max(x) - min(x)))))

d<-t(processed_data)

average_d<-apply(d,2,mean)
sd_d<-apply(d,2,sd)
min_d<-apply(d,2,min)
max_d<-apply(d,2,max)

#e<-1:length
e<-1:100
f<-rbind(average_d,sd_d,e)
f<-t(f)
f<-as.data.frame(f)
f$e<-as.numeric(f$e)

which(f$average_d>0.0008)

# 绘制密度分布曲线和误差范围
plot_1<-ggplot(f,aes(x=e, y=average_d,)) + geom_errorbar(aes(ymin=average_d-sd_d, ymax=average_d+sd_d),color="grey",alpha = 0.2)+  
  theme(axis.ticks=element_blank())+xlab("position")+
  theme_minimal()+geom_line(color = "black") +  # 添加密度分布曲线
  geom_smooth(aes(group = 1), method = "loess", se = F, color = "darkgreen") +  # 添加平滑曲线和误差范围
  theme_minimal()+
  #geom_hline(yintercept = average_d, linetype = "dashed", color = "red") +  # 添加均值水平线
  #geom_segment(aes(x = 1, y = average_d - sd_d, xend = 1000, yend = average_d + sd_d), color = "red") +  # 添加误差范围线
  labs(title = "Urine_Gene_coverage", x = "Gene body percentile(5'->3')", y = "Coverage")

plot_1

ggsave("coverage_N_filtered.png", width = 4, height = 4, dpi = 300,units='in')
ggsave("coverage_N_filtered.pdf", width = 4, height = 4, dpi = 300,units='in')



data_for_plot<-coverage[idx_X,]
a<-as.data.frame(t(data_for_plot))
rownames(a)<-1:100

processed_data <-as.data.frame(apply(a, 2, function(x) ((x - min(x))/(max(x) - min(x)))))

d<-t(processed_data)

average_d<-apply(d,2,mean)
sd_d<-apply(d,2,sd)
min_d<-apply(d,2,min)
max_d<-apply(d,2,max)

#e<-1:length
e<-1:100
f<-rbind(average_d,sd_d,e)
f<-t(f)
f<-as.data.frame(f)
f$e<-as.numeric(f$e)

which(f$average_d>0.0008)

# 绘制密度分布曲线和误差范围
plot_1<-ggplot(f,aes(x=e, y=average_d,)) + geom_errorbar(aes(ymin=average_d-sd_d, ymax=average_d+sd_d),color="grey",alpha = 0.2)+  
  theme(axis.ticks=element_blank())+xlab("position")+
  theme_minimal()+geom_line(color = "black") +  # 添加密度分布曲线
  geom_smooth(aes(group = 1), method = "loess", se = F, color = "darkgreen") +  # 添加平滑曲线和误差范围
  theme_minimal()+
  #geom_hline(yintercept = average_d, linetype = "dashed", color = "red") +  # 添加均值水平线
  #geom_segment(aes(x = 1, y = average_d - sd_d, xend = 1000, yend = average_d + sd_d), color = "red") +  # 添加误差范围线
  labs(title = "Plasma_Gene_coverage", x = "Gene body percentile(5'->3')", y = "Coverage")

plot_1

ggsave("coverage_X_filtered.png", width = 4, height = 4, dpi = 300, units='in')
ggsave("coverage_X_filtered.pdf", width = 4, height = 4, dpi = 300, units='in')
