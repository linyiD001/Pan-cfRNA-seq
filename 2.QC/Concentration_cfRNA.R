setwd("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/QC/starstandard")
rm(list = ls())
library(readxl)
library(ggplot2)
library(ggsignif)
library(readxl)
library(ggplot2)
library(ggsignif)
#install.packages("wesanderson")
library(wesanderson)
library(ggbreak)

# Load
RNA_CON_cf <- as.data.frame(read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/QC/con/RNA_CON_cf.xlsx"))

colnames(RNA_CON_cf)[8]<-"Con"

#### 排除QC 不合格的
QC_list<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/QC/starstandard/QC_starstandard_list_withoutR.csv",row.names = 1)
RNA_CON_cf<-RNA_CON_cf[which(RNA_CON_cf$NAME %in% QC_list[,1]),]


idx_S_counts<-grep("S",RNA_CON_cf$NAME)
idx_H_counts<-grep("H",RNA_CON_cf$NAME)
idx_C_counts<-grep("C",RNA_CON_cf$NAME)

RNA_CON_cf$Group[idx_H_counts]<-"HC"
RNA_CON_cf$Group[idx_S_counts]<-"BS"
RNA_CON_cf$Group[idx_C_counts]<-"BC"

idx_X<-grep("X",RNA_CON_cf$CATE)
RNA_CON_cf$CATE[idx_X]<-"Plasma"
idx_N<-grep("N",RNA_CON_cf$CATE)
RNA_CON_cf$CATE[idx_N]<-"Urine"


aggregate(Con ~ Group + CATE, data = RNA_CON_cf, FUN = median)

RNA_CON_cf$AC<-paste0(RNA_CON_cf$Group,RNA_CON_cf$CATE)

summary_stats <- aggregate( Con ~ Group + CATE, data = RNA_CON_cf, 
                            FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))

#SX-CN SX-CX SX-HN SX-SN 
ggplot(data = summary_stats, aes(x = CATE, y = Con[,1], fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Con[,1] - Con[,2], ymax = Con[,1] + Con[,2]), 
                width = 0.2, position = position_dodge(width = 0.9)) +
  labs( x = "Sample Type", y = "Con of cfRNA / (ng/ml)", fill = "patient") + 
  #scale_fill_discrete_qualitative(palette = "Dark 2") + 
  #geom_signif(y_position = c(7,7.5), xmin = c(0.65,0.88), 
  #            xmax = c(0.87,1.9), annotation = c("***","***"),
  #            tip_length = 0)+
  theme_light(base_size=10)+
  theme_classic()


# Create the violin plot
# 自定义颜色
custom_colors <- c("HC" = "#379FB4", "BS" = "#FCAE59", "BC" = "#DB498E")

# 画图
p <- ggplot(data = RNA_CON_cf, aes(x = CATE, y = Con, fill = Group)) +
  geom_boxplot(outliers = F, alpha = 0.6,position = position_dodge(width = 0.9)) +
  geom_jitter(aes(color = Group),
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9),
              size = 1, alpha = 0.7)+
  labs(x = "Sample Type", y = "cfRNA (ng/ml)") +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  scale_y_log10() +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    axis.line = element_line(color = "black"),  # 显示 xy 轴线
    axis.text.x = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   # 设置图例项的水平间距
    legend.spacing.y = unit(0.5, 'cm'),   # 设置图例项的垂直间距
    legend.key.height = unit(0.5, "cm"),    # 设置图例色块的高度
    legend.key.spacing.y = unit(0.1, "cm"),    # 设置图例色块的高度
    legend.key.width = unit(0.5, "cm"),     # 设置图例色块的宽度
    legend.box.spacing = unit(0.5, 'cm'),
    #panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    # 图例背景和边框
    #legend.background = element_rect(fill = "white", color = "black", size = 0.5), # 图例背景和边框
    legend.key = element_rect(color = "black", size = 0.1)  # 每个图例项的边框
  )
p

kruskal_result <- kruskal.test(Con ~ AC, data = RNA_CON_cf)
kruskal_result

# 加载 dunn.test 包
library(dunn.test)

# Dunn’s test
dunn_result <- dunn.test(RNA_CON_cf$Con, RNA_CON_cf$AC, method = "bonferroni")
dunn_result<-as.data.frame(dunn_result)
significant_results <- dunn_result$comparisons[which(dunn_result$P.adjusted < 0.05) ]
print(significant_results)

# 添加统计学检验并在图上标注显著性
q <- p + 
  #BSPlasma-BCPlasma
  geom_segment(aes(x = 0.7, xend = 1, y = 1.45^10, yend = 1.45^10),linewidth=0.3) +  
  geom_text(aes(x = 0.85, y = 1.46^10, label = "****"), size = 4,family = "") +
  #BCPlasma-BSPlasma
  geom_segment(aes(x = 0.7, xend = 1.3,y = 1.49^10, yend = 1.49^10), linewidth=0.3) +  
  geom_text(aes(x = 1, y = 1.50^10, label = "***"), size = 4,family = "") +
  #BCUrine-BSUrine
  geom_segment(aes(x = 1.7, xend = 2, y = 1.45^10, yend = 1.45^10),linewidth=0.3) + 
  geom_text(aes(x = 1.85, y = 1.46^10, label = "***"), size = 4,family = "") +
  #BCPlasma-BCUrine
  geom_segment(aes(x = 0.7, xend = 1.7, y = 1.53^10, yend = 1.53^10),linewidth=0.3) +  
  geom_text(aes(x = 1.2, y = 1.54^10, label = "****"), size = 4,family = "") +   
  #BSPlasma-BSUrine
  geom_segment(aes(x = 1, xend = 2, y = 1.57^10, yend = 1.57^10),linewidth=0.3) +  
  geom_text(aes(x = 1.5, y = 1.58^10, label = "***"), size = 4,family = "") +
  #HCPlasma-HCUrine
  geom_segment(aes(x = 1.3, xend = 2.3, y = 1.61^10, yend = 1.61^10),linewidth=0.3) +  
  geom_text(aes(x = 1.8, y = 1.62^10, label = "***"), size = 4,family = "")  
q    
ggsave("concentration_of_cfRNA.png",plot = q, height=4,width=4,dpi=1200, units = "in")
ggsave("concentration_of_cfRNA.pdf", plot = q, height=4,width=4,dpi=1200, units = "in")


