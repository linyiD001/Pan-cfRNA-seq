####做一个血尿配对的？

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






# 只保留血尿配对样本
paired_samples <- RNA_CON_cf %>%
  group_by(NUMBER) %>%
  filter(all(c("Plasma", "Urine") %in% CATE)) %>%
  ungroup()



library(dplyr)

# 为配对绘图准备数据
paired_plot <- paired_samples %>%
  dplyr::select(NUMBER, CATE, Con, Group) %>%
  spread(key = CATE, value = Con) %>%
  drop_na() # 确保血尿都有值

# 配对统计：血 vs 尿
paired_test <- wilcox.test(paired_plot$Plasma, paired_plot$Urine, paired = TRUE)

# 绘图
paired_long <- paired_plot %>%
  pivot_longer(cols = c("Plasma","Urine"), names_to = "CATE", values_to = "Con")

ggplot(paired_long, aes(x = CATE, y = Con, group = NUMBER)) +
  geom_point(aes(color = Group), size = 3) +
  geom_line(aes(group = NUMBER, color = Group), alpha = 0.5) +
  #geom_boxplot(aes(fill = CATE), alpha = 0.3, width = 0.4, outlier.shape = NA) +
  scale_fill_manual(values = wes_palette("Darjeeling1")[1:2]) +
  scale_color_manual(values = wes_palette("Darjeeling1")[3:5]) +
  theme_classic(base_size = 14) +
  labs(title = "Paired Plasma vs Urine RNA Concentration",
       y = "RNA Concentration",
       x = "") +
  geom_signif(comparisons = list(c("Plasma", "Urine")),
              test = "wilcox.test", paired = TRUE, map_signif_level = TRUE)
# 假设 paired_plot 已经生成（每个个体都有 Plasma 和 Urine）

paired_plot_clean <- paired_plot %>% 
  filter(!is.na(Plasma) & !is.na(Urine))

# 全部样本的 Pearson 相关系数
cor_all <- cor.test(log10(paired_plot_clean$Plasma), log10(paired_plot_clean$Urine), method = "spearman")
R_all <- cor_all$estimate
P_all <- cor_all$p.value


# 绘图
ggplot(paired_plot_clean, aes(x = Plasma, y = Urine, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("BC"="#DB498E", "BS"="#FCAE59", "HC"="#379FB4")) +
  scale_x_continuous(trans = "log10", labels = scales::label_number()) +
  scale_y_continuous(trans = "log10", labels = scales::label_number()) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  #scale_color_manual(values = wes_palette("Darjeeling1")[3:5]) +
  theme_classic(base_size = 14) +
  labs(title = "Plasma vs Urine RNA Concentration",
       x = "Plasma RNA Concentration (log10 ng/ml)",
       y = "Urine RNA Concentration (log10 ng/ml)",
       color = "Group") +
  annotate("text",
           x = max((paired_plot_clean$Plasma))*0.2,
           y = max((paired_plot_clean$Urine))*0.6,
           label = paste0("R = ", round(R_all, 2), "\nP = ", signif(P_all, 3)),
           size = 4) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title = element_text(size = 15, family = 'Arial',face = "bold")
  )
# 保存图
ggsave("Plasma_vs_Urine_log10_Scatter_with_RP.png", width = 6, height = 6, dpi = 600)
ggsave("Plasma_vs_Urine_log10_Scatter_with_RP.pdf", width = 6, height = 6, units = 'in', device = cairo_pdf)






aggregate(Con ~ NAME + CATE, data = RNA_CON_cf, FUN = median)

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
  geom_text(aes(x = 1.5, y = 1.62^10, label = "***"), size = 4,family = "")  
q    
ggsave("concentration_of_cfRNA.png",plot = q, height=4,width=4,dpi=1200, units = "in")
ggsave("concentration_of_cfRNA.pdf", plot = q, height=4,width=4,dpi=1200, units = "in")


