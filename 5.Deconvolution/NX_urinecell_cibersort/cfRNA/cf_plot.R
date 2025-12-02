library(pheatmap)
setwd("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/cell_type_marker/NX_urinecell_cibersort/cfRNA")
data<-read.table("CIBERSORTx_Adjusted.txt",row.names=1,header=T)
data<-data[grep("N",rownames(data)),]
#data<-data[-which(rownames(data)=='C107_N'),]

pheatmap(data[,1:14], cluster_rows=T, show_rownames=T,cluster_cols=T,show_colnames = T, 
         #annotation_col=coldata_all,
         #annotation_colors = ann_colors,na_col = "gray",
         scale = "none")

coldata<-matrix(ncol = 2,nrow = 106)
rownames(coldata)<-rownames(data)

idx_N_counts<-grep("N",rownames(data))
idx_S_counts<-grep("S",rownames(data))
idx_H_counts<-grep("H",rownames(data))
idx_C_counts<-grep("C",rownames(data))

coldata[idx_N_counts,1]<-"Urine"
coldata[idx_H_counts,2]<-"HC"
coldata[idx_S_counts,2]<-"BS"
coldata[idx_C_counts,2]<-"BC"


colnames(coldata)<-c("category","Group")

library(readxl)

BC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BC")
BC_cd<-as.data.frame(cbind(BC_clinicial_data$Patient,BC_clinicial_data$MIBC,BC_clinicial_data$AGE,BC_clinicial_data$SEX,
                           BC_clinicial_data$Tumor_size,BC_clinicial_data$`T`,BC_clinicial_data$Tumor_number,BC_clinicial_data$Grade))
colnames(BC_cd)<-c("Patients","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade')

BS_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BS")
BS_cd<-as.data.frame(cbind(BS_clinicial_data$Patient,NA,BS_clinicial_data$AGE,BS_clinicial_data$SEX,NA,NA,NA,NA))
colnames(BS_cd)<-c("Patients","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade')

HC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "HC")
HC_cd<-as.data.frame(cbind(HC_clinicial_data$Patient,NA,HC_clinicial_data$AGE,HC_clinicial_data$SEX,NA,NA,NA,NA))
colnames(HC_cd)<-c("Patients","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade')


clinicial_design<-rbind(BC_cd,BS_cd,HC_cd)

coldata<-as.data.frame(coldata)
coldata$sample_num<-rownames(coldata)
coldata$Patients<-rownames(coldata)
coldata$Patients <- sub(paste0("*", "_X"), "", coldata$Patients)
coldata$Patients <- sub(paste0("*", "_N"), "", coldata$Patients)


coldata_all<-merge(coldata,clinicial_design,by="Patients",all.x=T)
rownames(coldata_all)<-coldata_all$sample_num
coldata_all<-coldata_all[,-1]
coldata_all<-coldata_all[,-3]

t_coldata_all<-as.data.frame(t(coldata_all))
t_coldata_all <- t_coldata_all %>% dplyr::select(rownames(data))

coldata_all<-as.data.frame(t(t_coldata_all))
coldata_all$Tumor_size<-as.numeric(coldata_all$Tumor_size)
coldata_all$Age<-as.numeric(coldata_all$Age)

idx_2<-grep('2',coldata_all$StageT)
coldata_all$StageT[idx_2]<-'2'
idx_4<-grep('4',coldata_all$StageT)
coldata_all$StageT[idx_4]<-'4'

idx_single<-grep('1',coldata_all$Tumor_number)
coldata_all$Tumor_number[idx_single]<-'single'

idx_multiple<-grep('2',coldata_all$Tumor_number)
coldata_all$Tumor_number[idx_multiple]<-'multiple'


coldata_all<-coldata_all[,-1]
ann_colors <- list(
  Sex = c("Male" = "#406A93", "Female" = "#EFA9AE",'0'='white'),
  MIBC = c("MIBC" = "#65AE00", "NMIBC" = "#D3B356",'0'='white'),
  Group = c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4"),
  Tumor_size=colorRampPalette(c("white", "#2C1714"))(100),
  Age = colorRampPalette(c("#EEE0CC", "#91191E"))(100),
  StageT=c('4'="#FF3300",'3'='#FF6633','2'='#FF9966','1'='#FFCC66','a'='#FFFF00','0'='white'),
  Grade=c('Low'='#7A84D0','High'="#D78203",'0'='white'),
  sample=c('T'='#8E5D2A','NAT'="#64BE52",'0'='white'),
  Tumor_number=c('single'='#90D4A4','multiple'="#C49406",'0'='white')
)
coldata_all[is.na(coldata_all)] <- 0

data_heatmap<-as.matrix(t(data))
data_heatmap<-data_heatmap[1:14,]

#CairoPDF(file = "cfRNA_urine_cibersort.pdf", width = 16, height = 8,family='Arial')
#CairoPNG(file = "cfRNA_urine_cibersort.png", width = 16, height = 8,res=600,units='in',family='Arial')


pheatmap(data_heatmap,
         cluster_rows = TRUE,
         show_rownames = T,
         cluster_cols = T,
         show_colnames = F,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = coldata_all, 
         annotation_colors = ann_colors, 
         fontsize_row =10,
         #fontsize = 10,
         na_col = "gray",
         scale = "none",
         #color = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
         color = colorRampPalette(c("blue", "white", "red"))(100),
         legend=TRUE,
         name = 'Z-score')

#dev.off()

data_stack_pre<-as.data.frame(t(data_heatmap))
data_stack_pre<-cbind(data_stack_pre,coldata_all)
data_stack_S<-as.data.frame(apply(data_stack_pre[grep("S",rownames(data_stack_pre)),1:14],2,mean))
data_stack_S$Cell_Type<-rownames(data_stack_S)
colnames(data_stack_S)[1]<-'Value'
data_stack_C<-as.data.frame(apply(data_stack_pre[grep("C",rownames(data_stack_pre)),1:14],2,mean))
data_stack_C$Cell_Type<-rownames(data_stack_C)
colnames(data_stack_C)[1]<-'Value'
data_stack_H<-as.data.frame(apply(data_stack_pre[grep("H",rownames(data_stack_pre)),1:14],2,mean))
data_stack_H$Cell_Type<-rownames(data_stack_H)
colnames(data_stack_H)[1]<-'Value'

# 对数据按 Value 从大到小排序并计算百分比
data_stack_S <- data_stack_S %>%
  arrange(desc(Value)) %>%  # 按 Value 降序排序
  mutate(
    Percentage = Value / sum(Value) * 100,  # 计算百分比
    Cell_Type = factor(Cell_Type, levels = Cell_Type)  # 设置因子顺序为排序后的顺序
  )
gradient_colors <- c(
  "#E8AEB0",  # 红粉色
  "#F4C28B",  # 橙色
  "#E8DE7D",  # 浅黄色
  "#F5F58C",  # 柔和黄（稍微前移）
  "#D4C89A",  # 米色
  "#E6CFC2",  # 海贝壳色
  "#A9E0A9",  # 浅绿色
  "#C0DDE0",  # 淡青绿色（补充）
  "#A9DCE3",  # 浅青色
  "#A7CCE5",  # 淡蓝色
  "#C0DDEE",  # 蓝灰色
  "#B4B4DC",  # 淡紫色
  "#D0B7E1",  # 柔紫色
  "#F2B5D4",  # 柔粉紫色
  "#F5C0DA",  # 粉紫渐变（新增）
  "#EED4EA",  # 浅紫渐变（新增）
  "#F5E6F0",  # 淡紫粉（新增）
  "#CFCFCF",  # 浅灰色
  "#D9D9D9",  # 灰色过渡（新增）
  "#F8F7F0",  # 柔和米色（新增）
  "#FFF2CC",  # 浅黄色过渡（新增）
  "#FFE5B4"   # 柔橙色（新增）
)

# 绘制饼图
p1<-ggplot(data_stack_S, aes(x = "", y = Value, fill = Cell_Type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +  # 绘制条形图（圆环）
  coord_polar(theta = "y") +  # 转换为饼图
  labs(title = "", x = NULL, y = NULL) +
  scale_fill_manual(values = gradient_colors) + 
  theme_void() +  # 去掉背景和坐标轴
  geom_text(
    aes(label = ifelse(Percentage >= 4, paste0(round(Percentage, 1), "%"), "")), 
    position = position_stack(vjust = 0.5), size = 3
  )  # 添加大于等于 2% 的百分比标签
p1
dev.off()
#ggsave("BS_cell_type.png",height=6,width=6,dpi=600)


# 对数据按 Value 从大到小排序并计算百分比
data_stack_H <- data_stack_H %>%
  arrange(desc(Value)) %>%  # 按 Value 降序排序
  mutate(
    Percentage = Value / sum(Value) * 100,  # 计算百分比
    Cell_Type = factor(Cell_Type, levels = Cell_Type)  # 设置因子顺序为排序后的顺序
  )

# 绘制饼图
ggplot(data_stack_H, aes(x = "", y = Value, fill = Cell_Type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +  # 绘制条形图（圆环）
  coord_polar(theta = "y") +  # 转换为饼图
  labs(title = "Cell Type Proportions in Urine cfRNA in HC", x = NULL, y = NULL) +
  #scale_fill_manual(values = gradient_colors) +  # 使用渐变色
  theme_void() +  # 去掉背景和坐标轴
  geom_text(
    aes(label = ifelse(Percentage >= 4, paste0(round(Percentage, 1), "%"), "")), 
    position = position_stack(vjust = 0.5), size = 3
  )  # 添加大于等于 2% 的百分比标签
#ggsave("HC_cell_type.png",height=6,width=6,dpi=600)


# 对数据按 Value 从大到小排序并计算百分比
data_stack_C <- data_stack_C %>%
  arrange(desc(Value)) %>%  # 按 Value 降序排序
  mutate(
    Percentage = Value / sum(Value) * 100,  # 计算百分比
    Cell_Type = factor(Cell_Type, levels = Cell_Type)  # 设置因子顺序为排序后的顺序
  )

# 绘制饼图
ggplot(data_stack_C, aes(x = "", y = Value, fill = Cell_Type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +  # 绘制条形图（圆环）
  coord_polar(theta = "y") +  # 转换为饼图
  labs(title = "Cell Type Proportions in Urine cfRNA in HC", x = NULL, y = NULL) +
  #scale_fill_manual(values = gradient_colors) +  # 使用渐变色
  theme_void() +  # 去掉背景和坐标轴
  geom_text(
    aes(label = ifelse(Percentage >= 4, paste0(round(Percentage, 1), "%"), "")), 
    position = position_stack(vjust = 0.5), size = 3
  )  # 添加大于等于 2% 的百分比标签
#ggsave("HC_cell_type.png",height=6,width=6,dpi=600)

#ggsave("BC_cell_type.png",height=6,width=6,dpi=600)





gradient_colors <- c(
  "#E8AEB0",  # 红粉色
  "#F4C28B",  # 橙色
  "#E8DE7D",  # 浅黄色
  "#F5F58C",  # 柔和黄
  "#D4C89A",  # 米色
  "#E6CFC2",  # 海贝壳色
  "#A9E0A9",  # 浅绿色
  "#C0DDE0",  # 淡青绿色
  "#A9DCE3",  # 浅青色
  "#A7CCE5",  # 淡蓝色
  "#C0DDEE",  # 蓝灰色
  "#B4B4DC",  # 淡紫色
  "#D0B7E1",  # 柔紫色
  "#F2B5D4",  # 柔粉紫色
  "#F5C0DA",  # 粉紫渐变
  "#EED4EA",  # 浅紫渐变
  "#F5E6F0",  # 淡紫粉
  "#CFCFCF",  # 浅灰色
  "#D9D9D9",  # 灰色
  "#F8F7F0",  # 柔和米色
  "#FFF2CC",  # 浅黄色
  "#FFE5B4"   # 柔橙色
)

all_cell_types <- unique(c(
  as.character(data_stack_H$Cell_Type),
  as.character(data_stack_C$Cell_Type),
  as.character(data_stack_S$Cell_Type)
))
cell_type_colors <- setNames(gradient_colors[1:length(all_cell_types)], all_cell_types)
make_pie <- function(data, title = "") {
  data <- data %>%
    arrange(desc(Value)) %>%
    mutate(
      Percentage = Value / sum(Value) * 100,
      Cell_Type = factor(Cell_Type, levels = unique(Cell_Type))  # 按当前排序设定因子顺序
    )
  
  ggplot(data, aes(x = "", y = Value, fill = Cell_Type)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y", start = 0) +  # 从12点方向开始顺时针
    labs(title = title, x = NULL, y = NULL) +
    scale_fill_manual(values = cell_type_colors) +  # 使用统一颜色
    theme_void() +
    
    geom_text(
      aes(label = ifelse(Percentage >= 4, paste0(round(Percentage, 1), "%"), "")),
      position = position_stack(vjust = 0.5), size = 3
    )+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 10, margin = margin(t = 5)),
      axis.title.y = element_text(margin = margin(r = 0)),
      #axis.text.x = element_text(size = 10),
      plot.title = element_text(hjust = 0.5),
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
}
p1 <- make_pie(data_stack_H, "Healthy Controls")+ theme(legend.position = "none")
p2 <- make_pie(data_stack_C, "Bladder Cancer Patients")+ theme(legend.position = "none")
p3 <- make_pie(data_stack_S, "Bladder Stone Patients")

# 合并图，保留第一个图的图注
combined_plot <- p1 + p2 + p3 +plot_layout(ncol = 3)

combined_plot


ggsave("ALL_cell_type.png",height=4,width=10,units = 'in',dpi=600)
ggsave("ALL_cell_type.pdf",height=4,width=10,units = 'in', device = cairo_pdf)

# ====== Step 1. 整理数据为长格式 ======
data_stack_long <- data_stack_pre %>%
  dplyr::select(1:14, Group) %>%
  mutate(Sample = rownames(data_stack_pre)) %>%
  pivot_longer(cols = 1:14, names_to = "Cell_Type", values_to = "Value") %>%
  filter(!is.na(Group))

# ====== Step 2. 绘制带显著性标记的箱线图 ======
# 设置比较组
data_stack_long$Group <- factor(data_stack_long$Group, levels = c("HC","BC","BS"))
data_stack_long$Cell_Type <- factor(data_stack_long$Cell_Type)
my_comparisons <- list(c("HC", "BC"), c("HC", "BS"),c("BS","BC"))
library(ggpubr)

# 绘图
p <- ggplot(data_stack_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  #geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4")) +
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns = FALSE) +
  facet_wrap(~ Cell_Type, scales = "free_y") +  # 每个 Cell_Type 独立比较
  theme_bw(base_size = 13) +
  labs(x = "Group",
       y = "Estimated Proportion",
       title = "Urinary cfRNA Cell Type Composition across Groups") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank()
  )

print(p)

ggsave("comparing_cell_type.png",height=16,width=16,dpi=600,unit='in')
ggsave("comparing_cell_type.pdf",height=16,width=16)


selected_cells <- c("Macrophage", "Mast.Cell","Prostate.Epithelial","B","Endothelial","Neuronal")
data_subset <- data_stack_long %>% 
  filter(!Cell_Type %in% selected_cells)

# 确保 Group 是 factor
data_subset$Group <- factor(data_subset$Group, levels = c("HC","BC","BS"))

# 比较组
my_comparisons <- list(c("HC","BC"), c("HC","BS"),c("BS","BC"))

# 绘图
p <- ggplot(data_subset, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.7, width = 0.6) +
  #geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4")) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns = FALSE) +
  facet_wrap(~ Cell_Type, scales = "free_y",ncol = 4) +
  theme_bw(base_size = 13) +
  labs(x = "Group",
       y = "Estimated Proportion",
       #title = "Urinary cfRNA Cell Type Composition (Selected Types)"
       )+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank()
  )

print(p)

ggsave("comparing_cell_type.png",height=8,width=12,dpi=600,unit='in')
ggsave("comparing_cell_type.pdf",height=8,width=16)

str(data_stack_long)
table(data_stack_long$Group)
table(data_stack_long$Cell_Type)




# ====== Step 1. 整理数据为长格式 ======
data_stack_long <- data_stack_pre %>%
  dplyr::select(1:14, Sex) %>%
  mutate(Sample = rownames(data_stack_pre)) %>%
  pivot_longer(cols = 1:14, names_to = "Cell_Type", values_to = "Value") %>%
  filter(!is.na(Sex))

# ====== Step 2. 绘制带显著性标记的箱线图 ======
# 设置比较组
data_stack_long$Sex <- factor(data_stack_long$Sex, levels = c("Male","Female"))

#data_stack_long$Group <- factor(data_stack_long$Group, levels = c("HC","BC","BS"))
data_stack_long$Cell_Type <- factor(data_stack_long$Cell_Type)
my_comparisons <- list(c("Male", "Female"))
library(ggpubr)
data_stack_long <- data_stack_long %>% 
  filter(Cell_Type %in% c("Prostate.Epithelial"))
# 绘图
p <- ggplot(data_stack_long, aes(x = Sex, y = Value, fill = Sex)) +
  geom_boxplot(alpha = 0.7, width = 0.6) +
  #geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = c("Male" = "#406A93", "Female" = "#EFA9AE")) +
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns = FALSE) +
  facet_wrap(~ Cell_Type, scales = "free_y") +  # 每个 Cell_Type 独立比较
  theme_bw(base_size = 13) +
  labs(x = "Sex",
       y = "Estimated Proportion",
       title = "Urinary cfRNA Cell Type Composition across Groups") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank()
  )

print(p)
data_stack_long %>% 
  group_by(Cell_Type, Sex) %>% 
  summarise(mean_value = mean(Value, na.rm = TRUE),
            n = n()) %>% 
  print()
ggsave("sex_comparing_cell_type.png",height=6,width=6,dpi=600,unit='in')
ggsave("sex_comparing_cell_type.pdf",height=6,width=6)




# ====== Step 1. 整理数据为长格式 ======
data_stack_long <- data_stack_pre %>%
  dplyr::select(1:14, MIBC) %>%
  mutate(Sample = rownames(data_stack_pre)) %>%
  pivot_longer(cols = 1:14, names_to = "Cell_Type", values_to = "Value") %>%
  filter(MIBC %in% c("MIBC","NMIBC"))

# ====== Step 2. 绘制带显著性标记的箱线图 ======
# 设置比较组
data_stack_long$MIBC <- factor(data_stack_long$MIBC, levels = c("MIBC","NMIBC"))

#data_stack_long$Group <- factor(data_stack_long$Group, levels = c("HC","BC","BS"))
data_stack_long$Cell_Type <- factor(data_stack_long$Cell_Type)
my_comparisons <- list(c("MIBC", "NMIBC"))
library(ggpubr)
data_stack_long <- data_stack_long %>% 
  filter(Cell_Type %in% c("DC"))
#data_stack_long$Group <- factor(data_stack_long$Group, levels = c("HC","BC","BS"))

# 绘图
p <- ggplot(data_stack_long, aes(x = MIBC, y = Value, fill = MIBC)) +
  geom_boxplot(alpha = 0.7, width = 0.6) +
  #geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = c("MIBC" = '#65AE00', "NMIBC" = "#D3B356")) +
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns = FALSE) +
  facet_wrap(~ Cell_Type, scales = "free_y") +  # 每个 Cell_Type 独立比较
  theme_bw(base_size = 13) +
  labs(x = "Group",
       y = "Estimated Proportion",
       title = "Urinary cfRNA Cell Type Composition across Groups") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank()
  )

print(p)
ggsave("MIBC_comparing_cell_type.png",height=6,width=6,dpi=600,unit='in')
ggsave("MIBC_comparing_cell_type.pdf",height=6,width=6)
