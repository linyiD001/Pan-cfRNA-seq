#BiocManager::install("ggplot2")
library(stringr)
BC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BC")
BS_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BS")
HC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "HC")

clinic <- bind_rows(BC_clinicial_data,BS_clinicial_data,HC_clinicial_data)
clinic$group<-NA
clinic[grep('C',clinic$Patient),15]<-"BC"
clinic[grep('S',clinic$Patient),15]<-"BS"
clinic[grep('H',clinic$Patient),15]<-"HC"

filtered_namelist<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/QC/starstandard/QC_starstandard_list_withoutR.csv",header=T,row.names = 1)

N_data_filter<-filtered_namelist[grep("N",filtered_namelist[,1]),]
  
filtered_namelist_1<-t(as.data.frame(str_split(N_data_filter,"_")))[,1]

clinic<-clinic[which(clinic$Patient %in% filtered_namelist_1),]
ggplot(clinic, aes(x = group, y =  AGE, fill = group)) +
  geom_boxplot() +
  theme_light(base_size=10)+
  theme_classic()

custom_colors <- c("BC" = "#DB498E", "BS" = "#FCAE59", "HC" = "#379FB4")

# 绘图
p <- ggplot(clinic, aes(x = group, y = AGE, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_point(aes(color = group),
             position = position_dodge(width = 0.9),
             size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = custom_colors, name = "Group") +
  scale_color_manual(values = custom_colors, name = "Group") +
  annotate("text", x = 2, y = 110, label = "ns", size = 5) +
  theme_minimal() +
  labs(x = "Sample Type", y = "Age (year)") +
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


ggsave("N_age.png",height=4,width=4,dpi=1200, units = "in")
ggsave("N_age.pdf",height=4,width=4,dpi=1200, units = "in",device = cairo_pdf)

X_data_filter<-filtered_namelist[grep("X",filtered_namelist[,1]),]

filtered_namelist_2<-t(as.data.frame(str_split(X_data_filter,"_")))[,1]

clinic<-clinic[which(clinic$Patient %in% filtered_namelist_2),]
ggplot(clinic, aes(x = group, y =  AGE, fill = group)) +
  geom_boxplot() +
  theme_light(base_size=10)+
  theme_classic()

custom_colors <- c("BC" = "#DB498E", "BS" = "#FCAE59", "HC" = "#379FB4")

# 绘图
p <- ggplot(clinic, aes(x = group, y = AGE, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_point(aes(color = group),
             position = position_dodge(width = 0.9),
             size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = custom_colors, name = "Group") +
  scale_color_manual(values = custom_colors, name = "Group") +
  annotate("text", x = 2, y = 110, label = "ns", size = 5) +
  theme_minimal() +
  labs(x = "Sample Type", y = "Age (year)") +
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


ggsave("X_age.png",height=4,width=4,dpi=1200, units = "in")
ggsave("X_age.pdf",height=4,width=4,dpi=1200, units = "in",device = cairo_pdf)