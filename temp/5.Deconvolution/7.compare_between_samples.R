###cf_compare

data<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/cell_type_marker/NX_urinecell_cibersort/cf_cibersort_urinecell/CIBERSORTx_Job2_Results.csv",row.names=1)
data<-data[grep("N",rownames(data)),]
data<-data[-which(rownames(data)=='C107_N'),]

data<-data[,-c(15:17)]
data$immuno<-data$B+data$CD4T+data$CD8T+data$DC+data$Macrophage+data$Mono+data$Mast.Cell+data$Neutrophil+data$NK


library(readxl)

BC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BC")
BC_cd<-as.data.frame(cbind(BC_clinicial_data$Patient,BC_clinicial_data$MIBC,BC_clinicial_data$AGE,BC_clinicial_data$SEX,
                           BC_clinicial_data$Tumor_size,BC_clinicial_data$`T`,BC_clinicial_data$Tumor_number,BC_clinicial_data$Grade))
colnames(BC_cd)<-c("Patients","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade')

BS_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BS")
BS_cd<-as.data.frame(cbind(BS_clinicial_data$Patient,NA,BS_clinicial_data$AGE,BS_clinicial_data$SEX,NA,NA,NA,NA))
colnames(BS_cd)<-c("Patients","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade')

HC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "HC")
HC_cd<-as.data.frame(cbind(HC_clinicial_data$Patient,NA,HC_clinicial_data$AGE,HC_clinicial_data$SEX,NA,NA,NA,NA))
colnames(HC_cd)<-c("Patients","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade')


clinicial_design<-rbind(BC_cd,BS_cd,HC_cd)

idx_S_counts<-grep("S",clinicial_design$Patients)
idx_H_counts<-grep("H",clinicial_design$Patients)
idx_C_counts<-grep("C",clinicial_design$Patients)

clinicial_design$Group[idx_H_counts]<-"HC"
clinicial_design$Group[idx_S_counts]<-"BS"
clinicial_design$Group[idx_C_counts]<-"BC"

df_long <-data %>%
  rownames_to_column(var = "sample") %>% 
  pivot_longer(cols = -sample,           
               names_to = "celltype",
               values_to = "percentage")

df_long$Patients<- sub(paste0("*", "_N"), "", df_long$sample)


data_merge<-merge(df_long,clinicial_design,by='Patients')


df_line <- data_merge %>% 
  filter(celltype == "Prostate.Epithelial")
df_line_summary <- df_line %>%
  group_by(Sex) %>%
  dplyr::summarise(
    mean_pct = mean(percentage),
    sd_pct = sd(percentage)
  )
custom_colors <- c("Male" = "#406A93", "Female" = "#EFA9AE")

y_max <- max(df_line$percentage)

Sex_Prostate_plot<-ggplot()+geom_bar(
  data = df_line_summary,
  aes(x = Sex, y = mean_pct, fill = Sex),
  stat = "identity", width = 0.6, color = "black",
  alpha = 0.6
) +
  geom_errorbar(
    data = df_line_summary,
    aes(x = Sex, ymin = mean_pct,ymax = mean_pct + sd_pct),
    width = 0.2
  ) +
  geom_jitter(
    data = df_line,
    aes(x = Sex, y = percentage, color = Sex),
    width = 0.15, size = 2, alpha = 0.8,
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(
    #title = "Mean LINE Percentage by Group",
    x = "Sex",
    y = "Prostate Epithelial Proportion in Urine CfRNA"
    #fill = "Sex",
    #color = "Sex"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
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
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.position = "none"  
  )+
  scale_y_continuous(expand = c(0, 0))+
  annotate("text", x = 1.5, y = y_max*1.2, label = "****",size=5) +
  expand_limits(y = y_max*1.3) +
  annotate("segment", x = 1, xend = 2, y = y_max*1.18, yend = y_max*1.18,linewidth=0.3)  

Sex_Prostate_plot


df_line <- data_merge %>% 
  filter(celltype == "immuno")
df_line_summary <- df_line %>%
  group_by(Group) %>%
  dplyr::summarise(
    mean_pct = mean(percentage),
    sd_pct = sd(percentage)
  )
custom_colors <- c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4")

y_max <- max(df_line$percentage)

Immuno_Group_plot<-ggplot()+geom_bar(
  data = df_line_summary,
  aes(x = Group, y = mean_pct, fill = Group),
  stat = "identity", width = 0.6, color = "black",
  alpha = 0.6
) +
  geom_errorbar(
    data = df_line_summary,
    aes(x = Group, ymin = mean_pct,ymax = mean_pct + sd_pct),
    width = 0.2
  ) +
  geom_jitter(
    data = df_line,
    aes(x = Group, y = percentage, color = Group),
    width = 0.15, size = 2, alpha = 0.8,
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(
    #title = "Mean LINE Percentage by Group",
    x = "Group",
    y = "Immune Proportion in Urine CfRNA"
    #fill = "Sex",
    #color = "Sex"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
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
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.position = "none"  
  )+
  scale_y_continuous(expand = c(0, 0))+
  expand_limits(y = y_max*1.3) + 
  annotate("text", x = 2, y = y_max*1.2, label = "*",size=5) +
  annotate("segment", x = 1, xend = 3, y = y_max*1.18, yend = y_max*1.18,linewidth=0.3)+
  annotate("text", x = 2.5, y = y_max*1.12, label = "****",size=5) +
  annotate("segment", x = 2, xend = 3, y = y_max*1.1, yend = y_max*1.1,linewidth=0.3)+
  annotate("text", x = 1.5, y = y_max*1.08, label = "*",size=5) +
  annotate("segment", x = 1, xend = 2, y = y_max*1.06, yend = y_max*1.06,linewidth=0.3)  

Immuno_Group_plot



df_line <- data_merge %>% 
  filter(celltype == "DC")
df_line_summary <- df_line %>%
  group_by(Group) %>%
  dplyr::summarise(
    mean_pct = mean(percentage),
    sd_pct = sd(percentage)
  )
custom_colors <- c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4")

y_max <- max(df_line$percentage)

DC_Group_plot<-ggplot()+geom_bar(
  data = df_line_summary,
  aes(x = Group, y = mean_pct, fill = Group),
  stat = "identity", width = 0.6, color = "black",
  alpha = 0.6
) +
  geom_errorbar(
    data = df_line_summary,
    aes(x = Group, ymin = mean_pct,ymax = mean_pct + sd_pct),
    width = 0.2
  ) +
  geom_jitter(
    data = df_line,
    aes(x = Group, y = percentage, color = Group),
    width = 0.15, size = 2, alpha = 0.8,
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(
    #title = "Mean LINE Percentage by Group",
    x = "Group",
    y = "DC Proportion in Urine CfRNA"
    #fill = "Sex",
    #color = "Sex"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
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
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.position = "none"  
  )+
  scale_y_continuous(expand = c(0, 0))+
  expand_limits(y = y_max*1.3) + 
  annotate("text", x = 1.5, y = y_max*1.2, label = "***",size=5) +
  annotate("segment", x = 1, xend = 2, y = y_max*1.18, yend = y_max*1.18,linewidth=0.3)

DC_Group_plot


df_line <- data_merge %>% 
  filter(celltype == "DC", !is.na(MIBC))
df_line_summary <- df_line %>%
  group_by(MIBC) %>%
  dplyr::summarise(
    mean_pct = mean(percentage),
    sd_pct = sd(percentage)
  )
custom_colors <- c("MIBC" = '#65AE00', "NMIBC" = "#D3B356")

y_max <- max(df_line$percentage)

DC_MIBC_plot<-ggplot()+geom_bar(
  data = df_line_summary,
  aes(x = MIBC, y = mean_pct, fill = MIBC),
  stat = "identity", width = 0.6, color = "black",
  alpha = 0.6
) +
  geom_errorbar(
    data = df_line_summary,
    aes(x = MIBC, ymin = mean_pct,ymax = mean_pct + sd_pct),
    width = 0.2
  ) +
  geom_jitter(
    data = df_line,
    aes(x = MIBC, y = percentage, color = MIBC),
    width = 0.15, size = 2, alpha = 0.8,
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(
    #title = "Mean LINE Percentage by Group",
    x = "Subgroup",
    y = "DC Proportion in Urine CfRNA"
    #fill = "Sex",
    #color = "Sex"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
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
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.position = "none"  
  )+
  scale_y_continuous(expand = c(0, 0))+
  expand_limits(y = y_max*1.3) 
  annotate("text", x = 1.5, y = y_max*1.2, label = "*",size=5) +
  annotate("segment", x = 1, xend = 2, y = y_max*1.18, yend = y_max*1.18,linewidth=0.3)

DC_MIBC_plot

combined_plot <- Sex_Prostate_plot+Immuno_Group_plot + DC_Group_plot + DC_MIBC_plot
combined_plot

ggsave("Compare_group_celltype.png",height=10,width=10,units = 'in',dpi=600)
ggsave("Compare_group_celltype.pdf",height=10,width=10,units = 'in', device = cairo_pdf)