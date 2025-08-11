library(pheatmap)

data<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/cell_type_marker/NX_urinecell_cibersort/cf_cibersort_urinecell/CIBERSORTx_Job2_Results.csv",row.names=1)
data<-data[grep("N",rownames(data)),]
data<-data[-which(rownames(data)=='C107_N'),]

pheatmap(data[,1:14], cluster_rows=T, show_rownames=T,cluster_cols=T,show_colnames = T, 
         #annotation_col=coldata_all,
         #annotation_colors = ann_colors,na_col = "gray",
         scale = "none")

coldata<-matrix(ncol = 2,nrow = 104)
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
  Sex = c("Male" = "#406A93", "Female" = "#EFA9AE"),
  MIBC = c("MIBC" = "#65AE00", "NMIBC" = "#D3B356"),
  Group = c("BC"="#DB498E", "BS"= "#FCAE59", "HC"="#379FB4"),
  Tumor_size=colorRampPalette(c("gray", "#2C1714"))(100),
  Age = colorRampPalette(c("#EEE0CC", "#91191E"))(100),
  StageT=c('4'="#FF3300",'3'='#FF6633','2'='#FF9966','1'='#FFCC66','a'='#FFFF00'),
  Grade=c('Low'='#7A84D0','High'="#D78203"),
  sample=c('T'='#8E5D2A','NAT'="#64BE52"),
  Tumor_number=c('single'='#90D4A4','multiple'="#C49406")
)


data_heatmap<-as.matrix(t(data))
data_heatmap<-data_heatmap[1:14,]
#CairoPDF(file = "cfRNA_urine_cibersort.pdf", width = 12, height = 5,family='Arial')
CairoPNG(file = "cfRNA_urine_cibersort.png", width = 12, height = 5,res=600,units='in',family='Arial')


pheatmap(data_heatmap,
         cluster_rows = TRUE,
         show_rownames = T,
         cluster_cols = F,
         show_colnames = F,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = coldata_all, 
         annotation_colors = ann_colors, 
         fontsize_row =10,
         fontsize = 10,
         na_col = "gray",
         scale = "none",
         #color = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
         color = colorRampPalette(c("blue", "white", "red"))(100),
         legend=TRUE,
         name = 'Z-score')

dev.off()

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

data_stack_S <- data_stack_S %>%
  arrange(desc(Value)) %>% 
  mutate(
    Percentage = Value / sum(Value) * 100,  
    Cell_Type = factor(Cell_Type, levels = Cell_Type)  
  )

gradient_colors <- c(
  "#E8AEB0",  
           "#F4C28B", 
           "#E8DE7D",  
           "#A9E0A9", 
           "#B9E8C3",  
           "#B8E6B8",
           "#A9DCE3", 
           "#A7CCE5", 
           "#C0DDEE", 
           "#B4B4DC",  
           "#E6CFC2",  
           "#D4C89A",  
           "#CFCFCF",  
           "#F5F58C"   
)

p1<-ggplot(data_stack_S, aes(x = "", y = Value, fill = Cell_Type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +  
  coord_polar(theta = "y") +  
  labs(title = "", x = NULL, y = NULL) +
  scale_fill_manual(values = gradient_colors) + 
  theme_void() +  
  geom_text(
    aes(label = ifelse(Percentage >= 2, paste0(round(Percentage, 1), "%"), "")), 
    position = position_stack(vjust = 0.5), size = 3
  )  
dev.off()
ggsave("BS_cell_type.png",height=6,width=6,dpi=600)

data_stack_H <- data_stack_H %>%
  arrange(desc(Value)) %>%  
  mutate(
    Percentage = Value / sum(Value) * 100,  
    Cell_Type = factor(Cell_Type, levels = Cell_Type)  
  )

ggplot(data_stack_H, aes(x = "", y = Value, fill = Cell_Type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +  
  coord_polar(theta = "y") +  
  labs(title = "Cell Type Proportions in Urine cfRNA in HC", x = NULL, y = NULL) +
  scale_fill_manual(values = gradient_colors) +  
  theme_void() + 
  geom_text(
    aes(label = ifelse(Percentage >= 2, paste0(round(Percentage, 1), "%"), "")), 
    position = position_stack(vjust = 0.5), size = 3
  )  
ggsave("HC_cell_type.png",height=6,width=6,dpi=600)

data_stack_C <- data_stack_C %>%
  arrange(desc(Value)) %>%  
  mutate(
    Percentage = Value / sum(Value) * 100,  
    Cell_Type = factor(Cell_Type, levels = Cell_Type)  
  )

gradient_colors <- colorRampPalette(c("lightblue","red", "yellow",'green'))(14)


ggplot(data_stack_C, aes(x = "", y = Value, fill = Cell_Type)) +
  geom_bar(stat = "identity", width = 1, color = "white") + 
  coord_polar(theta = "y") +  
  labs(title = "Cell Type Proportions in Urine cfRNA in BC", x = NULL, y = NULL) +
  scale_fill_manual(values = gradient_colors) +  
  theme_void() +  
  geom_text(
    aes(label = ifelse(Percentage >= 2, paste0(round(Percentage, 1), "%"), "")), 
    position = position_stack(vjust = 0.5), size = 3
  )  
ggsave("BC_cell_type.png",height=6,width=6,dpi=600)


gradient_colors <- c(
           "#E8AEB0", 
           "#F4C28B", 
           "#E8DE7D", 
           "#A9E0A9",
           "#D0B7E1", 
           "#F2B5D4", 
           "#A9DCE3",  
           "#A7CCE5",  
           "#C0DDEE", 
           "#B4B4DC",  
           "#E6CFC2",  
           "#D4C89A",  
           "#CFCFCF",  
           "#F5F58C"   
)

gradient_colors <- rev(gradient_colors) 
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
      Cell_Type = factor(Cell_Type, levels = unique(Cell_Type))  
    )
  
  ggplot(data, aes(x = "", y = Value, fill = Cell_Type)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y", start = 0) +  
    labs(title = title, x = NULL, y = NULL) +
    scale_fill_manual(values = cell_type_colors) + 
    theme_void() +
    
    geom_text(
      aes(label = ifelse(Percentage >= 2, paste0(round(Percentage, 1), "%"), "")),
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
      legend.spacing.x = unit(0.5, 'cm'),   
      legend.spacing.y = unit(0.5, 'cm'),   
      legend.key.height = unit(0.5, "cm"),    
      legend.key.spacing.y = unit(0.1, "cm"),    
      legend.key.width = unit(0.5, "cm"),     
      legend.box.spacing = unit(0.5, 'cm'),
      #panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      #legend.background = element_rect(fill = "white", color = "black", size = 0.5), 
      legend.key = element_rect(color = "black", size = 0.1)  
    )
}
p1 <- make_pie(data_stack_H, "Urine cfRNA cell proportion in healthy controls")+ theme(legend.position = "none")
p2 <- make_pie(data_stack_C, "Urine cfRNA cell proportion in bladder cancer patients")+ theme(legend.position = "none")
p3 <- make_pie(data_stack_S, "Urine cfRNA cell proportion in bladder stone patients")

combined_plot <- p1 + p2 + p3 +plot_layout(ncol = 3)

combined_plot


ggsave("ALL_cell_type.png",height=6,width=20,units = 'in',dpi=600)
ggsave("ALL_cell_type.pdf",height=6,width=20,units = 'in', device = cairo_pdf)
