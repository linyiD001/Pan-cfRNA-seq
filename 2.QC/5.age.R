
BC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BC")
BS_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BS")
HC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "HC")

clinic <- bind_rows(BC_clinicial_data,BS_clinicial_data,HC_clinicial_data)
clinic$group<-NA
clinic[grep('C',clinic$Patient),15]<-"BC"
clinic[grep('S',clinic$Patient),15]<-"BS"
clinic[grep('H',clinic$Patient),15]<-"HC"


idx_out<-grep("C89",clinic$Patient)
clinic$group[idx_out]<-"R"
idx_out<-grep("C107",clinic$Patient)
clinic$group[idx_out]<-"R"
idx_out<-grep("C103",clinic$Patient)
clinic$group[idx_out]<-"R"
clinic<-clinic[-(which(clinic$group=="R")),]


filtered_namelist<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/QC/QC_repeataware_list.csv",header=T,row.names = 1)

N_data_filter<-filtered_namelist[grep("N",filtered_namelist$x),]
  
filtered_namelist_1<-t(as.data.frame(str_split(N_data_filter,"_")))[,1]

clinic<-clinic[which(clinic$Patient %in% filtered_namelist_1),]
ggplot(clinic, aes(x = group, y =  AGE, fill = group)) +
  geom_boxplot() +
  theme_light(base_size=10)+
  theme_classic()

custom_colors <- c("BC" = "#DB498E", "BS" = "#FCAE59", "HC" = "#379FB4")

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
    axis.line = element_line(color = "black"),  
    axis.text.x = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(0.5, 'cm'),   
    legend.key.height = unit(0.5, "cm"),    
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.width = unit(0.5, "cm"),     
    legend.box.spacing = unit(0.5, 'cm'),
    legend.key = element_rect(color = "black", size = 0.1)  
  )

p

ggsave("N_age.png",height=4,width=4,dpi=1200, units = "in")
ggsave("N_age.pdf",height=4,width=4,dpi=1200, units = "in",device = cairo_pdf)



ggplot(clinic, aes(x = group, fill = SEX)) +
  geom_bar(position = "dodge") +  
  theme_classic(base_size = 10) +
  labs(y = "Count", title = "Sample Count by Group and SEX") 

kruskal_result <- kruskal.test(AGE ~ group, data = clinic)
kruskal_result


table_data <- table(clinic$group, clinic$SEX)
chi_result <- chisq.test(table_data)
print(chi_result)

ggsave("N_gender.png",height=6,width=6,dpi=1200)