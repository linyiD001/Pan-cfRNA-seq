library(dplyr) 
library(readxl)
#BiocManager::install("ggnewscale")
library(ggnewscale)
setwd("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/QC/starstandard")

BC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BC")
BS_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BS")
HC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "HC")
clinic <- bind_rows(BC_clinicial_data,BS_clinicial_data,HC_clinicial_data)


clinic$Plasma<-NA
clinic$Urine<-NA
clinic$Tumor<-NA
clinic$NAT<-NA

cfRNA_count<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/cpm_matrix_from_star.csv",header=T,row.names = 1)
QC_list<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/QC/starstandard/QC_starstandard_list_withoutR.csv",row.names = 1)

tissue_count<-read.csv("/dssg/home/acct-dahan/share/BC/tissue/RESULT/MAPPING/star/unanno_count.csv",header=T,row.names = 1)

idx_T<-grep("T",colnames(tissue_count))
T_tissue_count<-tissue_count[,idx_T]
colnames(T_tissue_count)<-sub(paste0("*", "T"), "", colnames(T_tissue_count))

idx_P<-grep("P",colnames(tissue_count))
P_tissue_count<-tissue_count[,idx_P]
colnames(P_tissue_count)<-sub(paste0("*", "P"), "", colnames(P_tissue_count))

idx_N<-grep("N",colnames(cfRNA_count))
N_tissue_count<-cfRNA_count[,idx_N]
colnames(N_tissue_count)<-sub(paste0("*", "_N"), "", colnames(N_tissue_count))

idx_X<-grep("X",colnames(cfRNA_count))
X_tissue_count<-cfRNA_count[,idx_X]
colnames(X_tissue_count)<-sub(paste0("*", "_X"), "", colnames(X_tissue_count))

clinic[which(clinic$Patient %in% colnames(X_tissue_count)),15]<-'plasma'
clinic[which(clinic$Patient %in% colnames(N_tissue_count)),16]<-'urine'
clinic[which(clinic$Patient %in% colnames(T_tissue_count)),17]<-'Tumor'
clinic[which(clinic$Patient %in% colnames(P_tissue_count)),18]<-'NAT'
clinic<-as.data.frame(clinic)
idx_H<-grep("H",clinic$Patient)
clinic$group[idx_H]<-"HC"
idx_S<-grep("S",clinic$Patient)
clinic$group[idx_S]<-"BS"
idx_C<-grep("C",clinic$Patient)
clinic$group[idx_C]<-"BC"

idx_R<-grep("C89",clinic$Patient)
clinic$group[idx_R]<-"R"
idx_R<-grep("C103",clinic$Patient)
clinic$group[idx_R]<-"R"
idx_R<-grep("C107",clinic$Patient)
clinic$group[idx_R]<-"R"
clinic<-clinic[-which(clinic$group=="R"),]
write.csv(clinic,"filtered_sample_summary.csv")

library(ggplot2)
library(dplyr)
library(ggforce)

summary_data <- clinic%>%
  dplyr::group_by(group,SEX,MIBC,Cis,Grade,Plasma,Urine,Tumor,NAT) %>%
  dplyr::summarise(Count = n())

summary_data <- clinic %>%
  dplyr::mutate(group = factor(group, levels = c("C", "S","H","R")),
                Sex = factor(SEX, levels = c("Male","Female")),
                MIBC =factor(MIBC,levels=c ("MIBC","NMIBC","0","NA")) ,
                Plasma = factor(Plasma),
                Urine=factor(Urine),
                Tumor=factor(Tumor),
                NAT = factor(NAT))

aggregated_data <- clinic %>%
  dplyr::group_by(group, SEX,MIBC, Plasma, Urine,Tumor,NAT) %>%
  dplyr::summarise(Count = n()) %>%
  dplyr::ungroup()

library(gridExtra)



main<-ggplot(aggregated_data, aes(y = Count, group = group)) +
  geom_col(aes(fill = "White Center", x = 2.75), width = 2, fill = "white") +  # 添加白色背景
  geom_col(data=aggregated_data,aes(fill = group, x = 2.5), width = .4, color = "gray") + 
  scale_fill_manual(aesthetics = "fill", values = (c("BC" = '#DB498E', "BS" = "#FCAE59","HC" = "#379FB4")),
                    name = "Patient Group",guide = guide_legend(title.position = "top", order = 1),
                    breaks = c("BC", "BS","HC"))+
  new_scale_fill() +
  geom_col(aes(fill = SEX, x = 2.9), width = .4, color = "gray") + 
  scale_fill_manual(aesthetics = "fill", values = c('Male' = "#406A93", "Female" = "#EFA9AE"),
                    na.value = "white",name = "SEX",guide = guide_legend(title.position = "top", order = 2),
                    breaks = c("Male", "Female"))+
  new_scale_fill() +
  geom_col(aes(fill = MIBC, x = 3.2), width = .4, color = "gray") +
  scale_fill_manual(aesthetics = "fill", values = c("MIBC" = '#65AE00', "NMIBC" = "#D3B356"),
                    na.value = "white",name = "Pathology Type",guide = guide_legend(title.position = "top", order = 3),
                    breaks = c("MIBC", "NMIBC"))+
  new_scale_fill() +
  geom_col(aes(fill = Plasma, x = 3.4), width = .4, color = "gray") +
  geom_col(aes(fill = Urine, x = 3.6), width = .4, color = "gray") +
  geom_col(aes(fill = Tumor, x = 3.7), width = .2, color = "gray") +
  geom_col(aes(fill = NAT, x = 3.8), width = .1, color = "gray") +
  scale_fill_manual(aesthetics = "fill", values = c('plasma' = "#EEE0CC", 'urine' = "#B6AFAD", 
                                                    'Tumor' = "azure3", 'NAT' = "azure4"),
                    na.value = "white",name = "Sample Type",guide = guide_legend(title.position = "top", order = 4), breaks = c("plasma", "urine", "Tumor", "NAT"))+
  coord_polar(theta = 'y') +
  theme_void() + 
  theme(panel.border = element_blank(),
        legend.position = "right",
        legend.box = "vertical",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.box.spacing = unit(0.5, "cm"),
        legend.text.align = 0)

main
ggsave("filtered_sample_info.png",height=6,width=6,dpi=1200, units = "in")
ggsave("filtered_sample_info.pdf",height=6,width=6,dpi=1200, units = "in",device = cairo_pdf)

length(which(clinic$group=="BC"))
length(which(clinic$group=="BS"))
length(which(clinic$group=="HC"))


length(which(clinic$Plasma=='plasma' & clinic$group=="BC"))
length(which(clinic$Plasma=='plasma' & clinic$group=="BS"))
length(which(clinic$Plasma=='plasma' & clinic$group=="HC"))

length(which(clinic$Urine=='urine' & clinic$group=="BC"))
length(which(clinic$Urine=='urine' & clinic$group=="BS"))
length(which(clinic$Urine=='urine' & clinic$group=="HC"))

length(which(clinic$Tumor=='Tumor'))
length(which(clinic$NAT=='NAT'))

length(which(clinic$Urine=='urine' & clinic$Plasma=='plasma' 
             & clinic$Tumor=='Tumor' & clinic$NAT=='NAT'))

length(which(clinic$Urine=='urine'& clinic$Plasma=='plasma' & clinic$group=="BS"))
length(which(clinic$Urine=='urine'& clinic$Plasma=='plasma' & clinic$group=="HC"))
length(which(clinic$Urine=='urine'& clinic$Plasma=='plasma' & clinic$group=="BC"))
