##paired raw count cor gene level
library(dplyr)
library(VariableScreening)
library(biomaRt)
setwd("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/cell_type_marker/paired_cor")

CPM_cfRNA<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/cpm_matrix_from_star.csv",header=T,row.names = 1)
CPM_tissue<-read.csv("/dssg/home/acct-dahan/share/BC/tissue/RESULT/NORMALIZATION/cpm_matrix_from_star.csv",header=T,row.names = 1)

idx_T<-grep("T",colnames(CPM_tissue))
T_CPM_tissue<-CPM_tissue[,idx_T]
rownames(T_CPM_tissue)<-rownames(CPM_tissue)

idx_P<-grep("P",colnames(CPM_tissue))
P_CPM_tissue<-CPM_tissue[,idx_P]
rownames(P_CPM_tissue)<-rownames(CPM_tissue)

idx_N<-grep("N",colnames(CPM_cfRNA))
N_CPM_cfRNA<-CPM_cfRNA[,idx_N]
rownames(N_CPM_cfRNA)<-rownames(CPM_cfRNA)

colnames(T_CPM_tissue) <- sub(paste0("*", "T"), "", colnames(T_CPM_tissue))
colnames(P_CPM_tissue) <- sub(paste0("*", "P"), "", colnames(P_CPM_tissue))
colnames(N_CPM_cfRNA) <- sub(paste0("*", "_N"), "", colnames(N_CPM_cfRNA))

all_paired_list<-Reduce(intersect, list(colnames(T_CPM_tissue),
                                        colnames(P_CPM_tissue),
                                        colnames(N_CPM_cfRNA)))

all_paired_gene<-Reduce(intersect, list(rownames(T_CPM_tissue),
                                        rownames(P_CPM_tissue),
                                        rownames(N_CPM_cfRNA)))

#all_paired_list<-as.data.frame(all_paired_list)
#write.table(all_paired_list,"all_3_paired_list.txt",col.names = F,row.names=F,quote = F)

#expr_10 <- rowMeans(N_CPM_cfRNA > 10)
#keep_genes <- (expr_10 >= 0.5)
#keep_genes_names <- names(keep_genes[keep_genes == TRUE])

N_C_H_DEG<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/DEG_standardstar/N_type_C_vs_H_result.csv",header=T,row.names = 1)
padjthreshold<-0.05
log2FoldChangethreshold<-4
basemeanthreshold<-5
up_N_C_H_DEG<-rownames(N_C_H_DEG[which(N_C_H_DEG$padj<padjthreshold&N_C_H_DEG$log2FoldChange>log2FoldChangethreshold&N_C_H_DEG$baseMean>basemeanthreshold),])

#keep_genes_names <- c("IGF2","H19","IVL")
keep_genes_names <- up_N_C_H_DEG
keep_genes_names <- keep_genes_names[-grep("ENSG",keep_genes_names)]
keep_genes_names
#expr_10_T <- rowMeans(T_CPM_tissue > 10)
#keep_genes_T <- (expr_10_T >= 0.5)
expr_10_P <- rowMeans(P_CPM_tissue > 10)
#keep_genes_P <- (expr_10_P >= 0.5)

keep_genes_tissue <- (expr_10_T >= 0.5) | (expr_10_P >= 0.5)

keep_genes_names_tissue <- names(keep_genes_tissue[keep_genes_tissue == TRUE])
###24 all paired samples

#T_P <- read.csv("/dssg/home/acct-dahan/share/BC/tissue/RESULT/NORMALIZATION/type_T_vs_P.csv",header=T,row.names=1)
#Tumor_specific <- T_P[which(T_P$log2FoldChange>1&T_P$padj<0.05),]
#NAT_specific <- T_P[which(T_P$log2FoldChange<(-1)&T_P$padj<0.05),]

library(VariableScreening)
# 2️⃣ 初始化存储列表
all_gene_screen <- list()

#all_paired_gene<-"H19"
# 3️⃣ 遍历每个尿液基因
for (g in keep_genes_names) {
  paired_T_CPM_tissue<- T_CPM_tissue[which(rownames(T_CPM_tissue) %in% c(keep_genes_names_tissue,g)),which(colnames(T_CPM_tissue) %in% all_paired_list)]
  paired_P_CPM_tissue<- P_CPM_tissue[which(rownames(P_CPM_tissue) %in% c(keep_genes_names_tissue,g)),which(colnames(P_CPM_tissue) %in% all_paired_list)]
  paired_N_CPM_cfRNA <- N_CPM_cfRNA[which(rownames(N_CPM_cfRNA) %in% keep_genes_names),which(colnames(N_CPM_cfRNA) %in% all_paired_list)]
  
  log2_paired_T_CPM_tissue <- log2(paired_T_CPM_tissue + 1)
  log2_paired_P_CPM_tissue <- log2(paired_P_CPM_tissue + 1)
  log2_paired_N_CPM_cfRNA <- log2(paired_N_CPM_cfRNA + 1)
  
  rownames(log2_paired_T_CPM_tissue) <- paste0(rownames(log2_paired_T_CPM_tissue),"_T")
  rownames(log2_paired_P_CPM_tissue) <- paste0(rownames(log2_paired_P_CPM_tissue),"_P")
  
  X_tissue <- t(rbind(log2_paired_T_CPM_tissue, log2_paired_P_CPM_tissue))  # n x (2*genes)
  X_tissue <- as.data.frame(X_tissue)
  X_tissue <- X_tissue[, colSums(is.na(X_tissue)) == 0]  # 去掉含NA的列
  y <- t(log2_paired_N_CPM_cfRNA[g, ])
  
  # 跳过全 NA 或方差为0的基因
  if (all(is.na(y)) || sd(y, na.rm = TRUE) == 0) next
  
  # DC-SIS 筛选
  screen_data <- screenIID(X_tissue, y, method = "DC-SIS")
  screen_data <- as.data.frame(screen_data)
  screen_data$Gene <- g  # 标记对应的尿液基因
  screen_data$Feature <- rownames(screen_data)
  
  # 重排列顺序
  screen_data <- screen_data %>% dplyr::select(Gene, Feature, everything())
  
  all_gene_screen[[g]] <- screen_data
}
save(all_gene_screen, file = "all_gene_screen_df_select_all.RData")

# 4️⃣ 合并所有结果
all_gene_screen_df <- bind_rows(all_gene_screen)

#write.csv("all_gene_screen_df.csv")

# 初始化存储结果
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=110)
#ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",mirror='www')

geneid=rownames(T_CPM_tissue)
mart<-useMart("ensembl")
data=listDatasets(mart)
listFilters(ensembl)
symbols <- getBM(attributes=c('hgnc_symbol','gene_biotype'),
                 filters = 'hgnc_symbol', values = geneid,
                 mart = ensembl)
symbols
dedu_symbols<-symbols[!duplicated(symbols$hgnc_symbol),]
dedu_symbols
protein_symbols <- dedu_symbols[which(dedu_symbols$gene_biotype=='protein_coding'),]
library(stringr)

lasso_results <- list()
#g='IGF2'
# 遍历每个尿液基因进行 LASSO 拟合
for (g in keep_genes_names) {
  paired_T_CPM_tissue<- T_CPM_tissue[which(rownames(T_CPM_tissue) %in% c(keep_genes_names_tissue,g)),which(colnames(T_CPM_tissue) %in% all_paired_list)]
  paired_P_CPM_tissue<- P_CPM_tissue[which(rownames(P_CPM_tissue) %in% c(keep_genes_names_tissue,g)),which(colnames(P_CPM_tissue) %in% all_paired_list)]
  paired_N_CPM_cfRNA<- N_CPM_cfRNA[which(rownames(N_CPM_cfRNA) %in% keep_genes_names),which(colnames(N_CPM_cfRNA) %in% all_paired_list)]
  
  log2_paired_T_CPM_tissue <- log2(paired_T_CPM_tissue + 1)
  log2_paired_P_CPM_tissue <- log2(paired_P_CPM_tissue + 1)
  log2_paired_N_CPM_cfRNA <- log2(paired_N_CPM_cfRNA + 1)
  rownames(log2_paired_T_CPM_tissue) <- paste0(rownames(log2_paired_T_CPM_tissue),"_T")
  rownames(log2_paired_P_CPM_tissue) <- paste0(rownames(log2_paired_P_CPM_tissue),"_P")
  X_tissue <- t(rbind(log2_paired_T_CPM_tissue, log2_paired_P_CPM_tissue))  # n x (2*genes)
  X_tissue <- as.data.frame(X_tissue)
  X_tissue <- X_tissue[, colSums(is.na(X_tissue)) == 0] 
  
  rm(screen_df)
  screen_df <- all_gene_screen[[g]]
  gene_prefix <- str_extract(screen_df$Feature, "^[^_]+")
  screen_df <- screen_df[gene_prefix %in% protein_symbols$hgnc_symbol, ]  
  top_features <- screen_df %>%
    arrange(desc(-rank)) %>%
    slice_head(n = 10) %>%
    pull(Feature)
  top_features <- unique(c(top_features,paste0(g,"_T"),paste0(g,"_P")))
  top_features
  X <- X_tissue[, top_features]
  y <- as.numeric(log2_paired_N_CPM_cfRNA[g, ])
  X <- data.matrix(X)
  cv_fit <- cv.glmnet(X, y, alpha = 1, nfolds = 5,family = "gaussian")
  best_lambda <- cv_fit$lambda.min
  fit <- glmnet(X, y, alpha = 1, lambda = best_lambda)
  best_coef <- coef(fit, s = best_lambda)
  R2_val<- R2_val<- fit$dev.ratio
  coef_df <- as.data.frame(as.matrix(best_coef))
  coef_df <- coef_df[coef_df[,1] != 0, , drop = FALSE]
  coef_df <- tibble::rownames_to_column(coef_df, "Feature")
  colnames(coef_df)[2] <- "Coefficient"
  selected_features <- coef_df$Feature[coef_df$Feature != "(Intercept)"]
  selected_coef <- coef_df$Coefficient[coef_df$Feature != "(Intercept)"]
  n_selected <- length(selected_features)
  list(fit=cv_fit,selected_features=selected_features, n_selected=n_selected)
  # 保存结果
  lasso_results[[g]] <- list(
    model = fit,
    R2 = R2_val,
    n_selected = n_selected,
    selected_features = selected_features,
    coef= selected_coef
  )
}

## ===== 汇总并保存结果 =====
#save(lasso_results,  file = "lasso_results_top10_pred_selected.RData")
#save.image(file = "lasso_results_top10_pred_selected_Rimage.RData")

#save(lasso_results,  file = "lasso_results_top100_pred_NAT.RData")

lasso_summary <- data.frame(
  Gene = names(lasso_results),
  R2 = sapply(lasso_results, function(x) x$R2),
  n_selected = sapply(lasso_results, function(x) x$n_selected)
)

lasso_results$IGF2

lasso_summary <- lasso_summary[order(-lasso_summary$R2), ]
head(lasso_summary)

ggplot(lasso_summary, aes(x = R2)) +
  geom_density(fill = "#69b3a2", alpha = 0.5) +   # 填充颜色和透明度
  geom_vline(aes(xintercept = mean(R2, na.rm = TRUE)), 
             color = "red", linetype = "dashed", size = 1) + # 平均值参考线
  labs(title = "Density of R² from LASSO results",
       x = expression(R^2),
       y = "Density") +
  theme_minimal(base_size = 14)

#ggsave("lasso_multi_prediction.png",height=6,width=12,dpi=600,unit='in')
#ggsave("lasso_multi_prediction.pdf",height=6,width=12)

library(circlize)
library(tidyverse)
library(ComplexHeatmap)
library(Cairo)

# --------------------------
# 构建 edge list
# --------------------------

df_list <- map_df(names(lasso_results), function(gene) {
  feats <- lasso_results[[gene]]$selected_features
  coefs <- lasso_results[[gene]]$coef
  
  # ⚠️ 跳过没有选中特征的基因
  if (is.null(feats) || length(feats) == 0) return(NULL)
  
  # ⚠️ 如果 coef 数量与 features 不一致，尝试修正
  if (length(coefs) != length(feats)) {
    if (length(coefs) == 1) {
      coefs <- rep(coefs, length(feats))
    } else {
      coefs <- coefs[seq_len(min(length(coefs), length(feats)))]
      feats <- feats[seq_len(min(length(coefs), length(feats)))]
    }
  }
  
  data.frame(
    Target = gene,
    Feature = feats,
    Coef = coefs
  )
})

df_list_new<-df_list
colnames(df_list_new)[1]<-"Gene"
data_merge <- merge(df_list_new,lasso_summary,by='Gene',all=T)



namelist_tumor <- data_merge[which(paste0(data_merge$Gene,"_T") == data_merge$Feature),1]
namelist_NAT <- data_merge[which(paste0(data_merge$Gene,"_P") == data_merge$Feature),1]

name_both <- intersect(namelist_tumor,namelist_NAT)
name_tumor <- setdiff(namelist_tumor,namelist_NAT)
name_NAT <-  setdiff(namelist_NAT,namelist_tumor)

list_both <- lasso_summary[which(lasso_summary$Gene%in% name_both),]
list_tumor <- lasso_summary[which(lasso_summary$Gene%in% name_tumor),]
list_NAT <- lasso_summary[which(lasso_summary$Gene%in% name_NAT),]
  
colnames(dedu_symbols)[1] <- "Gene"
list_both_gene <- merge(list_both,dedu_symbols,by='Gene',all.x=T)
list_tumor_gene <- merge(list_tumor,dedu_symbols,by='Gene',all.x=T)
list_NAT_gene <- merge(list_NAT,dedu_symbols,by='Gene',all.x=T)

df_list<- df_list %>%
  filter(Target %in% c("NKAIN3","CCER1","TRIM72"))

# 给 Feature 添加类别
df_list <- df_list %>%
  mutate(Source = ifelse(grepl("_T$", Feature), "Tumor", "NAT"))

# --------------------------
# 设置扇区顺序
# --------------------------
targets <- unique(df_list$Target)
features_T <- df_list %>% filter(Source == "Tumor") %>% pull(Feature) %>% unique()
features_P <- df_list %>% filter(Source == "NAT") %>% pull(Feature) %>% unique()

features <- c(features_T, features_P)
sector_order <- c(targets, features)

# --------------------------
# 弦线颜色
# --------------------------
col_fun <- colorRamp2(
  c(min(df_list$Coef), 0, max(df_list$Coef)),
  c("#4B9CD3", "lightyellow", "#D73027")
)
link_colors <- col_fun(df_list$Coef)





# --------------------------
# 自定义 Target 颜色
# --------------------------
target_colors <- c(
  "NKAIN3" = "#FFCC66",
  "CCER1"  = "#66CC99",
  "TRIM72"  = "#CC99FF"
)

# --------------------------
# 绘制 chord 图
# --------------------------
pdf("chord_diagram_targets.pdf", width = 6, height = 5)  
circos.clear()
chordDiagram(
  x = df_list[, c("Target", "Feature")],
  order = sector_order,
  col = link_colors,
  transparency = 0.3,
  annotationTrack = c("grid")
)

# --------------------------
# 第一层：Target 背景颜色注释
# --------------------------
circos.track(
  track.index = 1,
  bg.border = NA,
  ylim = c(0, 1),
  track.height = 0.05,
  panel.fun = function(x, y) {
    sector <- get.cell.meta.data("sector.index")
    if (sector %in% names(target_colors)) {
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                  col = target_colors[sector], border = NA)
    }
  }
)
# --------------------------
# 第二层：Feature 分类注释
# --------------------------
features_T_all <- features[grepl("_T$", features)]
features_P_all <- features[grepl("_P$", features)]

targets_all <- unique(df_list$Target)

circos.track(
  track.index = 1,
  bg.border = NA,
  ylim = c(0, 1),
  track.height = 0.05,
  panel.fun = function(x, y) {
    sector <- get.cell.meta.data("sector.index")
    if (sector %in% features_T_all) {
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                  col = "#FFDDCC", border = NA)
    } else if (sector %in% features_P_all) {
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                  col = "#CCE5FF", border = NA)
    }
  }
)

circos.track(
  track.index = 2,
  bg.border = NA,
  ylim = c(0, 1),
  track.height = 0.05,
  panel.fun = function(x, y) {
    circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                col = "gray", border = NA)
    sector <- get.cell.meta.data("sector.index")
    if (sector %in% c("NKAIN3_T","NKAIN3_P","NKAIN3")) {
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                  col = "#FFCC66", border = NA)}
    if (sector %in% c("CCER1_T","CCER1_P","CCER1")) {
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                  col = "#66CC99", border = NA)}
    if (sector %in% c("TRIM72_T","TRIM72_P","TRIM72")) {
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                  col = "#CC99FF", border = NA)}
  }
  
)

# --------------------------
# 第三层：标签（去掉 _T/_P 后缀）
# --------------------------
circos.trackPlotRegion(
  track.index = 1,
  bg.border = NA,
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    sector <- get.cell.meta.data("sector.index")
    label <- gsub("_[TP]$", "", sector)
    circos.text(
      x = mean(get.cell.meta.data("xlim")),
      y = 0.5,
      labels = label,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.5
    )
  }
)

# --------------------------
# 图例
# --------------------------
lgd_coef <- Legend(
  col_fun = col_fun,
  title = "Coefficient",
  at = c(round(min(df_list$Coef), 2), 0, round(max(df_list$Coef), 2)),
  labels = c("Negative", "0", "Positive"),
  title_gp = gpar(fontsize = 10),
  labels_gp = gpar(fontsize = 9)
)

draw(lgd_coef, x = unit(1, "npc") - unit(10, "mm"), y = unit(0.5, "npc"))
dev.off()

