
gene_fit_df<-read.csv("Paired_Prediction.csv",header=T,row.names = 1)
#gene_fit_df_select <- gene_fit_df[which(gene_fit_df$P_overall<0.05),]

#load("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/cell_type_marker/paired_cor/lasso_results_top100_pred.RData")
load("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/cell_type_marker/paired_cor/lasso_results_top10_pred_selected.RData")
lasso_summary <- data.frame(
  Gene = names(lasso_results),
  R2 = sapply(lasso_results, function(x) x$R2),
  n_selected = sapply(lasso_results, function(x) x$n_selected)
)

lasso_results$IGF2


lasso_summary <- lasso_summary[order(-lasso_summary$R2), ]
# 给两个数据框添加来源标签
gene_fit_df$Source <- "Single gene"
lasso_summary$Source <- "Multiple genes"

# 只保留 R2 列和来源
plot_df <- bind_rows(
  gene_fit_df %>% select(R2, Source),
  lasso_summary %>% select(R2, Source)
)

# 绘图
ggplot(plot_df, aes(x = R2, fill = Source)) +
  geom_density(alpha = 0.5) +
  geom_vline(data = plot_df %>% group_by(Source) %>%
               summarise(mean_R2 = mean(R2, na.rm = TRUE)),
             aes(xintercept = mean_R2, color = Source),
             linetype = "dashed", size = 1) +
  labs(
    title = "Comparison of R² Distributions",
    x = expression(R^2),
    y = "Density",
    fill = "Source",
    color = "Mean R²"
  ) +
  scale_fill_manual(values = c("Single gene" = "#377EB8", 
                               "Multiple genes" = "#69b3a2")) +
  scale_color_manual(values = c("Single gene" = "#377EB8", 
                                "Multiple genes" = "#69b3a2")) +
  theme_minimal(base_size = 14)

ggsave("Tumor_specific_lasso_multi_prediction.png",height=6,width=12,dpi=600,unit='in')
ggsave("Tumor_specific_lasso_multi_prediction.pdf",height=6,width=12)

#ggsave("NAT_specific_lasso_multi_prediction.png",height=6,width=12,dpi=600,unit='in')
#ggsave("NAT_specific_lasso_multi_prediction.pdf",height=6,width=12)


gene <- "IGF2"

library(ggplot2)
# 提取并合并为数据框
df <- data.frame(
  selected_features = lasso_results[[gene]][5],
  coef = lasso_results[[gene]][6]
)

# 按系数大小排序
df$selected_features <- factor(df$selected_features, levels = df$selected_features[order(df$coef, decreasing = TRUE)])
R2_val = as.numeric(lasso_results[[gene]][3])
# 绘制柱状图
# 绘图
ggplot(df, aes(x = selected_features, y = coef, fill = coef)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_gradient(low = "skyblue", high = "red")  +
  geom_text(aes(label = round(coef, 2)), vjust = -0.5, size = 3.5) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("R² = ", round(R2_val, 2)),
    hjust = 1.2, vjust = 1.5,
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste("LASSO Coefficients for", gene, "Model"),
    x = "Features",
    y = "Coefficient"
  )+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(gene,"_coefficient.png"),height=6,width=12,dpi=600,unit='in')
ggsave(paste0(gene,"_coefficient.pdf"),height=6,width=12)
