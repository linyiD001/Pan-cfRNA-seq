library(tidymodels)
library(ggplot2)
library(dplyr)
library(caret)
library("glmnet")
rm(list=ls())
set.seed(250818)

setwd("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/model/star_standard/N_SH")
train_data<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/model/star_standard/train_data_N.csv",header=T,row.names = 1)
test_data<-read.csv("/dssg/home/acct-dahan/share/BC/cfRNA/RESULT/ANALYSIS/model/star_standard/test_data_N.csv",header=T,row.names = 1)

idx_C<-which(train_data$outcome=='3')
train_data<-train_data[-idx_C,]

idx_C<-which(test_data$outcome=='3')
test_data<-test_data[-idx_C,]

#idx_H<-grep('H',colnames(count_analysis))
#coldata$type[idx_H]<-'1'
#idx_C<-grep('C',colnames(count_analysis))
#coldata$type[idx_C]<-'3'

# Identify predictors and response
train_x <- as.matrix(train_data[,1:(ncol(train_data)-1)])
train_y <- train_data[,ncol(train_data)]
train_y<-as.numeric(train_y)
test_x <- as.matrix(test_data[,1:(ncol(test_data)-1)])
test_y <- test_data[,ncol(test_data)]
test_y<-as.numeric(test_y)

lambdas = NULL
for (i in 1:100)
{
  fit <- cv.glmnet(train_x, train_y, family = "binomial",type.measure = "class",nfolds=10)
  errors = data.frame(fit$lambda,fit$cvm)
  lambdas <- rbind(lambdas,errors)
}
# take mean cvm for each lambda
lambdas <- aggregate(lambdas[, 2], list(lambdas$fit.lambda), mean)

# select the best one
bestindex = which(lambdas[2]==min(lambdas[2]))
bestlambda = lambdas[bestindex,1]
fit<- glmnet(train_x, train_y, lambda=bestlambda, family = "binomial",type.measure = "class")


# Save the model
save.image("N_SH.RData")
saveRDS(fit, file = "glmnet_model.rds")

#cvfit<-glmnet_model

# Prepare training data
train_x <- as.matrix(train_data[, 1:(ncol(train_data) - 1)])
train_y <- as.factor(train_data$outcome)

# Predict probabilities for training set using the best lambda
train_pred <- predict(fit, newx = train_x,s = bestlambda, type = "response")

# Prepare testing data
test_x <- as.matrix(test_data[, 1:(ncol(test_data) - 1)])
test_y <- as.factor(test_data$outcome)

# Predict probabilities for testing set using the best lambda
test_pred <- predict(fit, newx = test_x,s = bestlambda,type = "response")


train <- roc(train_y, train_pred,ci=TRUE, boot.n=2000, ci.alpha=0.95, stratified=T)
test <- roc(test_y, test_pred,ci=TRUE, boot.n=2000, ci.alpha=0.95, stratified=T)

#ci_train <- ci.se(train, specificities=seq(0, 1, l=50))
#ci_test <- ci.se(test, specificities=seq(0, 1, l=50))

ci_train_auc <- ci.auc(train, boot.n=2000, ci.alpha=0.95, stratified=T)
ci_test_auc <- ci.auc(test, boot.n=2000, ci.alpha=0.95, stratified=T)

train_auc<-auc(train)
legend_train_text <- paste0("Train AUC = ", round(train_auc, 2), " (95% CI = ", round(ci_train_auc[1], 2), " - ", round(ci_train_auc[3], 2), ")")
test_auc<-auc(test)
legend_test_text <- paste0("Test AUC = ", round(test_auc, 2), " (95% CI = ", round(ci_test_auc[1], 2), " - ", round(ci_test_auc[3], 2), ")")

# Multiple curves:
g2 <- ggroc(list(train = train, test = test), legacy.axes = TRUE) +
  ylab("Sensitivity") + xlab("1-Specificity") +
  geom_segment(aes(x = 1, xend = 0, y = 1, yend = 0), color = "grey", linetype = "dashed") +
  annotate("text", x = 0.7, y = 0.35, label = legend_train_text, size = 3.7) +
  annotate("text", x = 0.7, y = 0.25, label = legend_test_text, size = 3.7) +
  scale_color_manual(values = c("train" = "#91191E", "test" = "#17734F")) +  # ← 修改颜色
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),
    legend.spacing.y = unit(1, 'cm'),
    legend.key.spacing.y = unit(0.1, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
g2
ggsave("combined_ROC_curve_train_test_with_AUC_ggplot_95_CI.jpg", g2, width = 5, height = 4, units = "in", dpi = 300)



# Prepare the data frame for training scores
train_scores_df <- data.frame(
  Score = as.numeric(train_pred),
  Group = as.factor(train_y),
  Dataset = "Training"
)

# Prepare the data frame for test scores
test_scores_df <- data.frame(
  Score = as.numeric(test_pred),
  Group = as.factor(test_y),
  Dataset = "Test"
)

# Combine both training and test scores into one data frame
scores_df <- bind_rows(train_scores_df, test_scores_df)

# Plot the boxplot
boxplot_plot <- ggplot(scores_df, aes(x = Group, y = Score, fill = Dataset)) +
  geom_boxplot(outlier.size = 1, outlier.shape = 21, alpha = 0.7) +
  labs(title = "Boxplot of Prediction Scores by Group",
       x = "Group",
       y = "Prediction Score",
       fill = "Dataset") +
  scale_fill_manual(values = c("blue", "red")) +
  theme_minimal() +
  theme(
    text = element_text(size = 8),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )
boxplot_plot
ggsave("boxplot_plot.jpg", boxplot_plot, width = 5, height = 4, units = "in", dpi = 300)

g2 <- ggroc(list(train = train, test = test), legacy.axes = TRUE) +
  ylab("Sensitivity") + xlab("1-Specificity") +
  geom_segment(aes(x = 1, xend = 0, y = 1, yend = 0), color = "grey", linetype = "dashed") +
  annotate("text", x = 0.7, y = 0.35, label = legend_train_text, size = 3.5) +
  annotate("text", x = 0.7, y = 0.25, label = legend_test_text, size = 3.5) +
  scale_color_manual(values = c("train" = "#91191E", "test" = "#17734F")) +  # ← 修改颜色
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),
    legend.spacing.y = unit(1, 'cm'),
    legend.key.spacing.y = unit(0.1, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
g2

ggsave("N_S_H_ROC_plot.png", g2, width = 5, height = 4, units = "in", dpi = 300)
ggsave("N_S_H_ROC_plot.pdf", g2, width = 5, height = 4, units = "in", dpi = 300)


