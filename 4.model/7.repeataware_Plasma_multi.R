library(tidymodels)
library(ggplot2)
library(h2o)
library(h2plots)
library(dplyr)
library(caret)
library("plyr")
library("dplyr")
library("base64enc")
library("labeling")
library("lattice")
library("Matrix")
library("caret")
library("class")
library("e1071")
library("evaluate")
library("ggplot2")
library("gbm")
library("glmnet")
library("gplots")
library("graphics")
library("grDevices")
library("grid")
library("kernlab")
library("KernSmooth")
library("klaR")
library("knitr")
library("labeling")
library("lmtest")
library("markdown")
library("ModelMetrics")
library("modeltools")
library("readr")
library("stats")
library("MASS")
library("corrplot")
library("Hmisc")
library("Boruta")

set.seed(250818)
rm(list=ls())
setwd("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/model/salmon_repeataware_stranded_allgenome/multi_X")

train_data<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/model/salmon_repeataware_stranded_allgenome/train_data_X.csv",header=T,row.names = 1)
test_data<-read.csv("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/model/salmon_repeataware_stranded_allgenome/test_data_X.csv",header=T,row.names = 1)

#idx_S<-grep('S',colnames(count_analysis))
#coldata$type[idx_S]<-'2'
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
fit <- glmnet(train_x, train_y, family = "multinomial", type.multinomial = "grouped")
plot(fit, xvar = "lambda", label = TRUE, type.coef = "2norm")
library(caret)
flds <- createFolds(train_x, k = 10, list = TRUE, returnTrain = FALSE)

cvfit <- cv.glmnet(train_x, train_y, family = "multinomial", type.multinomial = "grouped",nfolds=10)

plot(cvfit)

cvfit$lambda.min
cvfit$lambda.1se

# Extract the best lambda
best_lambda <- cvfit$lambda.min

# Save the model
saveRDS(cvfit, file = "glmnet_model.rds")

# Predict the classes for training and test datasets
train_predictions <- predict(cvfit, newx = train_x, s = best_lambda, type = "class")
test_predictions <- predict(cvfit, newx = test_x, s = best_lambda, type = "class")

# Convert predictions and true labels to factors for multi-class confusion matrix
train_predictions_factor <- as.factor(train_predictions)
test_predictions_factor <- as.factor(test_predictions)
train_y_factor <- as.factor(train_y)
test_y_factor <- as.factor(test_y)

# Create confusion matrices for the training and test sets
train_conf_matrix <- caret::confusionMatrix(train_predictions_factor, train_y_factor)
test_conf_matrix <- caret::confusionMatrix(test_predictions_factor, test_y_factor)

# Print confusion matrices
print("Training Confusion Matrix:")
print(train_conf_matrix)

print("Test Confusion Matrix:")
print(test_conf_matrix)
# Visualization function to show C, H, S instead of 0, 1, 2
visualize_confusion_matrix <- function(conf_matrix, title) {
  library(ggplot2)
  library(dplyr)
  
  # Extract the table from the confusion matrix
  cm_table <- as.data.frame(conf_matrix$table)
  
  # Map numeric values to diagnosis labels
  cm_table$Reference <- factor(cm_table$Reference, levels = c(1, 2, 3), 
                               labels = c("Healthy Control", "Bladder Stone", "Bladder Cancer"))
  cm_table$Prediction <- factor(cm_table$Prediction, levels = c(1, 2, 3), 
                                labels = c("Healthy Control", "Bladder Stone", "Bladder Cancer"))
  
  # Adjust text color for contrast
  cm_table <- cm_table %>%
    mutate(text_color = ifelse(Freq > max(Freq) / 2, "white", "black"))
  
  # Plot
  ggplot(cm_table, aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = Freq), color = "white") +
    scale_fill_gradient(low = "#A1D7E3", high = "#373F89") +
    geom_text(aes(label = Freq, color = text_color), 
              vjust = 0.5, fontface = "bold", size = 5, show.legend = FALSE) +
    scale_color_identity() +  # use text_color values directly
    labs(
      #title = title,
      x = "Actual Diagnosis",
      y = "Predicted Diagnosis",
      fill = "Count"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 10),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10)
    )
}
# Plot confusion matrices for training and test sets
train_cm_plot <- visualize_confusion_matrix(train_conf_matrix, "Training")
test_cm_plot <- visualize_confusion_matrix(test_conf_matrix, "Testing")

# Display the plots
print(train_cm_plot)
print(test_cm_plot)

# Save the plots to files if needed
ggsave("train_confusion_matrix_multiclass.png", train_cm_plot, width = 6, height = 5, units = "in", dpi = 300)
ggsave("train_confusion_matrix_multiclass.pdf", train_cm_plot, width = 6, height = 5, units = "in", dpi = 300)

ggsave("test_confusion_matrix_multiclass.png", test_cm_plot, width = 6, height = 5, units = "in", dpi = 300)
ggsave("test_confusion_matrix_multiclass.pdf", test_cm_plot, width = 6, height = 5, units = "in", dpi = 300)

# Predict the classes for training and test datasets
train_response <- predict(cvfit, newx = train_x, s = best_lambda, type = "response")
train_response<-as.data.frame(train_response)
test_response<- predict(cvfit, newx = test_x, s = best_lambda, type = "response")
test_response<-as.data.frame(test_response)

data_response<-rbind(train_response,test_response)
data_response$Patient<-sub(paste0("*", "_X"), "", rownames(data_response))
library(readxl)
BC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BC")
BC_cd<-as.data.frame(cbind(BC_clinicial_data$Patient,BC_clinicial_data$MIBC,BC_clinicial_data$AGE,BC_clinicial_data$SEX,
                           BC_clinicial_data$Tumor_size,BC_clinicial_data$`T`,BC_clinicial_data$Tumor_number,BC_clinicial_data$Grade,"BC"))
colnames(BC_cd)<-c("Patient","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade','group')

BS_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "BS")
BS_cd<-as.data.frame(cbind(BS_clinicial_data$Patient,NA,BS_clinicial_data$AGE,BS_clinicial_data$SEX,
                           NA,NA,NA,NA,"BS"))
colnames(BS_cd)<-c("Patient","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade','group')

HC_clinicial_data <- read_excel("/dssg/home/acct-dahan/share/cfRNA/RESULT/ANALYSIS/clinical/1109_clinical data.xlsx", sheet = "HC")
HC_cd<-as.data.frame(cbind(HC_clinicial_data$Patient,NA,HC_clinicial_data$AGE,HC_clinicial_data$SEX,
                           NA,NA,NA,NA,"HC"))
colnames(HC_cd)<-c("Patient","MIBC","Age","Sex","Tumor_size","StageT",'Tumor_number','Grade','group')

clinic_data<-rbind(BC_cd,BS_cd,HC_cd)

data_merge<-merge(data_response,clinic_data,by='Patient',all.x=T)

head(data_merge)

data_merge$StageT[grep('1',data_merge$StageT)]<-'BC_T1'
data_merge$StageT[grep('2',data_merge$StageT)]<-'BC_T2'
data_merge$StageT[grep('3',data_merge$StageT)]<-'BC_T3'
data_merge$StageT[grep('4',data_merge$StageT)]<-'BC_T4'
data_merge$StageT[grep('a',data_merge$StageT)]<-'BC_Ta'
data_merge$StageT[data_merge$group=="HC"] <- 'HC'
data_merge$StageT[data_merge$group=="BS"] <- 'BS'

library(plotly) 
colnames(data_merge)[2:4]<-c("Prob_HC","Prob_BS","Prob_BC")
data_merge$group <- factor(data_merge$group)

plot_ly(data_merge, 
        x = ~Prob_HC, 
        y = ~Prob_BS, 
        z = ~Prob_BC, 
        color = ~group, 
        colors = c(#"BC_Ta" = "#FBD0DF", 
          #"BC_T1" = "#F3ADCA",
          #"BC_T2" = "#EB8CB6",
          #"BC_T3" = "#E36BA2",
          "BC" = "#DB498E",
          "HC" = "#379FB4",
          "BS" = "#FCAE59"),
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 2),
        text = ~paste("Patient:", Patient, "<br>group:", group)) %>%
  layout(scene = list(
    xaxis = list(title = "Prob_HC"),
    yaxis = list(title = "Prob_BS"),
    zaxis = list(title = "Prob_BC")
  ))

group_colors <- c("BC" = "#DB498E", "HC" = "#379FB4", "BS" = "#FCAE59")

col <- group_colors[as.character(data_merge$group)]
x <- data_merge$Prob_HC
y <- data_merge$Prob_BS
z <- data_merge$Prob_BC

library(Cairo)
CairoPDF("X_repeataware_scatter3D_plot.pdf", width = 7, height = 6,family='Arial')
#CairoPNG("X_repeataware_scatter3D_plot.png", width = 7, height = 6, units = "in", res = 300, family = "Arial")
scatter3D(x, y, z,
          col = col,
          colvar = NULL,   #
          pch = 19,
          cex = 0.6,
          cex.lab = 0.6,
          cex.axis = 0.6,
          xlab = "Prob_HC",
          ylab = "Prob_BS",
          zlab = "Prob_BC",
          theta = 135,
          phi = 45,
          xlim = c(0, 1),
          ylim = c(0, 1),
          zlim = c(0, 1),
          #main = "3D Scatter Plot",
          ticktype = "detailed",
          bty = "g",
          colkey = FALSE
)
legend("topright",                   
       legend = names(group_colors),
       col = group_colors,           
       pch = 19,                    
       cex = 1,                     
       bty = "n")                   

dev.off()

data_merge$StageT <- factor(data_merge$StageT, levels = c("HC","BS","BC_Ta", "BC_T1", "BC_T2", "BC_T3", "BC_T4"), ordered = TRUE)


library(ggpubr)
library(ggplot2)
library(ggpubr)

p2<-ggplot(data_merge, aes(x = StageT, y = Prob_BC, fill = StageT)) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "HC" = "#379BF4",
    "BS" = "#FCAE59",
    "BC_Ta" = "#FBD0DF", 
    "BC_T1" = "#F3ADCA",
    "BC_T2" = "#EB8CB6",
    "BC_T3" = "#E36BA2",
    "BC_T4" = "#DB498E"
  )) +
  labs(
    x = "Patient_stage",
    y = "Probability of BC"
  ) +  scale_y_continuous(
    #limits = c(0, 1.2),
    breaks = seq(0, 1, by = 0.2)
  )+
  stat_compare_means(
    method = "wilcox.test",  
    label = "p.signif",
    comparisons = list(
      c("HC", "BS"),
      c("HC", "BC_Ta"),
      c("BS", "BC_Ta")
    ),
    hide.ns = TRUE,          
    paired = FALSE,
    label.size = 6,
    tip.length = 0.01,
    step.increase = 0.05,
    label.y = c(0.95, 0.98, 1.01)
  )+
  theme(plot.title = element_text(size = 8, hjust = 0.5))+
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 0)),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.x = unit(0.5, 'cm'),   
    legend.spacing.y = unit(1, 'cm'),   
    legend.key.spacing.y = unit(0.1, "cm"),    
    legend.key.height = unit(0.5, "cm"),    
    legend.key.width = unit(0.5, "cm"),    
    legend.box.spacing = unit(0.5, 'cm'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),) +
  theme(legend.position = "right")
p2
ggsave("X_repeataware_Stage.png",height=8,width=12,dpi=300)
ggsave("X_repeataware_Stage.pdf",height=8,width=12,dpi=300)