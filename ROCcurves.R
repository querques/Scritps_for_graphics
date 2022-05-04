install.packages("ROCR")
library(ROCR)
Base_respi_ordenada$Respiratorysymptoms <- as.numeric(as.character(Base_respi_ordenada$Respiratorysymptoms))
Base_respi_ordenada$sitting_FVC <- as.numeric(as.character(Base_respi_ordenada$sitting_FVC))
Base_respi_ordenada$Ventilatorysupport <- as.numeric(as.character(Base_respi_ordenada$Ventilatorysupport))
pred <- prediction(Base_respi_ordenada$sitting_FVC, Base_respi_ordenada$Ventilatorysupport)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE)
#predicting now ROC curve for respiratory symptoms
pred <- prediction(Base_respi_ordenada$sitting_FVC, Base_respi_ordenada$Respiratorysymptoms)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE)
#let´s try a different approach with pROC that provides AUC value
library(pROC)
pROC_obj <- roc(Base_respi_ordenada$Respiratorysymptoms,Base_respi_ordenada$sitting_FVC, smoothed = TRUE,ci=TRUE, ci.alpha=0.9, stratified=FALSE,plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE, show.thres=TRUE)
sens.ci <- ci.se(pROC_obj)
plot(sens.ci, type="shape", col="lightblue")
plot(sens.ci, type="bars")
#predicting new ROC curve for Ventilatory support
pROC_obj <- roc(Base_respi_ordenada$Ventilatorysupport,Base_respi_ordenada$sitting_FVC, smoothed = TRUE,ci=TRUE, ci.alpha=0.9, stratified=FALSE,plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE, show.thres=TRUE)
sens.ci <- ci.se(pROC_obj)
plot(sens.ci, type="shape", col="lightblue")
plot(sens.ci, type="bars")
#let´s try a different approach with plotROC
install.packages("plotROC")
library(plotROC)
library(gplots)
rocplot <- ggplot(Base_respi_ordenada, aes(m = sitting_FVC, d = Ventilatorysupport))+ geom_roc(n.cuts=20,labels=FALSE)
rocplot + style_roc(theme = theme_grey) + geom_rocci(fill="pink") 
#let´s try a different approach with ROCit thta provides the Youden point
install.packages("ROCit")
library(ROCit)
ROCit_obj <- rocit(score=Base_respi_ordenada$sitting_FVC,class=Base_respi_ordenada$Ventilatorysupport)
plot(ROCit_obj)
#predicting new ROC curve for Respiratory symptoms
ROCit_obj <- rocit(score=Base_respi_ordenada$sitting_FVC,class=Base_respi_ordenada$Respiratorysymptoms)
plot(ROCit_obj)
ksplot(ROCit_obj)

