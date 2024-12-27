#########################proteins############################
rm(list = ls())
data <- read.csv(file = "select_pro_agesex.csv", row.names = 1) #don't , check.names = FALSE here
data <- data[,-c(2,3)] #delete sex and age colums
# change "male"and "female"to numeric
#data$sex <- ifelse(data$sex == "male", 1, 0)
data$BMI <- as.factor(data$BMI)
set.seed(1009)
ind<-sample(2,nrow(data),replace = TRUE,prob = c(0.7,0.3))
traindata<-data[ind == 1,]
testdata<-data[ind == 2,]


set.seed(1009)
#default is 500 trees
randomForest_rf<-randomForest(BMI ~.,data = traindata, #BMI ~ . - age - sex
                              ntree = 500,
                              mtry = 4,
                              proximity = TRUE,
                              importance = TRUE)
#model validation
set.seed(1009)
irispred <-predict(randomForest_rf,newdata = testdata)
freq1 <- table(irispred,testdata$BMI)
#plot(margin(randomForest_rf,testdata$BMI))
#mixed matrix
caret::confusionMatrix(testdata$BMI, irispred) #AUC=0.7841

irispred1 <- predict(randomForest_rf,newdata = testdata, type = "pro")
rocc1 <- roc(testdata$BMI, as.matrix(irispred1[ , 1]),ci=T,auc=T, direction=">")
rocc1
library(e1071)
set.seed(1009)
svm <- svm(BMI ~ . , 
           data = traindata,
           kernel = 'radial',C = 1, 
           sigma = 0.067, probability = TRUE) #type ="C-classification"
# prediction
pre_svm <- predict(svm, newdata = testdata, probability = TRUE)
table(pre_svm,testdata$BMI)
#mixed matrix
caret::confusionMatrix(testdata$BMI, pre_svm) 
#roc
#predict can output different typrs, such as "decision.values",
#"probabilities"
pre_svm1 <- attr(pre_svm, "probabilities")
rocc2 <- roc(testdata$BMI,as.matrix(pre_svm1[ , 1]),ci=T,auc=T, direction=">")
rocc2
set.seed(1009)
pre <- kknn(BMI ~ . ,traindata, testdata,
            distance = 2, k =9, 
            kernel= "triangular")
summary(pre)
fit <- predict(pre,newdata = testdata)
t<-table(testdata$BMI, fit)
sum(diag(t))/sum(t)
caret::confusionMatrix(testdata$BMI, fit)
fit1 <- predict(pre,newdata = testdata, type = "prob")
head(fit1)
rocc3 <- roc(testdata$BMI, as.matrix(fit1[ , 1]),ci=T,auc=T, direction=">")
rocc3
auc(rocc3)
set.seed(1009)
library(adabag)
ada <-boosting(BMI ~.,data= traindata,mfinal=150, 
               control=rpart.control(maxdepth=3))
pred<-predict.boosting(ada,newdata=testdata, type = "class")
pred$confusion
ad <- as.factor(pred$class)
#mixed matrix
caret::confusionMatrix(testdata$BMI, ad)
rocc4 <- roc(testdata$BMI, data.frame(pred$prob)[,-1],ci=T,auc=T, direction="<")
rocc4

#output result
# Create a vector of prediction results
predictions <- list("irispred", "pre_svm", "fit", "ad")

# Create an empty dataframe to store results
result <- data.frame(
  Value = numeric(),
  Model = character()
)
#AccuracyNull is  NoInformationRate
# Loop through each prediction
for (pred in predictions) {
  # Get confusion matrix
  conf_matrix <- confusionMatrix(testdata$BMI, get(pred))
  
  # Extract metrics
  accuracy <- conf_matrix$overall['Accuracy']
  kappa <- conf_matrix$overall['Kappa']
  sensitivity <- conf_matrix$byClass['Sensitivity']
  specificity <- conf_matrix$byClass['Specificity']
  ppv <- conf_matrix$byClass['Pos Pred Value']
  npv <- conf_matrix$byClass['Neg Pred Value']
  f1 <- conf_matrix$byClass['F1']
  NoInformationRate <- conf_matrix$overall['AccuracyNull']
  
  # Add results to dataframe
  result <- rbind(result, data.frame(
    Value = c( sensitivity, specificity, NoInformationRate, accuracy, kappa,ppv, npv, f1),
    Model = rep(pred, 8)
  ))
}

result$Value <- signif(result$Value, digits = 3)
# Transpose the result dataframe
result <- t(result)

##########################phenotype#######################
rm(list = ls())
data <- read.csv(file = "pheno_11_pro_809.csv", row.names = 1)
data <- data[,-c(1,12:820)] #delete sex and age colums
# change "male"and "female"to numeric
#data$sex <- ifelse(data$sex == "male", 1, 0)
data$BMI_grade <- as.factor(data$BMI)
colnames(data)[1] <- "BMI"
set.seed(1009)
ind<-sample(2,nrow(data),replace = TRUE,prob = c(0.7,0.3))
traindata<-data[ind == 1,]
testdata<-data[ind == 2,]


set.seed(1009)
#default is 500 trees
randomForest_rf<-randomForest(BMI ~.,data = traindata, #BMI ~ . - age - sex
                              ntree = 500,
                              mtry = 2,
                              proximity = TRUE,
                              importance = TRUE)
table(predict(randomForest_rf),traindata$BMI)
randomForest_rf
#Plot error rates based on different numbers in the random forest
#plot(randomForest_rf)
#importance score
importance(randomForest_rf,type = 1)
#importance plot
varImpPlot(randomForest_rf)

#model validation
irispred<-predict(randomForest_rf,newdata = testdata)
freq1 <- table(irispred,testdata$BMI)
#plot(margin(randomForest_rf,testdata$BMI))
#mixed matrix
caret::confusionMatrix(testdata$BMI, irispred) #AUC=0.7841
set.seed(1009)
irispred1 <- predict(randomForest_rf,newdata = testdata, type = "pro")
rocc5 <- roc(testdata$BMI, as.matrix(irispred1[ , 1]),ci=T,auc=T, direction=">")
rocc5

library(e1071)
set.seed(1009)
svm <- svm(BMI ~ . , 
           data = traindata,
           kernel = 'radial',C = 1, 
           sigma = 0.111, probability = TRUE) #type ="C-classification"
# prediction
pre_svm <- predict(svm, newdata = testdata, probability = TRUE)
table(pre_svm,testdata$BMI)
#mixed matrix
caret::confusionMatrix(testdata$BMI, pre_svm) 
#roc
#predict can output different typrs, such as "decision.values",
#"probabilities"
pre_svm1 <- attr(pre_svm, "probabilities")
rocc6 <- roc(testdata$BMI,as.matrix(pre_svm1[ , 1]),ci=T,auc=T, direction=">")
rocc6
auc(rocc6)

library(caret)
library(kknn)
set.seed(1009)
pre <- kknn(BMI ~ . ,traindata, testdata,
            distance = 2, k =9, 
            kernel= "triangular")
summary(pre)
fit <- predict(pre,newdata = testdata)
t<-table(testdata$BMI, fit)
sum(diag(t))/sum(t)
caret::confusionMatrix(testdata$BMI, fit)
fit1 <- predict(pre,newdata = testdata, type = "prob")
head(fit1)
rocc7 <- roc(testdata$BMI, as.matrix(fit1[ , 1]),ci=T,auc=T, direction=">")
rocc7
auc(rocc7)

set.seed(1009)
library(adabag)
ada <-boosting(BMI ~.,data= traindata,mfinal=50, 
               control=rpart.control(maxdepth=3))
pred<-predict.boosting(ada,newdata=testdata, type = "class")
pred$confusion

ad <- as.factor(pred$class)

#mixed matrix
caret::confusionMatrix(testdata$BMI, ad)
rocc8 <- roc(testdata$BMI, data.frame(pred$prob)[,-1],ci=T,auc=T, direction="<")
rocc8

#output result
# Create a vector of prediction results
predictions <- list("irispred", "pre_svm", "fit", "ad")

# Create an empty dataframe to store results
result <- data.frame(
  Value = numeric(),
  Model = character()
)
#AccuracyNull is  NoInformationRate
# Loop through each prediction
for (pred in predictions) {
  # Get confusion matrix
  conf_matrix <- confusionMatrix(testdata$BMI, get(pred))
  
  # Extract metrics
  accuracy <- conf_matrix$overall['Accuracy']
  kappa <- conf_matrix$overall['Kappa']
  sensitivity <- conf_matrix$byClass['Sensitivity']
  specificity <- conf_matrix$byClass['Specificity']
  ppv <- conf_matrix$byClass['Pos Pred Value']
  npv <- conf_matrix$byClass['Neg Pred Value']
  f1 <- conf_matrix$byClass['F1']
  NoInformationRate <- conf_matrix$overall['AccuracyNull']
  
  # Add results to dataframe
  result <- rbind(result, data.frame(
    Value = c( sensitivity, specificity, NoInformationRate, accuracy, kappa,ppv, npv, f1),
    Model = rep(pred, 8)
  ))
}

result$Value <- signif(result$Value, digits = 3)
# Transpose the result dataframe
result <- t(result)
#############################extract AUC#####################
res <- data.frame(
  Model = c("rocc1", "rocc2", "rocc3", "rocc4"),#"rocc5", "rocc6", "rocc7", "rocc8"
  #"rocc1", "rocc2", "rocc3", "rocc4"
  AUC = numeric(4),
  CI_lower = numeric(4),
  CI_upper = numeric(4)
)

rocc_list <- list(rocc1, rocc2, rocc3, rocc4)#rocc5, rocc6, rocc7, rocc8
#rocc1, rocc2, rocc3, rocc4
for (i in 1:length(rocc_list)) {
  auc_value <- auc(rocc_list[[i]])
  ci <- ci(rocc_list[[i]])
  res$AUC[i] <- auc_value
  res$CI_lower[i] <- ci[1]
  res$CI_upper[i] <- ci[3]
}

res$AUC <- sprintf("%.3f", res$AUC)
res$CI_lower <- sprintf("%.3f", res$CI_lower)
res$CI_upper <- sprintf("%.3f", res$CI_upper)
res_with_format <- paste(res$AUC, " (", res$CI_lower, "-", res$CI_upper, ")", sep="")

res_table <- data.frame(Model = res$Model, AUC_CI = res_with_format)
######################plot##############
dev.new()
#png("plot/mc_auc.png", width=60, height=30, units = "cm", res=600)
par(mfrow = c(2, 4), cex.lab = 1.5)  

plot(rocc8, print.auc = TRUE,
     print.auc.adj = (0.05), 
     print.auc.x = 0.63,
     print.auc.y = 0.3,
     print.auc.col = "#E72929",
     print.auc.cex =2, 
     auc.polygon = TRUE,
     max.auc.polygon = TRUE,
     auc.polygon.col = "#F0F3FF",
     grid = c(0.1, 0.2),
     cex.lab = 2, ci.col = "red",
     grid.col = c("green", "red"))


dev.off()
