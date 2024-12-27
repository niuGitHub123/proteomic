rm(list = ls())
dev.new()
png("plot/mc.png", width=80, height=40, units = "cm", res=600) 
par(mfrow = c(2, 4))
mar = c(10, 8, 8, 8)
#prepare and devide data (7:3)
setwd("C:/project/pro/logistic")
data <- read.csv(file = "select_pro_agesex.csv", row.names = 1)
data <- data[,-c(2,3)] #delete sex and age colums
# change "male"and "female"to numeric
#data$sex <- ifelse(data$sex == "male", 1, 0)
set.seed(1009)
data$BMI <- as.factor(data$BMI)
ind<-sample(2,nrow(data),replace = TRUE,prob = c(0.7,0.3))
traindata<-data[ind == 1,]
testdata<-data[ind == 2,]

#random forest
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9789589/
#use 10 fold cross validation to choose the best parameters
library(caret)
library(randomForest)
library(pROC)
set.seed(1009)
#default is 500 trees
randomForest_rf<-randomForest(BMI ~.,data = traindata, #BMI ~ . - age - sex
                              ntree = 500,
                              mtry = 4,
                              proximity = TRUE,
                              importance = TRUE)
table(predict(randomForest_rf),traindata$BMI)
randomForest_rf

#model validation
irispred<-predict(randomForest_rf,newdata = testdata)
freq1 <- table(irispred,testdata$BMI)
#mixed matrix
caret::confusionMatrix(testdata$BMI, irispred) #AUC=0.7841

irispred <- predict(randomForest_rf,newdata = testdata, type = "pro")
head(irispred)
rocc <- roc(testdata$BMI, as.matrix(irispred[ , 1]),ci=T,auc=T, direction=">")
rocc

plot(rocc, print.auc = TRUE,
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

#SVM
#chose parameters
set.seed(1009)
library(e1071)
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
pre_svm <- predict(svm, newdata = testdata, probability = TRUE)
#predict can output different typrs, such as "decision.values",
#"probabilities"
pre_svm <- attr(pre_svm, "probabilities")
rocc <- roc(testdata$BMI,as.matrix(pre_svm[ , 1]),ci=T,auc=T, direction="<")
rocc
auc(rocc)

plot(rocc, print.auc = TRUE,
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

#KNN
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
fit <- predict(pre,newdata = testdata, type = "prob")
head(fit)
rocc <- roc(testdata$BMI, as.matrix(fit[ , 1]),ci=T,auc=T, direction=">")
rocc
auc(rocc)

plot(rocc, print.auc = TRUE,
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


#adaboost
set.seed(1009)
library(adabag)
ada <-boosting(BMI ~.,data= traindata,mfinal=50, 
               control=rpart.control(maxdepth=3))
pred<-predict.boosting(ada,newdata=testdata, type = "class")
pred$confusion

#mixed matrix
caret::confusionMatrix(testdata$BMI, as.factor(pred$class))
rocc <- roc(testdata$BMI, data.frame(pred$prob)[,-1],ci=T,auc=T, direction="<")
rocc
plot(rocc, print.auc = TRUE,
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
     grid.col = c("green", "red")) #,print.thres = TRUE


##################phenotype machine learning###################
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

#random forest
set.seed(1009)
# include ntree, mtry(waiting to be choosen), run two times
#default is 500 trees
randomForest_rf<-randomForest(BMI ~.,data = traindata, #BMI ~ . - age - sex
                              ntree = 500,
                              mtry = 2,
                              proximity = TRUE,
                              importance = TRUE)
table(predict(randomForest_rf),traindata$BMI)
randomForest_rf

#model validation
irispred<-predict(randomForest_rf,newdata = testdata)
freq1 <- table(irispred,testdata$BMI)
#mixed matrix
caret::confusionMatrix(testdata$BMI, irispred) #AUC=0.7841

irispred <- predict(randomForest_rf,newdata = testdata, type = "pro")
head(irispred)
rocc <- roc(testdata$BMI, as.matrix(irispred[ , 1]),ci=T,auc=T, direction=">")
rocc

plot(rocc, print.auc = TRUE,
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


#SVM
#chose parameters
set.seed(1009)
library(e1071)
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
pre_svm <- predict(svm, newdata = testdata, probability = TRUE)
#predict can output different typrs, such as "decision.values",
#"probabilities"
pre_svm <- attr(pre_svm, "probabilities")
rocc <- roc(testdata$BMI,as.matrix(pre_svm[ , 1]),ci=T,auc=T, direction="<")
rocc
auc(rocc)

plot(rocc, print.auc = TRUE,
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

#KNN
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
fit <- predict(pre,newdata = testdata, type = "prob")
head(fit)
rocc <- roc(testdata$BMI, as.matrix(fit[ , 1]),ci=T,auc=T, direction=">")
rocc
auc(rocc)

plot(rocc, print.auc = TRUE,
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


#adaboost
library(adabag)
set.seed(1009)
ada <-boosting(BMI ~.,data= traindata,mfinal=50, 
               control=rpart.control(maxdepth=3))
pred<-predict.boosting(ada,newdata=testdata, type = "class")
pred$confusion

#mixed matrix
caret::confusionMatrix(testdata$BMI, as.factor(pred$class))
rocc <- roc(testdata$BMI, data.frame(pred$prob)[,-1],ci=T,auc=T, direction="<")
rocc
plot(rocc, print.auc = TRUE,
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
     grid.col = c("green", "red")) #,print.thres = TRUE
dev.off()
