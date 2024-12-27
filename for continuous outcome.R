options(scipen = 1000000000)
options(stringsAsFactors = F)
set.seed(1234)
library(readxl)
library(dplyr)
library(corrplot)
library(corrr)
library(RColorBrewer)
library(psych)
library(GGally)

library(CBCgrps)
library(tibble)
library(ggrepel)

library(ggplot2)
library(hrbrthemes)

setwd("C:/project/pro/analyses")
#################phenotype preparation##############
#load proteomics phenotype 2132
pheno_pro_2132 <- read.csv("pheno_pro_2132.csv")
#1884 participants have both trans and pros 
phenotype_co<- read.csv("phenotype_co.csv")

#subjects filter and preparation for analysis
pheno_all <- phenotype_co[, c("u3n_expr_ff4" , "u3n_prot_ff4", "ID.Links", "zz_nr",
                              "lg_dnaAxiom_s4f4","u3n_meth850K_ff4",
                              #all IDs
                              "u3n_expr_batch_ff4", #batch information
                              "u3talteru", "u3csex", "u3tgewi", "u3tbmi",
                              "u3tgroe","u3tbmiwho","u3talkkon","u3ttumf","u3twhrat",
                              "u3tdm110y15","u3tglukfastn","u3tsysmm","u3tnuecht",
                              "u3tdiamm", "u3tcigsmk","u3tphys","u3lk_hdln",
                              "u3lk_ldln","u3lk_trin","u3lk_choln",
                              "u3h_crp", "u3tbfp_k")]
sum(is.na(pheno_all))
clinic <- na.omit(pheno_all)

clinic$u3csex <- factor(clinic$u3csex,
                        levels = c(1,2),
                        labels = c("male", "female"))
clinic$u3tbmiwho <- factor(clinic$u3tbmiwho,
                           levels = c(1,2,3,4,5,6),
                           labels = c("underweight", "normal weight", "pre-obese","obese grade 1", "obese grade 2","obese grade 3"))
clinic$u3tdm110y15 <- factor(clinic$u3tdm110y15,
                             levels = c(0,1,2,3,4,5,6,7,8,9),
                             labels = c("normal","IFG","IGT","IFG/IGT","new diagnosed", "known Typ 2 Diabetes", "Typ 1 Diabetes",  "known drug",  "unknown", "no validation"))
table(clinic$u3tdm110y15)
clinic$u3tcigsmk <- factor(clinic$u3tcigsmk,
                           levels = c(1,2,3),
                           labels = c("smoker", "ex-smoker", "never-smoker"))
clinic$u3tphys <- factor(clinic$u3tphys, levels = c(1,2),
                         labels = c("active", "inactive"))

#firstly filter underweight (BMI < 15 kg/m2) or all missing covariate values
summary(clinic$u3tbmiwho)

clinic <- filter(clinic, clinic$u3tbmiwho != "underweight")
summary(clinic)

#filer didn’t have at least 8 h of fasting 
table(clinic$u3tnuecht)
clinic <- filter(clinic, clinic$u3tnuecht != 2)

#filer out diabetes
summary(clinic$u3tdm110y15)
clinic$u3tdm110y15 <- factor(clinic$u3tdm110y15,
                             levels = c("normal","IFG","IGT","IFG/IGT","new diagnosed", "known Typ 2 Diabetes", "Typ 1 Diabetes",  "known drug",  "unknown", "no validation"),
                             labels = c("Non_T2D", "Non_T2D", "Non_T2D", "Non_T2D", "T2D", "T2D", "Non_T2D", "Non_T2D", "Non_T2D", "Non_T2D"))
clinic$u3tbmiwho <- factor(clinic$u3tbmiwho,
                           levels = c("normal weight", "pre-obese", "obese grade 1", "obese grade 2", "obese grade 3"),
                           labels = c("Non_Obese", "Non_Obese", "Obese", "Obese", "Obese"))
summary(clinic[ , 8:29])
write.csv(clinic, file = "phenotype_final.csv", row.names = FALSE)

###########################stastical analysis####################
clinic <- read.csv("phenotype_final.csv", header = 1, row.names = 1)
continous <- clinic[,c("u3tbmiwho","u3talteru", "u3tgewi", "u3tgroe","u3talkkon",
                       "u3ttumf","u3twhrat", "u3tglukfastn","u3tsysmm",
                       "u3tdiamm", "u3lk_hdln", "u3lk_ldln","u3lk_trin",
                       "u3lk_choln", "u3h_crp", "u3tbfp_k")]
category <- clinic[,c("u3tbmiwho", "u3csex", "u3tcigsmk", "u3tphys",
                      "u3tdm110y15")]
table1 <- twogrps(continous, gvar = "u3tbmiwho", norm.rd = 1, sk.rd =1 , cat.rd =1 ,
                  p.rd =3,ShowStatistic = T,  pnormtest = 0)
print(table1$Table, quote = T)
table2 <- twogrps(category, gvar = "u3tbmiwho", norm.rd = 1, sk.rd =1 , cat.rd =1 ,
                  p.rd =3,ShowStatistic = F,  pnormtest = 0)
print(table2$Table, quote = T)

write.csv(table1$Table, file = "result/continous_var_stastic.csv")
write.csv(table2$Table, file = "result/category_var_stastic.csv")

################linear regression model######################
rm(list = ls())
clinic <- read.csv("phenotype_final.csv", header = 1)
proteomic_2132 <- read.csv("proteomic_2132.csv", header = 1)
#adjust/match protein names and numbers
offical_name <- read_excel("offical_name.xlsx")
pro_name <- offical_name[, c("koraID", "Genenames")]
pro_name <- pro_name[complete.cases(pro_name), ]
proteomic_3 <- proteomic_2132[, 1:3]
selected_columns <- names(proteomic_2132[,4:845]) %in% pro_name$koraID
proteomic_809 <- proteomic_2132[, selected_columns]
colnames(proteomic_809)[4:812] <- pro_name$Genenames

com_data <- merge(clinic, proteomic_809, by = "zz_nr")
write.csv(com_data, file = "pheno_pro_1557.csv", row.names = FALSE)
write.csv(proteomic_809, file = "proteomic_809.csv", row.names = FALSE)

clinic_linear <- data.frame(row.names = com_data$zz_nr,
                            BMI = com_data$u3tbmi,
                            BMI_grade= com_data$u3tbmiwho,
                            age = com_data$u3talteru,
                            sex = com_data$u3csex,
                            activities = com_data$u3tphys,
                            SBP = com_data$u3tsysmm,
                            TG = log(com_data$u3lk_trin),
                            HDL = com_data$u3lk_hdln,
                            smoking = com_data$u3tcigsmk,
                            FPG = com_data$u3tglukfastn,
                            T2D = com_data$u3tdm110y15)
sum(is.na(clinic_linear))
co_data <- cbind(clinic_linear, com_data[, 32:840])
lm(BMI ~ + age + sex + activities +
     SBP + TG + smoking + HDL + FPG + T2D, data = clinic_linear) 

#basic analysis
set.seed(123)
outcome1 <- c()
for (i in 12:820){
  model_linear <- lm(BMI ~ co_data[,i] + age + sex, data = co_data) 
  outcome1 <- rbind(outcome1, c(colnames(co_data)[i], coef(summary(model_linear))[2, c(1,2,4)]))
}

outcome1
outcome1 <- data.frame(outcome1)
names(outcome1)[4]  <- "pvalue"
names(outcome1)[1] <- "Proteins"
for (i in 2:4) {outcome1[,i] <- as.numeric(outcome1[,i])}
outcome1 <- outcome1[order(outcome1$pvalue),]
outcome1$P.adjust <- p.adjust(outcome1$pvalue,method="bonferroni",n=length(outcome1$pvalue))

# export the outcome1
outcome1[,2:3] <- round(outcome1[,2:3], digits = 3)
outcome1$pvalue <- signif(outcome1$pvalue, digits = 3)
outcome1$P.adjust <- signif(outcome1$P.adjust, digits = 3)
outcome1$log10_adj_pval = -log10(outcome1$P.adjust)
outcome1$log10_pval = -log10(outcome1$pvalue)
#table(outcome1[,5] < 0.05)
outcome1_sig <- outcome1[outcome1$P.adjust < 0.05, ]

#save it
write.csv(outcome1, file = "linear_basic.csv", quote = FALSE, row.names = FALSE)
write.csv(outcome1_sig, file = "linear_basic_sig.csv", quote = FALSE, row.names = FALSE)

#sensitive analysis
set.seed(123)
outcome <- c()
for (i in 12:820){
  model_linear <- lm(BMI ~ co_data[,i] + age + sex +
                       SBP + FPG + T2D + TG  + HDL + smoking + activities, data = co_data) 
  outcome <- rbind(outcome, c(colnames(co_data)[i], coef(summary(model_linear))[2, c(1,2,4)]))
}

outcome
outcome <- data.frame(outcome)
names(outcome)[4]  <- "pvalue"
names(outcome)[1] <- "Proteins"
for (i in 2:4) {outcome[,i] <- as.numeric(outcome[,i])}
outcome <- outcome[order(outcome$pvalue),]
outcome$P.adjust <- p.adjust(outcome$pvalue,method="bonferroni",n=length(outcome$pvalue))

# export the outcome
outcome[,2:3] <- round(outcome[,2:3], digits = 3)
outcome$pvalue <- signif(outcome$pvalue, digits = 3)
outcome$P.adjust <- signif(outcome$P.adjust, digits = 3)
outcome$log10_adj_pval = -log10(outcome$P.adjust)
outcome$log10_pval = -log10(outcome$pvalue)
#table(outcome[,5] < 0.05)
outcome_sig <- outcome[outcome$P.adjust < 0.05, ]

#save it
write.csv(outcome, file = "linear_full.csv", quote = FALSE, row.names = FALSE)
write.csv(outcome_sig, file = "linear_full_sig.csv", quote = FALSE, row.names = FALSE)

sum(outcome_sig$Proteins %in% outcome1_sig$Proteins)
#plot
outcome <- read.csv("linear_basic.csv", header=TRUE, row.names=1)
table(outcome[,4] < 0.05)
#filtered_outcome <- outcome1[outcome1$P.adjust < 0.05, ]
plot_log = data.frame(row.names = rownames(outcome),
                      Estimate = outcome$Estimate, 
                      log10_pval = outcome$log10_pval)
#filter upregulated and down regulated proteins and use ggrepel function to put it in volcano plot
sig_linear <- outcome[outcome[,4] < 0.05, ]
sig_linear$labels <- rownames(sig_linear)
head(sig_linear,20)
relative_size <- (outcome$log10_pval)/3

#x:c(-1.1, 1.1),y: c(-1,19) full
#x:c(-1.8, 1.7),y: c(-1,35) basic
linear_full <- ggplot(plot_log, aes(x = Estimate, y=log10_pval,
                                    color = ifelse(outcome$P.adjust < 0.05, "Association",
                                                   "Unassociation")))+
  geom_point(shape = 16, alpha = 0.7, size = 3.5) + 
  scale_x_continuous(limits=c(-1.8, 1.7), expand=c(0,0)) +
  scale_y_continuous(limits=c(-1,35), expand=c(0,0)) + 
  labs(x= "Beta coefficient", y=bquote('P ('*-Log[10]*')'),
       color = "Association Status") +
  ggtitle("BMI association") +
  geom_hline(yintercept = -log10(0.05/842), linetype="dashed")+
  theme_classic() +
  scale_color_manual(values=c("Association"="red", "Unassociation" = "grey")) + 
  geom_text_repel(data = sig_linear,
                  aes(label = labels), color = "black",
                  size = 3, family = "Arial", force = 1,
                  segment.color = "black", segment.size = 0.8, segment.alpha = 1,
                  hjust = 0.5, show.legend=F)+ 
  theme(legend.position = "none")#delete legend
linear_full
ggsave("plot/linear_basic.png", linear_full, dpi=600, width= 15, height= 15, unit="cm")

###################model type ################################
a <- data.frame(
  x = factor(c("+age,sex(BMI)", "+BP", "+FPG,T2D", "+lipids", "+lifestyle"), 
             levels = c("+age,sex(BMI)", "+BP", "+FPG,T2D", "+lipids", "+lifestyle")),
  y = as.numeric(c(33, 30, 28, 17, 16))
)

model_type <- ggplot(a, aes(x = x, y = y, group = 1)) +
  geom_line( color="#FF8F8F",size= 2, alpha = 1) +
  geom_point(shape=21, color="black", fill="#FF8F8F", size=5) +
  theme_ipsum()+
  labs(x = "Model Type", y = "Number of Significant Proteins")+
  theme_bw() +
  geom_text(aes(label = y), vjust = -0.8, hjust = -0, size = 4)+ 
  theme(axis.title.x = element_text(size = 14),  # 调整 X 轴标签字体大小
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),  # 调整 X 轴标签字体大小
        axis.text.y = element_text(size = 12)) 
model_type
ggsave("plot/model_type.png", model_type, dpi=600, width= 15, height= 15, unit="cm")

#################priority lasso###############################
#https://cran.r-project.org/web/packages/prioritylasso/vignettes/prioritylasso_vignette.html
library(prioritylasso)
library(pROC)
library(dplyr)

linear_basic_sig <- read.csv("linear_basic_sig.csv", header=TRUE, row.names=1)
pheno_pro_1557 <- read.csv("pheno_pro_1557.csv")
pro <- pheno_pro_1557[,32:840]
clinic<- data.frame(row.names = pheno_pro_1557$zz_nr,
                    BMI = pheno_pro_1557$u3tbmi,
                    age = pheno_pro_1557$u3talteru,
                    sex = pheno_pro_1557$u3csex)
data<- cbind(clinic, pro[,colnames(pro) %in% rownames(linear_basic_sig)])
# change "male"and "female"to numeric
data$sex <- ifelse(data$sex == "male", 1, 0)

data <- as.matrix(data)
pl_out <- as.numeric(data[,1])  # outcome
pl_pred <- data[,2:35]  # age+sex+ proteins
blocks <- list(bp1=1:2, bp2=3:34) #in 2:25, 1:2 is age+sex, 3:24 is proteins

#left:lambda.min, right:lambda.1se
pl1<- prioritylasso(X = pl_pred, Y = pl_out, family = "gaussian",
                    type.measure = "mse",  block1.penalization = FALSE, 
                    blocks = blocks, max.coef = c(Inf,Inf),standardize = FALSE,
                    nfolds = 10,lambda.type = "lambda.1se",
                    cvoffset = FALSE)
pl1$block1unpen

pl1$nzero
#the proteins results
pl1$coefficients
coeff1 <- pl1$coefficients
coeff1 <- coeff1[coeff1 != 0]
print(round(coeff1, 4))
names(coeff1)

selected_proteins_list <- list()
# repeat 1000 times
for (i in 1:1000) {
  pl <- prioritylasso(X = pl_pred, Y = pl_out, family = "gaussian",
                      type.measure = "mse",  block1.penalization = FALSE, 
                      blocks = blocks, max.coef = c(Inf,Inf),standardize = FALSE,
                      nfolds = 10, lambda.type = "lambda.1se")
  
  #extract proteins in each model with non-zero cofficident
  selected_proteins <- names(pl$coefficients[pl$coefficients != 0])
  
  # add the extracted proteins name in the list
  selected_proteins_list[[i]] <- selected_proteins
}

# extracted the proteins that frequency >20%
frequent_pro <- table(unlist(selected_proteins_list))
a <- data.frame(frequent_pro)
write.csv(a, file = "supplement/lasso-frequency.csv")
selected_pro <- names(frequent_pro[frequent_pro >= 200])

select_pro_agesex<- data[ ,c("BMI", "age", "sex", "A2M", "ADIPOQ", "AFM","APCS",
                             "APOD", "APOF", "APOM","AZGP1", "CFH",
                             "CRP", "CFI", "GPX3", "LGALS3BP","SERPINA6", "SERPINF1")]
write.csv(select_pro_agesex, file = "select_pro_agesex.csv")

####################machine learning#####################
####################RF###########################
data <- read.csv(file = "select_pro_agesex.csv", row.names = 1)
library(randomForest)

ind<-sample(2,nrow(data),replace = TRUE,prob = c(0.7,0.3))
traindata<-data[ind == 1,]
testdata<-data[ind == 2,]
#use 10 fold cross validation to choose the best parameters
library(caret)
# include ntree, mtry(waiting to be choosen), run two times
grid <- expand.grid(mtry = c(2, 3, 4,5,6,7,8,9, 10))
grid <- expand.grid(mtry = c(100, 200, 300, 400 ,500))
# select by train()
names(getModelInfo())
set.seed(1009)
model <- train(BMI ~ . - age - sex, #model design
               data = traindata,method = "rf",# method RF
               trControl = trainControl(method = "cv",  # method
                                        number = 10) , 
               tuneGrid = grid) 
model[["results"]][["Rsquared"]]
print(model$bestTune)

# Initialize matrices to store raw results
importance_type1_matrices <- list()
importance_type2_matrices <- list()
rmse_values <- numeric(100)
mae_values <- numeric(100)
r_squared_values <- numeric(100)

# Repeat the process 100 times
for (i in 1:100) {
  # Split data into train and test sets
  set.seed(i)  # Set seed for reproducibility
  train_indices <- sample(1:nrow(traindata), size = 0.7 * nrow(traindata))
  traindata_i <- traindata[train_indices, ]
  testdata_i <- traindata[-train_indices, ]
  
  # Train random forest model
  randomForest_rf <- randomForest(BMI ~ . - age - sex, data = traindata_i,
                                  ntree = 400, mtry = 4, trControl = ctrl,
                                  proximity = TRUE, importance = TRUE)
  
  # Predictions
  irispred_i <- predict(randomForest_rf, newdata = testdata_i)
  
  # Variable importance
  importance_type1_i <- randomForest_rf$importance[, "IncNodePurity"]
  importance_type2_i <- randomForest_rf$importance[, "%IncMSE"]
  
  # Store raw importance matrices
  importance_type1_matrices[[i]] <- importance_type1_i
  importance_type2_matrices[[i]] <- importance_type2_i
  
  # Calculate performance metrics
  rmse_i <- sqrt(mean((irispred_i - testdata_i$BMI)^2))
  mae_i <- mean(abs(irispred_i - testdata_i$BMI))
  r_squared_i <- cor(irispred_i, testdata_i$BMI)^2 
  # Store performance metrics
  rmse_values[i] <- rmse_i
  mae_values[i] <- mae_i
  r_squared_values[i] <- r_squared_i}

performance_df <- data.frame(RMSE = rmse_values, MAE = mae_values, R_squared = r_squared_values)
# Create a data frame to store the importance values
importance_df_type1 <- data.frame(matrix(NA, nrow = 100, ncol = ncol(traindata)-1 - 2))
importance_df_type2 <- data.frame(matrix(NA, nrow = 100, ncol = ncol(traindata)-1 - 2))

# Populate the data frame with importance values
for (i in 1:100) {
  importance_df_type1[i, ] <- importance_type1_matrices[[i]]
  importance_df_type2[i, ] <- importance_type2_matrices[[i]]
}

# Set column names(protein names)
colnames(importance_df_type1) <- colnames(data)[4:18]
colnames(importance_df_type2) <- colnames(data)[4:18]
# Calculate column-wise mean
mean_importance_type1 <- colMeans(importance_df_type1)
mean_importance_type2 <- colMeans(importance_df_type2)
# Sort column names based on mean values
sorted_columns <- names(mean_importance_type1)[order(-mean_importance_type1)]
importance_df_type1 <- importance_df_type1[, sorted_columns]
sorted_columns2 <- names(mean_importance_type2)[order(-mean_importance_type2)]
importance_df_type2 <- importance_df_type2[, sorted_columns2]
colMeans(performance_df)

write.csv(importance_df_type1, file = "accuracy_Rf_type1.csv", row.names = FALSE)
write.csv(importance_df_type2, file = "importance_Rf_type2.csv", row.names = FALSE)
write.csv(performance_df, file = "performance_Rf.csv", row.names = FALSE)

library(ggplot2)
# prepare data to match ggplot2
df <- reshape2::melt(importance_df_type1)
#type 1:Accuracy type2:Importance
#Accuracy Importance
p1 <- ggplot(df, aes(x = variable, y = value)) +
  geom_boxplot(width = 0.7, fill = "blue", color = "black", alpha = 0.6) +
  geom_jitter(width = 0.14, alpha = 0.2, color = "red", size = 0.9) +
  labs(x = "", y = "Accuracy", title = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        axis.title.x = element_text(size = 18),  # 调整 X 轴标签字体大小
        axis.title.y = element_text(size = 18),
        plot.title = element_text(hjust = 0.5))
p1
# prepare data to match ggplot2
df <- reshape2::melt(importance_df_type2)
#type 1:Accuracy type2:Importance
#Accuracy Importance
p2 <- ggplot(df, aes(x = variable, y = value)) +
  geom_boxplot(width = 0.7, fill = "blue", color = "black", alpha = 0.6) +
  geom_jitter(width = 0.14, alpha = 0.2, color = "red", size = 0.9) +
  labs(x = "", y = "Importance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        axis.title.x = element_text(size = 18),  # 调整 X 轴标签字体大小
        axis.title.y = element_text(size = 18),
        plot.title = element_text(hjust = 0.5))
p2
rf <- ggarrange(p1,p2, labels = c("A", "B"), font.label = list(size = 20),
                ncol = 1, nrow = 2, common.legend = FALSE,widths = c(4, 4, 4),
                heights = c(4, 4, 4), align = "hv", hjust = 0)
rf
ggsave("plot/rf.png", rf, dpi=600, width= 45, height= 30, unit="cm")
#######################SVM###############################
data <- read.csv(file = "select_pro_agesex.csv", row.names = 1)
library(e1071)
ind<-sample(2,nrow(data),replace = TRUE,prob = c(0.7,0.3))
traindata<-data[ind == 1,]
testdata<-data[ind == 2,]
#chose parameters
set.seed(1009)
model <- train(BMI ~ . - age - sex, #model design
               data = traindata, method = "svmRadial",# method svm
               trControl = trainControl(method = "cv",  # method
                                        number = 10)) 
model[["results"]][["Rsquared"]]
print(model$bestTune)
plot(model)
#svm：nu-regression(quick)eps-regression(high quality)
r_squared_values <- c()
for (i in 1:100) {
  svm <- svm(BMI ~ . - age - sex, 
             data = traindata,
             type = 'eps-regression', kernel = 'radial',
             C = 0.5, sigma = 0.05)
  # prediction
  pre_svm <- predict(svm, newdata = testdata)
  # R squared
  r_squared <- cor(pre_svm, testdata$BMI)^2
  r_squared_values <- c(r_squared_values, r_squared)
}
r_squared <- data.frame(r_squared_values)
colMeans(r_squared)
#####################Neural Network##################
data <- read.csv(file = "select_pro_agesex.csv", row.names = 1)
library(brnn)
library(neuralnet)
ind<-sample(2,nrow(data),replace = TRUE,prob = c(0.7,0.3))
traindata<-data[ind == 1,]
testdata<-data[ind == 2,]
#chose parameters
set.seed(1009)
model <- train(BMI ~ . - age - sex, #model design
               data = traindata, method = "neuralnet",
               trControl = trainControl(method = "cv",  # method
                                        number = 10),
               rep = 1) 
model[["results"]][["Rsquared"]]
print(model$bestTune)

r_squared_values <- c()
for (i in 1:100) {
  network <- neuralnet(BMI ~ . - age - sex,
                       hidden = 1, data = traindata)
  # predict
  pre_net <- predict(network, newdata = testdata)
  # R squared
  r_squared <- cor(pre_net, testdata$BMI)^2
  r_squared_values <- c(r_squared_values, r_squared)
}
colMeans(data.frame(r_squared_values))

#########################################################