draw_confusion_matrix <- function(cm) {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)
  
  # create the matrix 
  rect(150, 430, 240, 370, col='#E89DA0')
  text(195, 435, 'Non-Obese', cex=1.2)
  rect(250, 430, 340, 370, col='#88CEE6')
  text(295, 435, 'Obese', cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='#88CEE6')
  rect(250, 305, 340, 365, col='#E89DA0')
  text(140, 400, 'Non-Obese', cex=1.2, srt=90)
  text(140, 335, 'Obese', cex=1.2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  
draw_confusion_matrix(confusionMatrix(testdata$BMI, irispred)) 


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

write.csv(result, file = "result/mc_pro_rresult.csv")

write.csv(result, file = "result/mc_pheno_rresult.csv")
confusionMatrix(testdata$BMI, ad)$overall['AccuracyNull']
