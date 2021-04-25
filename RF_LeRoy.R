library(randomForest)
library(pROC)
library(caTools)
library(pracma)

final_acc <- c()
final_mcc <- c()
final_prec <- c()
final_rec <- c()
final_f1 <- c()
final_guess_acc <- c()

for (metab in 1:136)
{
  set.seed(123)
  DataFrame <- HM_Data[138:188]
  DataFrame$target <- HM_Data[metab]
  DataFrame <- as.data.frame(DataFrame)
  df <- data.frame(DataFrame$target)
  DataFrame$target = as.factor(df$Target)
  
  DataFrame_test <- LeRoy[138:188]
  DataFrame_test$target <- LeRoy[metab]
  DataFrame_test <- as.data.frame(DataFrame_test)
  df_test <- data.frame(DataFrame_test$target)
  DataFrame_test$target = as.factor(df_test$Target)
  
  trainDF <- DataFrame
  testDF <- DataFrame_test
  
  # Best Mtry
  trainDF2 = trainDF[,1:51]
  best_mtry <- tuneRF(x = trainDF2, 
                      y = trainDF$target, 
                      mtryStart = 17, 
                      ntreeTry = 500, 
                      stepFactor = 1.2, 
                      improve = 0.0001, 
                      trace = T, 
                      plot = T, 
                      doBest = T,
                      importance = T)
  
  # Model
  modelRandom <- randomForest(target~., data = trainDF, mtry = best_mtry$mtry, ntree = 500)
  #modelRandom
  #my_list[k] <- importance(modelRandom)
  
  # Prediction
  predWithClass <- predict(modelRandom, testDF, type = 'class')
  t <- table(predictions = predWithClass, actual = testDF$target)
  
  # MCC Calculation
  TP <- t[1]
  TN <- t[4]
  FP <- t[3]
  FN <- t[2]
  N <- TN + TP + FN + FP
  S <- (TP + FN)/N
  P <- (TP + FP)/N
  
  # Performance Parameters
  final_acc[metab] <- sum(diag(t))/sum(t)
  final_mcc[metab] <- ((TP/N) - (S*P))/(sqrt(P*S*(1-S)*(1-P)))
  final_prec[metab] <- TP/(TP+FP)
  final_rec[metab] <- TP/(TP+FN)
  final_f1[metab] <- 2*final_rec[metab]*final_prec[metab]/(final_rec[metab]+final_prec[metab])
  
  guess_acc <- rep(0,100)
  for (j in 1:100)
  {
    Y_fac <- DataFrame$target
    guess <- Y_fac[randperm(length(Y_fac))]
    guess_acc[j] =sum(Y_fac == guess)/length(Y_fac)
  }
  
  final_guess_acc[metab] <- mean(guess_acc, na.rm = TRUE)
  
  print(paste("Finished processing metabolite -", metab))
}
  
  
# write the final performance parameter values for each metabolite in text file for storage.
write.table(final_acc, "LeRoy_RF_acc_-0.5.txt", sep="\t")
write.table(final_mcc, "LeRoy_RF_mcc_-0.5.txt", sep="\t")
write.table(final_prec, "LeRoy_RF_prec_-0.5.txt", sep="\t")
write.table(final_rec, "LeRoy_RF_rec_-0.5.txt", sep="\t")
write.table(final_f1, "LeRoy_RF_f1_-0.5.txt", sep="\t")
write.table(final_guess_acc, "LeRoy_RF_guess_acc_-0.5.txt", sep="\t")
