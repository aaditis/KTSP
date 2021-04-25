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
  
  indices <- sample(10,nrow(HM_Data), replace = TRUE, prob = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))
  acc <- rep(10,1)
  MCC <- rep(10,1)
  prec <- rep(10,1)
  rec <- rep(10,1)
  F1 <- rep(10,1)
  
  for (k in 1:10)
  {
    test <- indices == k
    train <- indices != k
    trainDF <- DataFrame[train,]
    testDF <- DataFrame[test,]
    
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
    acc[k] <- sum(diag(t))/sum(t)
    MCC[k] <- ((TP/N) - (S*P))/(sqrt(P*S*(1-S)*(1-P)))
    prec[k] <- TP/(TP+FP)
    rec[k] <- TP/(TP+FN)
    F1[k] <- 2*rec[k]*prec[k]/(rec[k]+prec[k])
  }
  
  guess_acc <- rep(0,100)
  for (j in 1:100)
  {
    Y_fac <- DataFrame$target
    guess <- Y_fac[randperm(length(Y_fac))]
    guess_acc[j] =sum(Y_fac == guess)/length(Y_fac)
  }
  
  # Final performance parameter values for each metabolite
  final_acc[metab] <- mean(acc)
  final_mcc[metab] <- mean(MCC)
  final_prec[metab] <- mean(prec)
  final_rec[metab] <- mean(rec)
  final_f1[metab] <- mean(F1)
  final_guess_acc[metab] <- mean(guess_acc)
}

# write the final performance parameter values for each metabolite in text file for storage.
write.table(final_acc, "C:/Users/Vivek/Desktop/MICHIGAN/SEM4_COURSES/Independent Study/RandomForest/RF_acc_0.txt", sep="\t")
write.table(final_mcc, "C:/Users/Vivek/Desktop/MICHIGAN/SEM4_COURSES/Independent Study/RandomForest/RF_mcc_0.txt", sep="\t")
write.table(final_prec, "C:/Users/Vivek/Desktop/MICHIGAN/SEM4_COURSES/Independent Study/RandomForest/RF_prec_0.txt", sep="\t")
write.table(final_rec, "C:/Users/Vivek/Desktop/MICHIGAN/SEM4_COURSES/Independent Study/RandomForest/RF_rec_0.txt", sep="\t")
write.table(final_f1, "C:/Users/Vivek/Desktop/MICHIGAN/SEM4_COURSES/Independent Study/RandomForest/RF_f1_0.txt", sep="\t")
write.table(final_guess_acc, "C:/Users/Vivek/Desktop/MICHIGAN/SEM4_COURSES/Independent Study/RandomForest/RF_guess_acc_0.txt", sep="\t")

######################################################
# Additional
#predWithProb <- predict(modelRandom, testDF, type = 'prob')
#auc <- auc(testDF$target, predWithProb[,2])
#plot(roc(testDF$target, predWithProb[,2]))

