library(dplyr)
library(tidyr)
library(caret)
library(pracma)
require(switchBox)

# read excel file
# HM_Data <- read_excel("path/hist_metabolite_Z.xlsx")

# Indices of training data 
# Histone Markers which are used as training data set are present in excel sheet from columns 138:188
X <- HM_Data %>% select(138:188)
indices <- sample(10,nrow(HM_Data), replace = TRUE, prob = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))
final_acc <- c()
final_mcc <- c()
final_prec <- c()
final_rec <- c()
final_f1 <- c()
final_guess_acc <- c()
p_val <- c()

# The classification is repeated for 136 metabolites which are present in columns 1:136 in excel file.
for (metab in 1:136)
{
  Y <- HM_Data[metab]
  acc <- rep(10,1)
  prec <- rep(10,1)
  rec <- rep(10,1)
  F1 <- rep(10,1)
  MCC <- rep(10,1)
  
  for (k in 1:10)  #10-fold Cross Validation
  {
    test <- indices == k 
    train <- indices != k
    Xtrain <- X[train,]
    Ytrain <- Y[train,]
    Xtest <- X[test,]
    Ytest <- Y[test,]
    
    # Taking transpose to bring data to correct format
    Xtrain <- t(Xtrain)
    Xtest <- t(Xtest)
    
    aa <- data.frame(Ytrain)
    aa$Target <- as.factor(aa$Target)
    Ytrain_fac <- aa$Target
    
    aa_test <- data.frame(Ytest)
    aa_test$Target <- as.factor(aa_test$Target)
    Ytest_fac <- aa_test$Target
    
    # Train a classifier using default filtering function based on the Wilcoxon test
    classifier <- SWAP.Train.KTSP(Xtrain, Ytrain_fac, krange=c(3:15))
    
    testPrediction <- SWAP.KTSP.Classify(Xtest, classifier)
    
    # Resubstitution performance in the TEST set
    tab <- table(testPrediction, Ytest_fac)
    acc[k] <- sum(Ytest_fac == testPrediction)/length(Ytest_fac)
    prec[k] <- sum(Ytest_fac == 1 & testPrediction == 1)/sum(testPrediction == 1)
    rec[k] <- sum(Ytest_fac == 1 & testPrediction == 1)/sum(Ytest_fac == 1)
    F1[k] <- 2*rec[k]*prec[k]/(rec[k]+prec[k])
    
    # MCC Calculation
    TP <- tab[1]
    TN <- tab[4]
    FP <- tab[3]
    FN <- tab[2]
    N <- TN + TP + FN + FP
    S <- (TP + FN)/N
    P <- (TP + FP)/N
    MCC[k] <- ((TP/N) - (S*P))/(sqrt(P*S*(1-S)*(1-P)))
  }
  
  # Random guess accuracy
  guess_acc <- rep(0,100)
  for (j in 1:100)
  {
    guess <- Ytest_fac[randperm(length(Ytest_fac))]
    guess_acc[j] =sum(Ytest_fac == guess)/length(Ytest_fac)
  }
  
  # Final performance parameter values for each metabolite
  final_acc[metab] <- mean(acc)
  final_mcc[metab] <- mean(MCC)
  final_prec[metab] <- mean(prec)
  final_rec[metab] <- mean(rec)
  final_f1[metab] <- mean(F1)
  final_guess_acc[metab] <- mean(guess_acc)
  p_val[metab] <- t.test(acc, guess_acc)$p.value
  
}

# write the final performance parameter values for each metabolite in text file for storage.
write.table(final_acc, "path/acc_-1.txt", sep="\t")
write.table(final_mcc, "path/mcc_-1.txt", sep="\t")
write.table(final_prec, "path/prec_-1.txt", sep="\t")
write.table(final_rec, "path/rec_-1.txt", sep="\t")
write.table(final_f1, "path/f1_-1.txt", sep="\t")
write.table(final_guess_acc, "path/guess_acc_-1.txt", sep="\t")
write.table(p_val, "path/pval_-1.txt", sep="\t")






