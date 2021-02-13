library(dplyr)
library(tidyr)
library(caret)
library(pracma)
require(switchBox)

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

# Index of metabolite whose top pairs are to be calculated.
Y <- HM_Data[116]
acc <- rep(10,1)
prec <- rep(10,1)
rec <- rep(10,1)
F1 <- rep(10,1)
MCC <- rep(10,1)
my_list_of_lists <- list()

for (k in 1:10)  #10-fold Cross Validation
{
  test <- indices == k 
  train <- indices != k
  Xtrain <- X[train,]
  Ytrain <- Y[train,]
  Xtest <- X[test,]
  Ytest <- Y[test,]
  
  # Taling transpose to bring data in right format
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
  my_list_of_lists <- append(my_list_of_lists, list(classifier))
}

# Random guess accuracy
guess_acc <- rep(0,100)
for (j in 1:100)
{
  guess <- Ytest_fac[randperm(length(Ytest_fac))]
  guess_acc[j] =sum(Ytest_fac == guess)/length(Ytest_fac)
}

# The performance parameters of all 10 classifiers are checked and 
# the classifier with optimal performance parameters is selected. 
# The top scoring pairs are calculated for that model.

# Digit in bracket will change according to best classifier  
ind_max_classifier <- 10
write.table(my_list_of_lists[ind_max_classifier][[1]]$score, "path/topPairs.txt", sep="\t")






