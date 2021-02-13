library(ggplot2)

# read excel file
# data <- read_excel("path/KTSP_Results.xlsx")
metabolites <- data$METABOLITES
acc <- data$ACCURACY
mcc <- data$MCC

#1
ggplot(data, aes(x=acc, y=mcc)) + geom_point()+ geom_smooth(method=lm, se=FALSE) + labs(title = "Accuracy vs MCC for Threshold -1", y = "MCC", x = "Accuracy")

#2 - Scatterplots based on top metabolite pairs
# Select the metabolite whose top pairs are to be plotted
target <- HM_Data$Target...48
target <- as.factor(target)
# Select histone 1 from the pair
H3K4me1 <- HM_Data$H3K4me1
# Select histone 2 from the pair 
H3K27me0K36me0 <- HM_Data$H3K27me0K36me0

p <- ggplot(HM_Data, aes(x=H3K4me1, y=H3K27me0K36me0, shape=target, color=target)) + geom_point() + geom_rug() 
p