# KTSP
Algorithm to implement 10-fold cross validation KTSP in R

The metabolite data comprises of 225 metabolites in 928 cell lines from more than 20 cancer types in the Cancer Cell Line Encyclopedia profiled using liquid chromatography-mass spectrometry (LC-MS). The histone dataset from global chromatin profiling consists of 7 distinct histone markers along with varying the methylation and acetylation in the histones bringing the total number of histones to 42 for 897 cell lines. Histone methylation modifies certain amino acids in a histone protein by addition of one, two or three methyl groups.

The goal of the analysis was to use the histone markers to classify CCLE metabolite data. Common cell lines were considered from the 928 cell lines in the metabolite data and 897 cell lines from the histone data to form the final dataset which comprised of 879 cell lines. The 42 histones were divided into 4 categories and sum and variation of histones in each category was calculated increasing the feature space by 8 features. Doubling time was an additional feature which brought the total feature space to 51 features. 


