setwd("C:/Users/Subho/Documents/GitHub/statchem2")
lta98 = read.csv("lta98.csv", header=T)

## separating variable categories
cats = factor(as.numeric(lta98[1,]))
lta = lta98[,which(cats!=0)]
Y = lta98[,which(colnames(lta98)=='ACT')]