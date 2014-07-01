rm(list=ls())
setwd("C:/Users/Subho/Documents/GitHub/statchem2")
lta98 = read.csv("lta98.csv", header=T)

## separating variable categories
cats = factor(as.numeric(lta98[1,]))
lta = lta98[,which(cats!=0)]
Y = lta98[,which(colnames(lta98)=='ACT')]

## separate different variable types
ltaTS = lta98[,which(cats==1)]
SMILES = ltaTS$SMILES
ltaTS = ltaTS[,-which(colnames(ltaTS)=="SMILES")]
ltaTC = lta98[,which(cats==2)]
lta3D = lta98[,which(cats==3)]
ltaQC = lta98[,which(cats==4)]

lta98 = list(Y=Y,SMILES=SMILES, ltaTS=ltaTS,ltaTC=ltaTC,lta3D=lta3D,ltaQC=ltaQC)
save(lta98, file='lta98.rda')

# load('lta98.rda')
