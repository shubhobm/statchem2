rm(list=ls())
setwd("C:/Users/Subho/Documents/GitHub/statchem2")
source('supSVD_functions.R')

library(MASS)
library(leaps)

## Load data
load('lta98.rda')
attach(lta98)
ltaTS = ltaTS[-1,]
ltaTC = ltaTC[-1,]
ltaQC = ltaQC[-1,]
lta3D = lta3D[-1,]
Y = Y[-1]

## supSVD estimation
modTS = supSVD(Y=scale(ltaTS), X=as.matrix(Y), r=1, quiet=T)

## classfication perfformance??
supSVD.analyze(X=scale(cbind(ltaTS,ltaTC)),Y=as.matrix(Y>0), quiet=T)

## all subsets regrssion
ltaTS1 = cbind(ltaTS[,1:70],Y)
bestsets = regsubsets(Y~., data=ltaTS1, nbest=5, really.big=T)