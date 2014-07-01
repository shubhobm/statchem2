rm(list=ls())
setwd("C:/Users/Subho/Documents/GitHub/statchem2")
source('supSVD_functions.R')

library(MASS)
library(leaps)

## Load data
load('lta98.rda')
attach(lta98)

## supSVD estimation
modTS = supSVD(X=scale(ltaTS), Y=as.matrix(Y), r=5, quiet=F)

## classfication perfformance??
supSVD.analyze(X=scale(ltaQC),Y=as.matrix(Y<1))

## all subsets regrssion
ltaTS1 = cbind(ltaTS[,1:70],Y)
bestsets = regsubsets(Y~., data=ltaTS1, nbest=5, really.big=T)