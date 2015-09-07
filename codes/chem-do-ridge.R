rm(list=ls())
setwd("C:/Users/Subho/Documents/GitHub/statchem2")
library(MASS)
library(parcor)
library(ridge)
###### required for parallel computing
library(parallel)
library(doSNOW)

## function for two-deep CV
cv.deep = function(formula, data, cores=1){
  n = nrow(data)
  parfun = function(i){
    require(ridge)
    imod = logisticRidge(formula, lambda="automatic", data=data, subset=-i)
    predict(imod, data[i,], type="response")
  }
  if(cores==1){
    preds = lapply(1:n, parfun)
  }
  else{
    cl = makeCluster(cores)
    registerDoSNOW(cl)
    preds = foreach(i=1:n) %dopar% parfun(i)
    stopCluster(cl)
  }
  as.numeric(preds)
}

## Load data
load('lta98.rda')
attach(lta98)
ltaTS = ltaTS[-1,]
ltaTC = ltaTC[-1,]
ltaQC = ltaQC[-1,]
lta3D = lta3D[-1,]
Y = Y[-1]
Y01 = (Y>0)

## histogram of Y
# hist(Y, breaks=15,
#      col="grey", border=gray(.5),
#      main="",xlab="log R",
#      xlim=c(-4,4))
# abline(v=0, lwd=4, col="red")

system.time(pred1 <- cv.deep(Y01~., cbind(Y01,ltaTS), cores=3))
table(Y01, pred1>mean(Y01))

system.time(pred2 <- cv.deep(Y01~., cbind(Y01,ltaTS,ltaTC), cores=3))
table(Y01, pred2>mean(Y01))

system.time(pred3 <- cv.deep(Y01~., cbind(Y01,ltaTS,ltaTC,lta3D), cores=3))
table(Y01, pred3>mean(Y01))

system.time(pred4 <- cv.deep(Y01~., cbind(Y01,ltaTS,ltaTC,lta3D,ltaQC), cores=3))
table(Y01, pred4>mean(Y01))

## 508 data
data508 = read.table("data508.txt", header=T)
Y01 = data508[-1,309]
cats = data508[1,]
mut508 = list(Y01 = Y01,
              mutTS = data508[-1,which(cats==1)],
              mutTC = data508[-1,which(cats==2)],
              mut3D = data508[-1,which(cats==3)],
              mutQC = data508[-1,which(cats==4)])

attach(mut508)
system.time(pred508.1 <- cv.deep(Y01~., cbind(Y01,mutTS), cores=3))
table(Y01, pred508.1>mean(Y01))

system.time(pred508.2 <- cv.deep(Y01~., cbind(Y01,mutTS,mutTC), cores=3))
table(Y01, pred508.2>mean(Y01))

system.time(pred508.1 <- cv.deep(Y01~., cbind(Y01,mutTS,mutTC,mut3D), cores=3))
table(Y01, pred508.1>mean(Y01))

system.time(pred508.1 <- cv.deep(Y01~., cbind(Y01,mutTS,mutTC,mut3D,mutQC), cores=3))
table(Y01, pred508.1>mean(Y01))
