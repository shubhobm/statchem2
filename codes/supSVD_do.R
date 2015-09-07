library(MASS)
library(stats)
library(pracma)
library(psych)
setwd("C:/Study/UMN files/8932")
source("supSVD_functions.R")

## Apply on simulated data
set.seed(42914)
n = 1000; p = 10; q = 2; r = 3
Y = matrix(rbinom(2*n, size=1, prob=.5), ncol=2)
B = mvrnorm(n=q, mu=rep(0,r), Sigma=diag(r))
V = gramSchmidt(mvrnorm(n=p, mu=rep(0,r), Sigma=diag(r)))$Q
Sigma.f = diag(c(10,7,4)); sigma2 = 2
F = mvrnorm(n=n, mu=rep(0,r), Sigma=Sigma.f)
E = mvrnorm(n=n, mu=rep(0,p), Sigma=sigma2*diag(p))
X = Y%*%B%*%t(V) + F%*%t(V) + E

a=supSVD(X=X, Y=Y, r=3, quiet=T)

## estimation performance on simulated data
pr = p*r
err.list = matrix(rep(0,4*100),ncol=4)
for(i in 1:100){
  # supSVD case
  Y = matrix(rbinom(2*n, size=1, prob=.5), ncol=2)
  B = mvrnorm(n=q, mu=rep(0,r), Sigma=diag(r))
  V = gramSchmidt(mvrnorm(n=p, mu=rep(0,r), Sigma=diag(r)))$Q
  Sigma.f = diag(c(10,7,4)); sigma2 = 2
  F = mvrnorm(n=n, mu=rep(0,r), Sigma=Sigma.f)
  E = mvrnorm(n=n, mu=rep(0,p), Sigma=sigma2*diag(p))
  X = Y%*%B%*%t(V) + F%*%t(V) + E
  
  a=supSVD(X=X, Y=Y, r=3, quiet=T)
  pcmod=princomp(X)
  err.list[i,1] = norm(V-a$V, type="F")/pr
  err.list[i,3] = norm(V-pcmod$loadings[,1:r], type="F")/pr
  
  # PCA case
  X = F%*%t(V) + E
  a=supSVD(X=X, Y=Y, r=3, quiet=T)
  pcmod=princomp(X)
  err.list[i,2] = norm(V-a$V, type="F")/pr
  err.list[i,4] = norm(V-pcmod$loadings[,1:r], type="F")/pr
}
apply(err.list,2,mean)
apply(err.list,2,sd)

## prediction performance on simulated data
truth.list=matrix(rep(0,4*100),ncol=4)
for(i in 1:100){
  # creating data
  Y = matrix(c(rep(0,75),rep(1,125)),ncol=1)
  X1 = mvrnorm(n=75, mu=c(-40,30), Sigma=diag(c(40,1560)))
  X2 = mvrnorm(n=125, mu=c(-10,-30), Sigma=diag(c(55,35)))
  X = cbind(rbind(X1,X2),mvrnorm(n=200,mu=rep(0,8),Sigma=diag(8)))
  
  # doing only LDA
  m.lda = lda(Y~X, CV=T)
  tab = table(Y, m.lda$class)
  tab.prop = rbind(tab[1, ]/sum(tab[1, ]), tab[2, ]/sum(tab[2, ]))
  truth.list[i,1:2] = c(tab.prop[1,1],tab.prop[2,2])
  
  # doing supSVD and then LDA
  a=supSVD(X=X, Y=Y, r=2, quiet=T)
  pX = X%*%(a$V)
  sSVD.lda = lda(Y~pX, CV=T)
  tab = table(Y, sSVD.lda$class)
  tab.prop = rbind(tab[1, ]/sum(tab[1, ]), tab[2, ]/sum(tab[2, ]))
  truth.list[i,3:4] = c(tab.prop[1,1],tab.prop[2,2])
}
apply(truth.list,2,mean)
apply(truth.list,2,sd)

## Chemometrics data
data508 = read.csv("data508.csv", header=T)
X = logt(as.matrix(data508[-1,-c(1,309)]))
Y = as.matrix(data508[-1,309])
vartype = as.factor(as.matrix(data508[1,-c(1,309)]))

Xts = X[,vartype==1]
Xtc = X[,vartype==2]

## get outputs
# standardized data
supSVD.analyze(X=scale(Xts), Y=scale(Y), quiet=T)
supSVD.analyze(X=scale(Xtc), Y=scale(Y), quiet=T)
supSVD.analyze(X=scale(cbind(Xts,Xtc)), Y=scale(Y), quiet=T)
supSVD.analyze(X=scale(X), Y=scale(Y), quiet=T)

# scatterplot
sXts=scale(Xts)%*%(m$model$V)
pairs(sXts[,1:5], col=ifelse(vartype==1,'red','blue'), pch=19, cex=.8)