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
