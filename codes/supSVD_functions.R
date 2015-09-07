library(psych)
L2norm = function(x) return(sqrt(sum(x^2)))

## function obtain log likelihood given parameter values and data
loglik = function(n,p,X,Y,B,V,Sf,s2){
  
  f1 = V%*%Sf%*%t(V) + s2*diag(p)
  f2 = X - Y%*%B%*%t(V)
  return(-n*p*log(2*pi)/2 - n*log(det(f1))/2
         - tr(f2%*%f1%*%t(f2))/2)
}

## function for log transformation of data matrix
logt = function(X){
  n = nrow(X); p = ncol(X)
  trX = matrix(rep(0,n*p),ncol=p)
  for(i in 1:n){
    for(j in 1:p){
      if(X[i,j]<0) C = -floor(X[i,j])
      else C = 1
      trX[i,j] = log(X[i,j]+C)
    }
  }
  return(trX)
}

## function for implementation of the supSVD algorithm
supSVD = function(X, Y, r, tol=1e-5, maxit=1e3, quiet=T){
  n = nrow(X); p = ncol(X); q = ncol(Y)
  
  # initializing values
  B = diag(nrow=q, ncol=r)
  V = diag(nrow=p, ncol=r)
  Sig.f = diag(1, r); sig2e = 3
  iterating = TRUE; iter=1
  
  # creating variables for later use
  XtX = crossprod(X)
  
  while(iterating){
    # creating variables for later use
    dSig.f = diag(Sig.f)
    sig2.by.dSigf = sig2e/dSig.f
    
    # E step
    EU.cond = (Y%*%B%*%diag(sig2.by.dSigf, r) + X%*%V) %*% diag((1 + sig2.by.dSigf)^(-1), r)
    VarU.cond = diag((1/dSig.f+1/sig2e)^(-1), r)
    
    # M step
    Bh = as.matrix(solve(crossprod(Y))%*%t(Y)%*%EU.cond)
    
    EUtEU = crossprod(EU.cond)
    YBh = Y%*%Bh
    BtYtEU = t(YBh) %*% EU.cond
    
    Vh = (t(X)%*%EU.cond) %*% solve(n*VarU.cond + EUtEU)
    VEUt = V%*%t(EU.cond)
    
    Sig.fh = (n*VarU.cond + EUtEU + crossprod(YBh) - BtYtEU - t(BtYtEU))/n
    sig2eh = (tr(XtX) - 2*tr(t(VEUt)%*%t(X)) + n*tr(crossprod(Vh)%*%VarU.cond) + tr(crossprod(VEUt)))/(n*p)
    
    # S step
    e = eigen(Vh%*%Sig.fh%*%t(Vh))
    V1 = e$vectors[,1:r]
    dSig.f1 = e$values[1:r]
    B1 = Bh%*%t(Vh)%*%V1
    sig2e1 = sig2eh
    
    # reorder columns
    order.norm = order(apply(X%*%V1, 2, L2norm), decreasing=T)
    V11 = V1; B11 = B1; dSig.f11 = dSig.f1
    for(i in 1:r){
      whi = which(order.norm==i)
      V11[,i] = V1[,whi]
      B11[,i] = B1[,whi]
      dSig.f11[i] = dSig.f1[whi]
    }
    Sig.f11 = diag(dSig.f11)
    
    # stopping criterion
    l1 = loglik(n,p,X,Y,B11,V11,Sig.f11,sig2e1)
    l0 = loglik(n,p,X,Y,B,V,Sig.f,sig2e)
    if(!quiet)
      cat(iter, l1, l0, dSig.f11, "\n")
    err = l1 - l0
    if(abs(err)<tol || iter==maxit)
      iterating = FALSE
    else
      B = B11; V = V11; Sig.f = Sig.f11; sig2e = sig2e1; iter = iter+1
  }
  
  if(abs(err)>tol && iter==maxit)
    warning("Minimum tolerance limit reached before convergence")
  return(list(B=B1, V=V1, Sigma.f=diag(Sig.f11), sigma2 = sig2e1))  
}

## function for prediction performance of supSVD
supSVD.analyze = function(X, Y, tol=1e-5, maxit=1e3, quiet=T){
  
  # PCA of data matrix
  s2 = svd(X)$d; ss2 = sum(s2)
  r = min(which(cumsum(s2)>=0.90*ss2))
  s2.sel = s2[1:r]
  
    # fit supSVD models
  mod = supSVD(X=X, Y=as.matrix(Y), r=r, tol=tol, maxit=maxit, quiet=quiet)
  # project data onto loading vectors
  pX = X%*%(mod$V)
  # pairs(pX, col=ifelse(Y==1,'red','blue'),pch=19)
  # lda using first 7 PCA transformed predictors and sSVD transformed predictors
  sSVD.lda = lda(Y~pX, CV=T)
  tab = table(Y, sSVD.lda$class)
  tab.prop = rbind(tab[1, ]/sum(tab[1, ]), tab[2, ]/sum(tab[2, ]))
  
  return(list(model=mod, dim=r, VBeforeAfter = matrix(c(s2.sel,mod$Sigma.f),ncol=2, byrow=FALSE), TruthTable=tab.prop))
}