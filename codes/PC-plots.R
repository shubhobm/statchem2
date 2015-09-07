setwd("C:/Users/Subho/Documents/GitHub/statchem2")
source("./Codes/supSVD_functions.R")

## load 508 data
data508 = read.table("data508.txt", header=T)
p = ncol(data508)
X = logt(as.matrix(data508[-1,-c(1,p)]))
Y = as.matrix(data508[-1,p])
vartype = as.factor(as.matrix(data508[1,-c(1,p)]))
pch508 = ifelse(Y==0, 1, 8)

Xts = X[,vartype==1]
svd508 = svd(Xts)
tX508 = scale(Xts) %*% svd508$v

## load 95 data
lta98 = read.csv("lta98.csv", header=T)
minus = c(1,194,195,196)
X = logt(as.matrix(lta98[-1,-minus]))
Y = as.matrix(lta98[-1,195])>0
vartype = as.factor(as.matrix(lta98[1,-minus]))
pch98 = ifelse(Y==0, 1, 8)

Xts = X[,vartype==1]
svd98 = svd(Xts)
tX98 = scale(Xts) %*% svd98$v

defPar = par()
pdf('printsGreat.pdf', width = 8.5, height = 11)
par(cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5, mfrow=c(3,2), mar=rep(5,4))
plot(tX508[,1], tX508[,2], pch=pch508, cex=1.5, xlab="PC1", ylab="PC2",
     main="508 mutagen data")
plot(tX98[,1], tX98[,2], pch=pch98, cex=1.5, xlab="PC1", ylab="PC2",
     main="95 amine data")

plot(tX508[,2], tX508[,3], pch=pch508, cex=1.5, xlab="PC2", ylab="PC3")
plot(tX98[,2], tX98[,3], pch=pch98, cex=1.5, xlab="PC2", ylab="PC3")

plot(tX508[,1], tX508[,3], pch=pch508, cex=1.5, xlab="PC1", ylab="PC3")
plot(tX98[,1], tX98[,3], pch=pch98, cex=1.5, xlab="PC1", ylab="PC3")
dev.off()
par(defPar)
