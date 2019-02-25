# LARS-EN
Calculate the regression coefficients for the Elastic-Net family

The algorithm will compute the regression coefficients for the whole path for the penalization parameters 
in the Elastic-Net regression considering a variance-covariance matrix **C=X'X**/n for the (centered) predictors and a covariance (between the predictors and the response) vector **cov=X'y**/n.

The code is a modification to the 'lars' function from the R-package LARS by Hastie and Efron.

The following is an example which includes a comparation with the 'glmnet' package.
(so far it implements only the LARS-LASSO models)

```
rm(list=ls())

library(lars)
library(glmnet)

# Load lars2 function
source(url("https://raw.githubusercontent.com/MarcooLopez/LARS-EN/master/lars2.R"))

# Load data
data(diabetes)
X <- diabetes$x2
y <- diabetes$y
n <- nrow(X); p <- ncol(X)

# Standardizing variables.
# Predictors are standardized to have sum(x_i^2)/n equals to one.
y <- scale(y,center=T,scale=T)
X <- scale(X,center=T,scale=F)
normx <- apply(X,2,function(x)sum(x^2))
X <- scale(X,center=F,scale=sqrt(normx/(n)))

# Calculating the variance and covariance matrices
rhs <- drop(crossprod(X,y))/n
C <- crossprod(X)/n

fm1 <- lars2(C,rhs,type="lasso")
fm2 <- lars(X,y,type="lasso")
fm3 <- glmnet(X,y,lambda=fm1$lambda)

# Plot lambda vs coeffcients
out <- data.frame(lambda=c(fm1$lambda,0),fm1$beta)
tt <- melt(out,id=c("lambda"))
ggplot(tt,aes(x=-log(lambda),y=value,group=variable,color=variable))+
geom_line()+theme_bw()+theme(legend.position="none")+labs(y="coefficients")

#==============================================================
# Comparison of solutions for beta (LARS2 vs LARS)
for(k in 1:min(nrow(fm1$beta),nrow(fm2$beta))){
 tmp <- round(sum(abs(fm1$beta[k,]-fm2$beta[k,])),8)
 plot(fm1$beta[k,],fm2$beta[k,],main=paste0("b[",k,"]. D=",tmp),cex=0.7,xlab="LARS2",ylab="LARS")
 abline(a=0,b=1,col=2,lty=2); Sys.sleep(0.001)
}
head(cbind(fm1$beta[88,],fm2$beta[88,]))
k=2
(tt=cbind(fm1$beta[k,],fm2$beta[k,]))

#==============================================================
# Comparison of solutions for beta (LARS2 vs GLMNET)
for(k in 1:min(nrow(fm1$beta),ncol(fm3$beta))){
 tmp <- round(sum(abs(fm1$beta[k,]-fm3$beta[,k])),8)
 plot(fm1$beta[k,],fm3$beta[,k],main=paste0("b[",k,"]. D=",tmp),cex=0.7,xlab="LARS2",ylab="GLMNET")
 abline(a=0,b=1,col=2,lty=2); Sys.sleep(0.01)
}
head(cbind(fm2$beta[88,],fm3$beta[,88]))
k=4
(tt=cbind(fm2$beta[k,],fm3$beta[,k]))

#==============================================================
# Comparison of solutions for beta (LARS vs GLMNET)
for(k in 1:min(nrow(fm1$beta),ncol(fm3$beta))){
 tmp <- round(sum(abs(fm1$beta[k,]-fm3$beta[,k])),8)
 plot(fm1$beta[k,],fm3$beta[,k],main=paste0("b[",k,"]. D=",tmp),cex=0.7,xlab="myLARS",ylab="glmnet")
 abline(a=0,b=1,col=2,lty=2); Sys.sleep(0.01)
}
head(cbind(fm1$beta[50,],fm3$beta[,50]))
k=5
(tt=cbind(fm1$beta[k,],fm3$beta[,k]))

#==============================================================
# Comparison of ALL (LARS2, LARS and GLMNET)
head(cbind(LARS2=fm1$beta[50,],LARS=fm2$beta[50,],GLMNET=fm3$beta[,50]))
k=2
cbind(LARS2=fm1$beta[k,],LARS=fm2$beta[k,],GLMNET=fm3$beta[,k])

```
