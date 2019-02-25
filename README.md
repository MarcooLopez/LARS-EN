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
library(reshape)
library(ggplot2)

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
X <- scale(X,center=F,scale=sqrt(normx/n))

# Calculating the variance and covariance matrices
rhs <- drop(crossprod(X,y))/n
C <- crossprod(X)/n

fm1 <- lars2(C,rhs,type="lasso")
fm2 <- lars(X,y,type="lasso")
fm3 <- glmnet(X,y,lambda=fm1$lambda,thresh=1E-12)

#==============================================================
D <- c()
for(k in 1:min(nrow(fm1$beta),nrow(fm2$beta),ncol(fm3$beta))){
 tmp <- c(k,sum(abs(fm1$beta[k,]-fm2$beta[k,])),
        sum(abs(fm1$beta[k,]-fm3$beta[,k])),
        sum(abs(fm2$beta[k,]-fm3$beta[,k])))
 D <- rbind(D,round(tmp,6))
}
colnames(D) <- c("betaj","lars2 vs lars","lars2 vs glmnet","lars vs glmnet")
D <- melt(data.frame(D),id="betaj")
ggplot(D,aes(x=betaj,y=value))+geom_point()+facet_wrap(~variable)

# Print an specific subset
k=90
cbind(LARS2=fm1$beta[k,],LARS=fm2$beta[k,],GLMNET=fm3$beta[,k])

```
