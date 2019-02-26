# LARS-EN
Calculate the regression coefficients for the Elastic-Net family

The algorithm will compute the regression coefficients for the whole path for the penalization parameters 
in the Elastic-Net regression considering a variance-covariance matrix (upon a constant *k*) **C=X'X**/k for the (centered) predictors and a covariance (between the predictors and the response) vector **cov=X'y**/k.

The code is a modification to the 'lars' function from the R-package LARS by Hastie and Efron.

The following is an example which includes a comparison with the 'glmnet' package.
(so far it implements only the LARS-LASSO models).

Simulated data will be used using a homemade function 'simData'.

```r
rm(list=ls())

library(lars)
library(glmnet)
library(lattice)

# Load lars2 function
source(url("https://raw.githubusercontent.com/MarcooLopez/LARS-EN/master/lars2.R"))

# Load the function to simulate data
source(url("https://raw.githubusercontent.com/MarcooLopez/LARS-EN/master/simu_data.R"))

# Generate data with p=250 predictors and n=1000 observations
DATA <- simData(n=500,p=50,seed=123)
X0 <- DATA$X
y0 <- DATA$y
n <- nrow(X0); p <- ncol(X0)
```

### Using predictors and response having mean zero and norm equal to *n*
```r
y <- scale(y0,center=TRUE,scale=TRUE)
X <- scale(X0,center=TRUE,scale=TRUE)
normx <- apply(X,2,function(x)sum(x^2))
X <- scale(X,center=FALSE,scale=sqrt(normx/n))
y <- scale(y,center=FALSE,scale=sqrt(sum(y^2)/n))

# Calculating the variance and covariance matrices
rhs <- drop(crossprod(X,y))/n
C <- crossprod(X)/n

fm1 <- lars2(C,rhs,type="lasso")
fm2 <- glmnet(X,y,lambda=c(fm1$lambda,0),thresh=1E-10)

# Comparing both estimations
B1 <- fm1$beta
B2 <- t(as.matrix(fm2$beta))
round(cbind(B1[20,],B2[20,]),6)
D <- abs(B1-B2)
rownames(D) <- NULL
levelplot(t(D),xlab="predictor j",ylab="Step k",
 main=expression('|'*beta[j]^lars*'('*lambda[k]*')-'*beta[j]^glmnet*'('*lambda[k]*')|'))
```

### Using predictors and response having mean zero and norm equal to *n-1*
```r
y <- scale(y0,center=TRUE,scale=TRUE)
X <- scale(X0,center=TRUE,scale=TRUE)

# Calculating the variance and covariance matrices
rhs <- drop(crossprod(X,y))/(n-1)
C <- crossprod(X)/(n-1)

fm1 <- lars2(C,rhs,type="lasso")
fm2 <- glmnet(X,y,lambda=sqrt((n-1)/n)*c(fm1$lambda,0),thresh=1E-10)

# Comparing both estimations
B1 <- fm1$beta
B2 <- t(as.matrix(fm2$beta))
round(cbind(B1[20,],B2[20,]),6)
D <- abs(B1-B2)
rownames(D) <- NULL
levelplot(t(D),xlab="predictor j",ylab="Step k",
 main=expression('|'*beta[j]^lars*'('*lambda[k]*')-'*beta[j]^glmnet*'('*lambda[k]*')|'))
```
##
Levelplot indicating the absolute pointwise difference of the 'betas' obtained by lars2 and glmnet at each step of the LARS-LASSO.

<p align="center">
<img src="https://github.com/MarcooLopez/LARS-EN/blob/master/levelplot_abs_diff.png" width="450">
</b>
