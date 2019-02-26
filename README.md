# LARS-EN
Calculate the regression coefficients for the Elastic-Net family

The algorithm will compute the regression coefficients for the whole path for the penalization parameters 
in the Elastic-Net regression considering a variance-covariance matrix (upon a constant *k*) **C=X'X**/k for the (centered) predictors and a covariance (between the predictors and the response) vector **cov=X'y**/k.

The code is a modification to the 'lars' function from the R-package LARS by Hastie and Efron.

The following is an example which includes a comparation with the 'glmnet' package.
(so far it implements only the LARS-LASSO models).

Simulated data will be used using a homemade function 'simData'.

```
rm(list=ls())

library(lars)
library(glmnet)
library(reshape)
library(ggplot2)

# Load lars2 function
source(url("https://raw.githubusercontent.com/MarcooLopez/LARS-EN/master/lars2.R"))

# Load the function to simulate data
source(url("https://raw.githubusercontent.com/MarcooLopez/LARS-EN/master/simu_data.R"))

# Generate data with p=200 preductors and n=500 obsertvations
DATA <- simData(n=500,p=200)
X <- DATA$X
y <- DATA$y
n <- nrow(X); p <- ncol(X)
```

## Using predictors and response having mean zero and norm equal to n
```{r}
y <- scale(y,center=TRUE,scale=TRUE)
X <- scale(X,center=TRUE,scale=TRUE)
normx <- apply(X,2,function(x)sum(x^2))
X <- scale(X,center=FALSE,scale=sqrt(normx/n))
y <- scale(y,center=FALSE,scale=sqrt(sum(y^2)/n))

# Calculating the variance and covariance matrices
rhs <- drop(crossprod(X,y))/n
C <- crossprod(X)/n

fm1 <- lars2(C,rhs,type="lasso")
fm2 <- lars(X,y,type="lasso")
fm3 <- glmnet(X,y,lambda=fm1$lambda[fm1$lambda>0.001],thresh=1E-10)
```
#==============================================================
D <- c()
for(k in 1:min(nrow(fm1$beta),nrow(fm2$beta),ncol(fm3$beta))){
 tmp <- c(k,sum((fm1$beta[k,]-fm2$beta[k,])^2),
        sum((fm1$beta[k,]-fm3$beta[,k])^2),
        sum((fm2$beta[k,]-fm3$beta[,k])^2))
 D <- rbind(D,round(tmp,6))
}
colnames(D) <- c("betaj","lars2 vs lars","lars2 vs glmnet","lars vs glmnet")
D <- melt(data.frame(D),id="betaj")
ggplot(D,aes(x=betaj,y=value))+geom_point()+facet_wrap(~variable)+labs(y="sum(x-y)^2")

# Print an specific subset
k=50
cbind(LARS2=fm1$beta[k,],LARS=fm2$beta[k,],GLMNET=fm3$beta[,k])

```
