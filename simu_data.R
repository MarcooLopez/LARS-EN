# Simulates data having n observations containing one response and p predictors.
# 'meany' and 'sdy': mean and SD of the response
# 'meanx' and 'sdx': vectors with mean and SD of the predictors (if length=1: all the same)
# 'rhox': correlation between two consequtives predictors: cor(X[j],X[j+1])
# 'rhoxy': correlation between the signal=Xb and the response: cor(Xb,y)
simData <- function(n,p,meanx=rchisq(p,5),sdx=rchisq(p,1),meany=5,sdy=1.5,rhox=0.5,rhoxy=0.7,seed=123)
{
  if(length(meanx)==1)  meanx <- rep(meanx,p)
  if(length(sdx)==1)  sdx <- rep(sdx,p)
  set.seed(seed)
  X <- matrix(nrow=n,ncol=p,NA)
  X[,1] <- rnorm(n,mean=meanx[1],sd=sdx[1])
  for(j in 2:p){
    z <- rnorm(n,mean=meanx[j-1],sd=sdx[j-1])
    X[,j] <- rhox*X[,j-1] + sqrt(1-rhox^2)*z
    X[,j] <- sdx[j]*scale(X[,j])+meanx[j]
  }  
  b <- rnorm(p)
  signal <- X%*%b
  error <- rnorm(n=n,mean=mean(signal),sd=sd(signal))
  y <- rhoxy*signal+sqrt(1-rhoxy^2)*error
  y <- sdy*scale(y)+meany
  return(list(X=X,y=y))
}
