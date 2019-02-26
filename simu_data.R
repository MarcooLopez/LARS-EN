# Simulates data having n observations containing one response and p predictors.
# 'meany' and 'sdy': mean and SD of the response
# 'meanx' and 'sdx': mean and SD of the predictors (all the same)
# 'rhox': correlation between two consequtives predictors: cor(X[j],X[j+1])
# 'rhoxy': correlation between the signal=Xb and the response: cor(Xb,y)
simData <- function(n,p,meanx=3,sdx=2,meany=5,sdy=1.5,rhox=0.5,rhoxy=0.7)
{
  X <- matrix(nrow=n,ncol=p,NA)
  X[,1] <- rnorm(n,mean=meanx,sd=sdx)
  for(j in 2:p){
    z <- rnorm(n,mean=meanx,sd=sdx)
    X[,j] <- rhox*X[,j-1] + sqrt(1-rhox^2)*z
    X[,j] <- scale(X[,j],mean(X[,j])-meanx,scale=FALSE)
  }
  b <- rnorm(p)
  signal <- X%*%b
  error <- rnorm(n=n,mean=mean(signal),sd=sd(signal))
  y <- rhoxy*signal+sqrt(1-rhoxy^2)*error
  y <- (y-mean(y))/(sd(y)/sdy) + meany
  return(list(X=X,y=y))
}
