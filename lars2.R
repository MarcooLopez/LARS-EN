# Last modification Feb-27-2019.   00:33:00 hrs
# Adapted from 'lars' function in LARS package (Hastie & Efron)
# Dependending on some functions from LARS package
require(lars)

lars2 <- function(XtX,Xty, type = c("lasso","lar"),
verbose = FALSE, eps = .Machine$double.eps, max.steps,standardize=TRUE)
{
    call <- match.call()
    type <- match.arg(type)
    nm <- dim(XtX)
    p <- nm[2]
    im <- inactive <- seq(p)
    namesx <- dimnames(XtX)[[2]]
    s.p <- nchar(p)
    textPrint <- c(" Step","\tSec/Step","\tVariable")

    if(standardize){
      sdx <- sqrt(diag(XtX))
      W <- replicate(p,1/sdx)*t(replicate(p,1/sdx))
      XtX <- XtX*W
      Xty <- Xty/sdx
    }else sdx <- rep(1,p)

    ignores <- NULL
    if(missing(max.steps))  max.steps <- 8*p
    beta <- matrix(0,max.steps+1,p)
    lambda <- double(max.steps)
    Gamrat <- NULL
    arc.length <- NULL
    first.in <- integer(p)
    active <- NULL
    actions <- as.list(seq(max.steps))
    drops <- FALSE
    Sign <- NULL
    R <- NULL
    k <- 0
    time = proc.time()[3]
    while((k < max.steps) & (length(active) < (p-length(ignores))))
    {
        action <- NULL
        Cov <- Xty[inactive]
        Cmax <- max(abs(Cov))
        if(Cmax < eps*100){
            if(verbose) cat("Max |corr| = 0; exiting...\n")
            break
        }
        k <- k+1
        lambda[k] <- Cmax
        if(!any(drops)){
            new <- abs(Cov) >= Cmax-eps
            Cov <- Cov[!new]
            new <- inactive[new]
            for(inew in new){
                R <- updateR(XtX[inew,inew],R,drop(XtX[inew,active]),Gram=TRUE,eps=eps)
                if(attr(R,"rank")==length(active)){
                  nR <- seq(length(active))
                  R <- R[nR,nR,drop=FALSE]
                  attr(R,"rank") <- length(active)
                  ignores <- c(ignores,inew)
                  action <- c(action,-inew)
                  if(verbose){
                    cat("LARS Step",k,":\t Variable", inew,"\tcollinear; dropped for good\n")
                  }
                }else{
                  if(first.in[inew] == 0) first.in[inew] <- k
                  active <- c(active,inew)
                  Sign <- c(Sign,sign(Xty[inew]))
                  action <- c(action,inew)
                  if(verbose){
                    cat("--------------------------------------------------------------\n")
                    tmp <- proc.time()[3]
                    cat(paste(textPrint,"=",c(sprintf("%*d",s.p,k),sprintf('%.*f',4,tmp-time),
                      sprintf("%*d",s.p,inew))),"added\n")
                    time <- tmp
                  }
                }
            }
        }else action <- -dropid
        Gi1 <- backsolve(R,backsolvet(R,Sign))
        A <- 1/sqrt(sum(Gi1*Sign))
        w <- A*Gi1

        if((length(active) >= (p-length(ignores)))){
            gamhat <- Cmax/A
        }else{
            a <- drop(w %*% XtX[active, -c(active,ignores),drop=FALSE])
            gam <- c((Cmax-Cov)/(A-a),(Cmax+Cov)/(A+a))
            gamhat <- min(gam[gam > eps],Cmax/A)
        }
        if(type == "lasso"){
            dropid <- NULL
            b1 <- beta[k,active]
            z1 <- -b1/w
            zmin <- min(z1[z1 > eps],gamhat)
            if(zmin < gamhat){
                gamhat <- zmin
                drops <- z1 == zmin
            }else drops <- FALSE
        }
        beta[k+1,] <- beta[k,]
        beta[k+1,active] <- beta[k+1,active] + gamhat*w
        Xty <- Xty - gamhat*XtX[,active,drop=FALSE]%*%w
        Gamrat <- c(Gamrat,gamhat/(Cmax/A))
        arc.length <- c(arc.length,gamhat)
        if(type == "lasso" && any(drops)){
            dropid <- seq(drops)[drops]
            for(id in rev(dropid)){
                if(verbose){
                  cat("--------------------------------------------------------------\n")
                  tmp <- proc.time()[3]
                  cat(paste(textPrint,"=",c(sprintf("%*d",s.p,k+1),sprintf('%.*f',4,tmp-time),
                    sprintf("%*d",s.p,active[id]))),"dropped\n")
                  time <- tmp
                }
                R <- downdateR(R,id)
            }
            dropid <- active[drops]
            beta[k+1,dropid] <- 0
            active <- active[!drops]
            Sign <- Sign[!drops]
        }
        if(!is.null(namesx))
            names(action) <- namesx[abs(action)]
        actions[[k]] <- action
        inactive <- im[-c(active, ignores)]
    }
    beta <- beta[seq(k+1), ,drop = FALSE]
    lambda  <-  lambda[seq(k)]
    dimnames(beta) <- list(paste(0:k),namesx)
    beta <- scale(beta,FALSE,sdx)
    actions <- actions[seq(k)]
    netdf  <- sapply(actions,function(x)sum(sign(x)))
    df <-  cumsum(netdf)
    df <- c(Null=0,df)

    object <- list(call=call,type=type,df=df,lambda=lambda,
        actions=actions,entry=first.in,Gamrat=Gamrat,
        arc.length=arc.length,beta=beta,sdx=sdx)
    class(object) <- "lars"
    object
}
