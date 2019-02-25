
lars2 <- function(Gram,rhs, type = c("lasso","lar"),
trace = FALSE, eps = .Machine$double.eps, max.steps)
{
    call <- match.call()
    type <- match.arg(type)

    nm <- dim(Gram)
    p <- nm[2]
    im <- inactive <- seq(p)
    namesx <- dimnames(Gram)[[2]]

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
    while((k < max.steps) & (length(active) < (p-length(ignores))))
    {
        action <- NULL
        Cov <- rhs[inactive]
        Cmax <- max(abs(Cov))
        if(Cmax < eps*100){
            if(trace) cat("Max |corr| = 0; exiting...\n")
            break
        }
        k <- k+1
        lambda[k] <- Cmax
        if(!any(drops)){
            new <- abs(Cov) >= Cmax-eps
            Cov <- Cov[!new]
            new <- inactive[new]
            for(inew in new){
                R <- updateR(Gram[inew,inew],R,drop(Gram[inew,active]),Gram=TRUE,eps=eps)
                if(attr(R,"rank")==length(active)){
                  nR <- seq(length(active))
                  R <- R[nR,nR,drop=FALSE]
                  attr(R,"rank") <- length(active)
                  ignores <- c(ignores,inew)
                  action <- c(action,-inew)
                  if(trace)
                    cat("LARS Step",k,":\t Variable", inew,"\tcollinear; dropped for good\n")
                }else{
                  if(first.in[inew] == 0) first.in[inew] <- k
                  active <- c(active,inew)
                  Sign <- c(Sign,sign(rhs[inew]))
                  action <- c(action,inew)
                  if(trace)
                    cat("LARS Step", k, ":\t Variable",inew,"\tadded\n")
                }
            }
        }else action <- -dropid

        Gi1 <- backsolve(R,backsolvet(R,Sign)) # Ginv*1
        dropouts <- NULL
        A <- 1/sqrt(sum(Gi1*Sign))
        w <- A*Gi1

        if((length(active) >= (p-length(ignores)))){
            gamhat <- Cmax/A
        }else{
            a <- drop(w %*% Gram[active, -c(active,ignores),drop=FALSE])
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
        rhs <- rhs-gamhat*Gram[,active,drop=FALSE]%*%w
        Gamrat <- c(Gamrat,gamhat/(Cmax/A))
        arc.length <- c(arc.length,gamhat)
        if(type == "lasso" && any(drops)){
            dropid <- seq(drops)[drops]
            for(id in rev(dropid)){
                if(trace)
                  cat("Lasso Step",k+1,":\t Variable",active[id],"\tdropped\n")
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

    actions <- actions[seq(k)]
    netdf  <- sapply(actions,function(x)sum(sign(x)))
    df <-  cumsum(netdf)
    df <- c(Null=0,df)

    object <- list(call=call,type=type,df=df,lambda=lambda,
        actions=actions,entry=first.in,Gamrat=Gamrat,
        arc.length=arc.length,beta=beta)
    class(object) <- "lars"
    object
}
