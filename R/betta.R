betta <- function(chats,ses,X=NA) {
  if (isTRUE(is.na(X))) { X <- matrix(rep(1,length(chats)),ncol=1) }
  #else { X <- cbind(1,as.matrix(X)) }
  consider <- !(is.na(chats) | is.na(ses) | apply(is.na(X),1,sum))
  chats <- chats[consider]; ses <- ses[consider]; X <- as.matrix(X[consider,])
  n <- dim(X)[1]; p <- dim(X)[2]
  
  likelihood <- function(input) {
    ssq_u <- input[1]
    beta <- input[2:length(input)]
    W <- diag(1/(ssq_u+ses^2))
    -0.5*(sum(log(ssq_u+ses^2)+(chats-X%*%beta)^2/(ssq_u+ses^2)) + log(det(t(X)%*%W%*%X)))
  }
  
  mystart <- c(var(chats),solve(t(X)%*%X)%*%t(X)%*%chats)
  output <- optim(mystart,likelihood,hessian=FALSE,control=list(fnscale=-1)) # fnscale => maximises if -ve
  ssq_u <- output$par[1] 
  beta <- output$par[2:length(output$par)]
  
  W <- diag(1/(ssq_u+ses^2))
  vars <- 1/diag(t(X)%*%W%*%X)
  
  global <- t(beta)%*%(t(X)%*%W%*%X)%*%beta ## global test
  
  Q <- sum((chats-X%*%beta)^2/ses^2)
  
  R <- diag(ses^2); G <- diag(ssq_u,n)
  
  C <- matrix(NA,nrow=n+p,ncol=n+p)
  C[1:p,1:p] <- t(X)%*%solve(R)%*%X
  C[(p+1):(n+p),(p+1):(n+p)] <- solve(R)+solve(G)
  C[1:p,(p+1):(n+p)] <- t(X)%*%solve(R)
  C[(p+1):(n+p),1:p] <- solve(R)%*%X
  C <- solve(C)
  
  mytable <- list()
  mytable$table <- cbind("Estimates"=beta,"Standard Errors"=sqrt(vars),"p-values"=round(2*(1-pnorm(abs(beta/sqrt(vars)))),3))
  mytable$cov <- solve(t(X)%*%W%*%X)
  mytable$ssq_u <- ssq_u
  mytable$homogeneity <- c(Q,1-pchisq(Q,n-p))
  mytable$global <- c(global,1-pchisq(global,p-1))
  us <-  c(ssq_u*W%*%(chats-X%*%beta))
  mytable$blups <- c(X%*%beta + us)
  blupvars <- cbind(X,diag(1,n))%*%C%*%t(cbind(X,diag(1,n)))
  mytable$blupses <- c(sqrt(diag(blupvars)))
  return(mytable)  
}

betta_pic <- function(y,se,x,ylimu,myy=NA,mymain=NA,mycol=rep("black",length(y)),labs=NA,mypch=rep(16,length(y)),myxlim=c(0,1.1*max(x,na.rm=T))) {
  n <- length(y)
  par(xpd=NA)
  plot(0,0,type="n",xlim=myxlim,ylim=c(0,ylimu),xaxt="n",xlab="",bty="n",ylab=myy,main=mymain)
  for (i in 1:n) {
    if(!is.na(y[i]) & !is.na(x[i])) {
      points(x[i],y[i],pch=mypch[i],col=mycol[i])
      lines(c(x[i],x[i]),c(max(0,y[i]-1.96*se[i]),y[i]+1.96*se[i]),col=mycol[i])
    }
  }
  if(!is.na(labs)) axis(1,at=1:length(y),labels=labs,las=2,cex=0.8)
}