betta_basic <- function(Chats,seChats,X=NA) {  
  if (isTRUE(is.na(X))) { X <- matrix(rep(1,length(Chats)),ncol=1) }
  else { X <- cbind(1,as.matrix(X)) }
  p <- dim(X)[2]
  consider <- !(is.na(Chats) | is.na(seChats))
  Chats <- Chats[consider]; seChats <- seChats[consider]; X <- as.matrix(X[consider,])
  
  likelihood <- function(a) {
    sigsq_u <- a[1]; beta <- as.matrix(a[2:length(a)]); 
    return(-0.5*sum(log(sigsq_u + seChats^2) + (Chats-X%*%beta)^2/(sigsq_u+seChats^2)))
  }
  
  hess <- function(a) {
    sigsq_u <- a[1]; beta <- a[2:length(a)]; 
    dsds <- -0.5*sum(-1/(sigsq_u+seChats^2)^2+2*(Chats-X%*%beta)^2/(sigsq_u+seChats^2)^3)
    dbdb <- apply(-X^2/(sigsq_u+seChats^2),2,sum)
    return(c("sigsq_u"=dsds,"beta"=dbdb))
  }
  mystart <- c(var(Chats),solve(t(X)%*%X)%*%t(X)%*%Chats)
  output <- optim(mystart,likelihood,hess,hessian=FALSE,control=list(fnscale=-1))
  ests <- output$par
  vars <- 1/(-hess(output$par))
  sigsq_u_hat <- ests[1]
  
  mytable <- list()
  mytable$table <- cbind("Estimates"=ests,"Standard Errors"=sqrt(vars),"p-values"=c(1,rep(2,length(ests)-1))*(1-pnorm(abs(ests/sqrt(vars)))))
  rownames(mytable$table)[2] <- "Intercept"
  mytable$estimates <- ests
  mytable$variances <- vars
  mytable$blups <- c(sigsq_u_hat*diag(1/(sigsq_u_hat+seChats^2))%*%(Chats-X%*%ests[2:(p+1)]))
  mytable$y <- Chats
  mytable$yvars <- seChats^2
  mytable$condfits <- c(X%*%ests[2:(p+1)]+mytable$blups)
  mytable$interval <- c(sigsq_u_hat-1.96*sqrt(vars[1]),sigsq_u_hat+1.96*sqrt(vars[1]))
  return(mytable)
}
