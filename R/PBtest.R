#############
## PBmap() ##
#############

PBmap <- function(XX, Sigma, Delta)
{
  Sigma <- unname(Sigma)
  
  ## check: input data
  if((all.equal(Sigma, t(Sigma))) != TRUE)
  {stop("Error in input data: Sigma should be a symmetric and positive semi-definit matrix.")}
  
  ## Step 1: standardize Sigma
  
  Sigma.inv <- solve(Sigma)
  sigma.sq <- 1/sum(Sigma.inv)
  Sigma0 <- Sigma/sigma.sq
  Sigma0.inv <- solve(Sigma0)
  ## check
  if((all.equal(sum(Sigma-sigma.sq*Sigma0), 0)
      & all.equal(sum(Sigma0.inv), 1)) != TRUE)
  {stop("Error in step 1: standardize Sigma matrix.")}
  
  ## Step 2: estimate mu.hat (weighted average, i.e. MLE of XX ~ mu.hat*vec(1))
  
  mu.hat <- sum(Sigma0.inv%*%XX)
  
  ## Step 3: substract mu.hat from XX ==> XX.tilde
  
  XX.tilde <- XX - mu.hat
  
  n <- length(XX)
  JJ <- matrix(1,n,n)
  Sigma0.tilde <- Sigma0 - JJ
  Sigma.tilde <- sigma.sq*Sigma0.tilde
  
  ## Step 4: eigen-decomposition (transform the problem to an equivalent problem with identity correlation structure) ==> YY
  
  eig <- eigen(Sigma0.tilde)
  ## check: last eigenvalue is zero
  if((all.equal(eig$values[length(eig$values)], 0)
      # & all.equal(sigma.sq*(sum(eig$values)+sum(eigen(JJ)$val))-sum(eigen(Sigma)$val), 0)
  ) != TRUE)
  {stop("Error in step 4: eigen-decomposition. Last eigenvalue is not zero.")}
  
  LL <- diag(eig$values[-n]) #get rid of the last "zero" eigenvalue
  TT <- eig$vectors[,-n]
  ## check:
  if((all.equal(dim(TT), c(n,n-1))
      & all.equal(dim(LL), c(n-1,n-1))
      & all.equal(TT %*% LL %*% t(TT), Sigma0.tilde)) != TRUE)
  {stop("Error in step 4: eigen-decomposition. (Problem with LL and TT matrix).")}
  
  BB <- sqrt(LL) %*% t(TT) %*% Sigma0.inv
  ## check
  if((all.equal(dim(BB), c(n-1,n))
      & all.equal(BB%*%Sigma0%*%t(BB), diag(1,n-1))) != TRUE)
  {stop("Error in step 4: construct BB matrix.")}
  
  ## YY
  YY <- BB%*%XX
  
  ## Step 5: QR decomposition and rotation (transform the problem to an equivalent problem with constant mean difference under H1) ==> YY.tilde
  
  cc <- BB%*%Delta
  AA <- cbind(rep(1,n-1),cc) ## column order!!!
  
  qrstr <- qr(AA)
  QQ <- qr.Q(qrstr)
  RR <- qr.R(qrstr)
  ## check
  if((all.equal(abs(RR[1]) - sqrt(n-1), 0)
      & all.equal(as.numeric(abs(RR[1,2]) - abs(AA[,1]%*%AA[,2]/sqrt(n-1))), 0)
      & all.equal(as.numeric(abs(RR[2,2]) - sqrt(AA[,2]%*%AA[,2]-(RR[1,2])^2)), 0)
      & all.equal(as.numeric(AA[,2]%*%AA[,2] - sum(cc^2)), 0)) != TRUE)
  {stop("Error in step 5: QR decomposition.")}
  
  Rot <- rbind(RR[,2], c(-1,1)*rev(RR[,2])) / sqrt(sum(cc^2))
  PP <- diag(1,n-1) - QQ%*%t(QQ) + QQ%*%Rot%*%t(QQ)
  if(QQ[1]<0){PP=-PP} ## to ensure that PP%*%cc is parallel to rep(+1,n-1)!!!
  
  ## check
  if((all.equal(dim(PP), c(n-1,n-1))
      & all.equal(as.numeric(PP%*%cc * sqrt(n-1) / sqrt(sum(cc^2))), rep(1,n-1))
      & all.equal(t(PP)%*%PP, diag(1,n-1))
      & all.equal(PP%*%t(PP), diag(1,n-1))) != TRUE)
  {stop("Error in step 5: construct rotation matrix / PP matrix.")}
  if(sum(PP%*%cc < 0)>0){stop("Error in step 5: PP%*%cc returns negative scalar.")}
  
  ## YY.tilde
  YY.tilde <- drop(PP%*%YY)
  
  ### output ###
  return(list("sigma.sq"=sigma.sq, "Sigma0"=Sigma0, "mu.hat"=mu.hat, "XX.tilde"=XX.tilde, "Sigma0.tilde"=Sigma0.tilde, "LL"=LL, "TT"=TT, "BB"=BB, "YY"=YY, "cc"=cc, "QQ"=QQ, "RR"=RR, "Rot"=Rot, "PP"=PP, "YY.tilde"=YY.tilde))
}

##############
## PBtest() ##
##############

PBtest <- function(XX, Sigma, Delta, test="t", ...){
  YY.tilde <- PBmap(XX, Sigma, Delta)$YY.tilde
  if(test=="t"){pval <- t.test(YY.tilde, ...)$p.value}
  if(test=="wilcox"){pval <- wilcox.test(YY.tilde, ...)$p.value}
  return(pval)
}







