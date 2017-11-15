#############
## PBmap() ##
#############

PBmap <- function(Delta, Sigma=NULL)
{
  ## check: input data
  if(is.null(Sigma)){Sigma <- diag(length(Delta))}
  Sigma <- unname(Sigma)
  if(isSymmetric.matrix(Sigma)!=TRUE){stop("Sigma must be a symmetric matrix.")}
  if(nrow(Sigma)!=length(Delta)){stop("nrow(Sigma) must be equal to length(Delta).")}
  if(is.numeric(Delta)!=TRUE){stop("Delta must be a numeric vector.")}
  
  ## Step 1: standardize Sigma
  Sigma.inv <- solve(Sigma)
  sigma.sq <- 1/sum(Sigma.inv)
  Sigma0 <- Sigma/sigma.sq
  Sigma0.inv <- solve(Sigma0)
  ## Sigma0.tilde
  n <- nrow(Sigma)
  JJ <- matrix(1,n,n)
  Sigma0.tilde <- Sigma0 - JJ
  
  ## Step 4: eigen-decomposition
  eig <- eigen(Sigma0.tilde, symmetric=TRUE)
  LL <- diag(eig$values[-n]) #get rid of the last "zero" eigenvalue
  TT <- eig$vectors[,-n]
  ## B-map
  BB <- sqrt(LL) %*% t(TT) %*% Sigma0.inv
  
  ## Step 5: QR decomposition and rotation
  cc <- BB%*%Delta
  AA <- cbind(rep(1,n-1),cc) ## column order!!!
  qrstr <- qr(AA)
  QQ <- qr.Q(qrstr)
  RR <- qr.R(qrstr)
  ## P-map
  Rot <- rbind(RR[,2], c(-1,1)*rev(RR[,2])) / sqrt(sum(cc^2))
  PP <- diag(1,n-1) - QQ%*%t(QQ) + QQ%*%Rot%*%t(QQ)
  if(QQ[1]<0){PP=-PP} ## to ensure that PP%*%cc is parallel to rep(+1,n-1)!!!
  
  ### output ###
  return(list("sigma.sq"=sigma.sq, "Sigma0.inv"=Sigma0.inv, "Sigma0.tilde"=Sigma0.tilde, "LL"=LL, "TT"=TT, "BB"=BB, "cc"=cc, "QQ"=QQ, "RR"=RR, "Rot"=Rot, "PP"=PP))
}

##############
## PBtest() ##
##############

PBtest <- function(XX, Delta, Sigma=NULL, test="fast.ttest", debug=FALSE, ...)
{
  ## check: input data
  XX <- as.matrix(XX)
  if(is.numeric(XX) != TRUE){stop("XX must be a numeric matrix.")}
  if(nrow(XX) != length(Delta)){stop("nrow(XX) must be equal to length(Delta).")}

  ## Step 0: calculate PBmap
  out.PBmap <- PBmap(Delta, Sigma)

  ## Step 2: estimate mu.hat
  mu.hat <- colSums(out.PBmap$Sigma0.inv%*%XX)
  
  ## Step 3: substract mu.hat from XX ==> XX.tilde
  XX.tilde <- sweep(XX, 2, mu.hat, "-")

  ## Step 4: transform the problem to an equivalent problem with identity correlation structure ==> YY
  out.PBmap$BB <- sqrt(out.PBmap$LL) %*% t(out.PBmap$TT) %*% out.PBmap$Sigma0.inv
  YY <- out.PBmap$BB%*%XX
  
  ## Step 5: transform the problem to an equivalent problem with constant mean difference under H1 ==> YY.tilde
  YY.tilde <- out.PBmap$PP%*%YY
  
  ## TESTS
  if(test=="fast.ttest"){pval <- rowttests(t(YY.tilde))$p.value; names(pval) <- colnames(XX)}
  if(test=="ttest"){pval <- apply(YY.tilde, 2, FUN = function(z) t.test(z, ...)$p.value)}
  if(test=="wilcox"){pval <- apply(YY.tilde, 2, FUN = function(z) wilcox.test(z, ...)$p.value)}
  
  ### output ###
  if(debug){return(list("p.value"=pval, "sigma.sq"=out.PBmap$sigma.sq, "Sigma0.inv"=out.PBmap$Sigma0.inv, "mu.hat"=mu.hat, "XX.tilde"=XX.tilde, "Sigma0.tilde"=out.PBmap$Sigma0.tilde, "LL"=out.PBmap$LL, "TT"=out.PBmap$TT, "BB"=out.PBmap$BB, "YY"=YY, "cc"=out.PBmap$cc, "QQ"=out.PBmap$QQ, "RR"=out.PBmap$RR, "Rot"=out.PBmap$Rot, "PP"=out.PBmap$PP, "YY.tilde"=YY.tilde))}
  else return("p.value"=pval)
}



