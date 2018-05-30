
#######################
## FUNCTION: PBmap() ##
#######################

PBmap <- function(xx, Sigma)
{
  ## Sanity check: input data
  Sigma <- unname(Sigma)
  stopifnot(isSymmetric.matrix(Sigma))
  stopifnot(nrow(Sigma)==length(xx))
  ## if categorical variable
  if(is.factor(xx)){
    stopifnot(nlevels(xx)==2)
    xx <- model.matrix( ~ xx - 1 )[,1]
  }

  ## Step 1: standardize Sigma
  Sigma.inv <- solve(Sigma)
  sigma.sq <- 1/sum(Sigma.inv)
  SS <- Sigma/sigma.sq
  SS.inv <- solve(SS)
  ## SS.tilde
  n <- nrow(Sigma)
  JJ <- matrix(1,n,n)
  SS.tilde <- SS - JJ
  
  ## Step 4: eigen-decomposition
  eig <- eigen(SS.tilde, symmetric=TRUE)
  LL <- diag(eig$values[-n]) #get rid of the last "zero" eigenvalue
  TT <- eig$vectors[,-n]
  ## B-map
  BB <- sqrt(LL) %*% t(TT) %*% SS.inv
  
  ## Step 5: QR decomposition and rotation
  zz <- BB %*% xx
  AA <- cbind(rep(1,n-1),zz) ## column order!!!
  qrstr <- qr(AA)
  QQ <- qr.Q(qrstr)
  RR <- qr.R(qrstr)
  ## P-map
  Rot <- rbind(RR[,2], c(-1,1)*rev(RR[,2])) / sqrt(sum(zz^2))
  PP <- diag(1,n-1) - QQ%*%t(QQ) + QQ%*%Rot%*%t(QQ)
  if(QQ[1]<0){PP=-PP} ## to ensure that PP%*%zz is parallel to rep(+1,n-1)!!!
  
  ### output ###
  return(list("sigma.sq"=sigma.sq, "SS"=SS, "SS.inv"=SS.inv, "SS.tilde"=SS.tilde, "LL"=LL, "TT"=TT, "BB"=BB, "zz"=zz, "QQ"=QQ, "RR"=RR, "Rot"=Rot, "PP"=PP))
}

########################
## FUNCTION: PBtest() ##
########################

PBtest <- function(YY, xx, Sigma=NULL, test="t", id=NULL, weights=NULL, debug=FALSE)
{
  ## Sanity check: input data
  YY <- as.matrix(YY)
  stopifnot(is.numeric(YY) == TRUE)
  stopifnot(nrow(YY) == length(xx))
  ## if categorical variable
  if(is.factor(xx)){
    stopifnot(nlevels(xx)==2)
    xx <- model.matrix( ~ xx - 1 )[,1]
  }
  
  ## estimate Sigma if NULL
  if(is.null(Sigma)){
    if(!is.null(c(id,weights))){
      Sigma <- getSigma(YY, xx, id, weights) #YY is a MATRIX
    }
    else Sigma <- diag(length(xx))
  }

  ## Step 0: calculate PBmap
  out.PBmap <- PBmap(xx, Sigma)

  # ## Step 2: estimate mu.hat
  # mu.hat <- colSums(out.PBmap$SS.inv %*% YY)
  # 
  # ## Step 3: substract mu.hat from YY ==> YY.1
  # YY.1 <- sweep(YY, 2, mu.hat, "-")

  ## Step 4: B-map transforms YY ==> YY.2
  out.PBmap$BB <- sqrt(out.PBmap$LL) %*% t(out.PBmap$TT) %*% out.PBmap$SS.inv
  YY.2 <- out.PBmap$BB %*% YY
  
  ## Step 5: P-map transforms YY.2 ==> YY.tilde
  YY.tilde <- out.PBmap$PP %*% YY.2
  
  ###### TESTS ######
  
  ## calculate degrees of freedom using KR-approximation
  adj <- KRapprox(xx, Sigma) 
  dff <- as.numeric(adj$df)
  scaling <- adj$scaling
  
  ## parametric test
  if(test=="t"){
    stat <- rowttests(t(YY.tilde), tstatOnly=TRUE)$statistic
    stat.adj <- sqrt(scaling)*stat #adjust t-statistic using sqrt of the scaling factor from KRapprox
  }
  
  ## nonparametric test
  if(test=="wilcox"){
    stat <- apply(as.matrix(YY.tilde), 2, FUN = my.wilcox)
    nn <- apply(as.matrix(YY.tilde), 2, FUN = function(z) sum(z!=0))
    mu <- nn*(nn+1)/4
    varr <- nn*(nn+1)*(2*nn+1)/24
    stat.adj <- (stat-mu)/sqrt(varr) #normal approximation
  }
  
  ## p-values
  pval <- sapply(stat.adj, FUN=function(z){2*pt(abs(z), df = dff, lower.tail = F)}) #two-sided
  
  ## assine names
  names(pval) <- names(stat.adj) <- colnames(YY)
  
  ### output ###
  if(debug){
    return(c(list("p.value"=pval, "statistic"=stat.adj, "df"=dff, "YY.tilde"=YY.tilde, "Sigma"=Sigma), out.PBmap))
  }
  else return(list("p.value"=pval, "statistic"=stat.adj, "df"=dff))
}

##################
## some useful function
##################
my.wilcox <- function(x){
  x <- x[x!=0]
  SN <- x>0
  RK <- rank(abs(x))
  sum(SN*RK)
}

