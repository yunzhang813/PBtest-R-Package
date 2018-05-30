
##########################
## FUNCTION: KRapprox() ##
##########################

###########
## INPUTS: xx, Sigma (covariance matrix), L (testing contrast)
###########

KRapprox <- function(xx, Sigma, L=c(0,1)){
  nn <- length(xx)
  myX <- model.matrix(~xx)
  myPhi <- solve(t(myX)%*%solve(Sigma)%*%myX)
  
  mySigmaG <- list()
  mySigmaG$Sigma <- .makeSparse(Sigma)
  G1 <- .makeSparse(as.matrix(Sigma!=0))
  G2 <- .makeSparse(diag(nn))
  G <- list(G1,G2)
  mySigmaG$G <- G
  mySigmaG$n.ggamma <- 2
  
  myPhiA <- vcovAdj16_internal(myPhi, mySigmaG, myX, details=0)
  
  ## ajusted degrees of freedom and scaling factor
  .KR_adjust(myPhiA, myPhi, L) #values are returned in this function
}



#######################################################################
## The following functions are from the R package 'pbkrtest'
#######################################################################

## --------------------------------------------------------------------
## utility functions
## --------------------------------------------------------------------

.makeSparse<-function(X) {
  X <- as.matrix( X )
  w <- cbind( c(row(X)), c(col(X)), c(X))
  w <- w[ abs( w[,3] ) > 1e-16, ,drop = FALSE]
  Y <- sparseMatrix( w[,1], w[,2], x=w[,3], dims=dim(X))
}

.spur<-function(U){
  sum(diag(U))
}

.divZero = function(x,y,tol=1e-14){
  ## ratio x/y is set to 1 if both |x| and |y| are below tol
  x.y  =  if( abs(x)<tol & abs(y)<tol) {1} else x/y
  x.y
}

.indexSymmat2vec <- function(i,j,N) {
  ## S[i,j] symetric N times N matrix
  ## r the vector of upper triangular element  in row major order:
  ## r= c(S[1,1],S[1,2]...,S[1,j], S[1,N], S[2,2],...S[N,N]
  ##Result: k: index of k-th element of r
  k <-if (i<=j) {
    (i-1)*(N-i/2)+j
  } else {
    (j-1)*(N-j/2)+i
  }
}


## --------------------------------------------------------------------
## derive PhiA
## --------------------------------------------------------------------

vcovAdj16_internal <- function(Phi, SigmaG, X, details=0){

  #    save(SigmaG, file="SigmaG.RData")
  #    return(19)

  details=0
  DB <- details > 0 ## debugging only
  t0 <- proc.time()

  ##Sigma <- SigmaG$Sigma
  n.ggamma <- SigmaG$n.ggamma

  M <- cbind(do.call(cbind, SigmaG$G), X)
  if(DB)cat(sprintf("dim(M) : %s\n",      toString(dim(M)))) ## M can have many many columns
  if(DB)cat(sprintf("dim(SigmaG) : %s\n", toString(dim(SigmaG))))

  if(DB){cat(sprintf("M etc:   %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
  ##SinvM <- solve(SigmaG$Sigma, M, sparse=TRUE)
  SinvM <- chol2inv(chol( forceSymmetric( SigmaG$Sigma ))) %*% M
  ##SigmaInv <- chol2inv( chol( forceSymmetric(SigmaG$Sigma) ) )

  if(DB){cat(sprintf("SinvM etc:   %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}

  v   <- c(rep(1:length(SigmaG$G), each=nrow(SinvM)), rep(length(SigmaG$G)+1, ncol(X)))
  idx <- lapply(unique.default(v), function(i) which(v==i))

  SinvG <- lapply(idx, function(z) SinvM[,z])  ## List of SinvG1, SinvG2,... SinvGr, SinvX
  SinvX <- SinvG[[length(SinvG)]]              ## Kaldes TT andre steder
  SinvG[length(SinvG)] <- NULL                 ## Er HH^t

  if(DB){cat(sprintf("SinvG etc:   %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}

  ##stat <<- list(SigmaG=SigmaG, X=X, M=M)

  OO <- lapply(1:n.ggamma, function(i) {
    SigmaG$G[[i]] %*% SinvX  ## G_i \Sigma\inv X; n \times p
  })

  if(DB){cat(sprintf("Finding OO:   %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}

  PP <- vector("list", n.ggamma)
  QQ <- vector("list", n.ggamma * (n.ggamma + 1) / 2 )
  index <- 1
  for (r in 1:n.ggamma) {
    OOt.r <- t( OO[[ r ]] )
    #str(list("dim(OOt.r)"=dim(OOt.r), "dim(SinvX)"=dim(SinvX)))
    ##PP[[r]] <- forceSymmetric( -1 * OOt.r %*%  SinvX) ## PP : p \times p
    PP[[r]] <- -1 * (OOt.r %*%  SinvX) ## PP : p \times p

    for (s in r:n.ggamma) {
      QQ[[index]] <- OOt.r %*% ( SinvG[[s]] %*% SinvX )
      index <- index + 1;
    }
  }
  ##stat16 <<- list(Phi=Phi, OO=OO, PP=PP,QQ=QQ)

  if(DB){cat(sprintf("Finding PP,QQ:   %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}


  Ktrace <- matrix( NA, nrow=n.ggamma, ncol=n.ggamma )
  for (r in 1:n.ggamma) {
    HHr <- SinvG[[r]]
    for (s in r:n.ggamma){
      Ktrace[r,s] <- Ktrace[s,r] <- sum( HHr * SinvG[[s]] )
    }}

  if(DB){cat(sprintf("Finding Ktrace:   %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}

  ## Finding information matrix
  IE2 <- matrix(0, nrow=n.ggamma, ncol=n.ggamma )
  for (ii in 1:n.ggamma) {
    Phi.P.ii <- Phi %*% PP[[ii]]
    for (jj in c(ii:n.ggamma)) {
      www <- .indexSymmat2vec( ii, jj, n.ggamma )
      IE2[ii,jj]<- IE2[jj,ii] <- Ktrace[ii,jj] -
        2 * sum(Phi * QQ[[ www ]]) + sum( Phi.P.ii * ( PP[[jj]] %*% Phi))
    }}
  if(DB){cat(sprintf("Finding IE2:      %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}

  eigenIE2 <- eigen( IE2, only.values=TRUE )$values
  condi    <- min( abs( eigenIE2 ) )
  WW <- if ( condi > 1e-10 )
    forceSymmetric(2 * solve(IE2))
  else
    forceSymmetric(2 * ginv(IE2))

  ## print("vcovAdj")
  UU <- matrix(0, nrow=ncol(X), ncol=ncol(X))
  ## print(UU)
  for (ii in 1:(n.ggamma-1)) {
    for (jj in c((ii+1):n.ggamma)) {
      www <- .indexSymmat2vec( ii, jj, n.ggamma )
      UU <- UU + WW[ii,jj] * (QQ[[ www ]] - PP[[ii]] %*% Phi %*% PP[[jj]])
    }}
  ## print(UU)

  UU <- UU + t( UU )
  for (ii in 1:n.ggamma) {
    www <- .indexSymmat2vec( ii, ii, n.ggamma )
    UU  <- UU + WW[ii,ii] * (QQ[[ www ]] - PP[[ii]] %*% Phi %*% PP[[ii]])
  }
  if(DB){cat(sprintf("Finding UU:      %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
  ## print(UU)
  GGAMMA <-  Phi %*% UU %*% Phi
  PhiA   <-  Phi + 2 * GGAMMA
  attr(PhiA, "P")     <- PP
  attr(PhiA, "W")     <- WW
  attr(PhiA, "condi") <- condi
  PhiA
}


## --------------------------------------------------------------------
## This is the function that calculates the Kenward-Roger approximation
## --------------------------------------------------------------------
.KR_adjust <- function(PhiA, Phi, L){
  if (!is.matrix(L))
    L = matrix(L, nrow = 1)
  Theta  <-  t(L) %*% solve( L %*% Phi %*% t(L), L)
  P <- attr( PhiA, "P" )
  W <- attr( PhiA, "W" )
  
  A1 <- A2 <- 0
  ThetaPhi <- Theta %*% Phi
  n.ggamma <- length(P)
  for (ii in 1:n.ggamma) {
    for (jj in c(ii:n.ggamma)) {
      e  <- ifelse(ii==jj, 1, 2)
      ui <- ThetaPhi %*% P[[ii]] %*% Phi
      uj <- ThetaPhi %*% P[[jj]] %*% Phi
      A1 <- A1 + e* W[ii,jj] * (.spur( ui ) * .spur( uj ))
      A2 <- A2 + e* W[ii,jj] * sum(ui * t(uj))
    }
  }
  
  q <- rankMatrix(L)
  B <- (1/(2*q)) * (A1+6*A2)
  g <- ( (q+1)*A1 - (q+4)*A2 )  / ((q+2)*A2)
  c1<- g/(3*q+ 2*(1-g))
  c2<- (q-g) / (3*q + 2*(1-g))
  c3<- (q+2-g) / ( 3*q+2*(1-g))
  ##  cat(sprintf("q=%i B=%f A1=%f A2=%f\n", q, B, A1, A2))
  ##  cat(sprintf("g=%f, c1=%f, c2=%f, c3=%f\n", g, c1, c2, c3))
  ###orgDef: E<-1/(1-A2/q)
  ###orgDef: V<- 2/q * (1+c1*B) /  ( (1-c2*B)^2 * (1-c3*B) )
  
  ##EE     <- 1/(1-A2/q)
  ##VV     <- (2/q) * (1+c1*B) /  ( (1-c2*B)^2 * (1-c3*B) )
  EE     <- 1 + (A2/q)
  VV     <- (2/q)*(1+B)
  EEstar <- 1/(1-A2/q)
  VVstar <- (2/q)*((1+c1*B)/((1-c2*B)^2 * (1-c3*B)))
  ##  cat(sprintf("EE=%f VV=%f EEstar=%f VVstar=%f\n", EE, VV, EEstar, VVstar))
  V0<-1+c1*B
  V1<-1-c2*B
  V2<-1-c3*B
  V0<-ifelse(abs(V0)<1e-10,0,V0)
  ##  cat(sprintf("V0=%f V1=%f V2=%f\n", V0, V1, V2))
  
  ###orgDef: V<- 2/q* V0 /(V1^2*V2)
  ###orgDef: rho <-  V/(2*E^2)
  
  rho <- 1/q * (.divZero(1-A2/q,V1))^2 * V0/V2
  df2 <- 4 + (q+2)/ (q*rho-1)          ## Here are the adjusted degrees of freedom.
  
  ###orgDef: F.scaling <-  df2 /(E*(df2-2))
  ###altCalc F.scaling<- df2 * .divZero(1-A2/q,df2-2,tol=1e-12)
  ## this does not work because df2-2 can be about 0.1
  F.scaling <- ifelse( abs(df2 - 2) < 1e-2, 1 , df2 * (1 - A2 / q) / (df2 - 2))
  
  return(list("df"=df2, "scaling"=F.scaling))
}

