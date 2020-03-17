
##########################
## FUNCTION: getSigma() ##
##########################

###########
## INPUTS: YY (vector or matrix), xx, id, weights, correction=TRUE
###########

getSigma <- function(YY, xx, id, weights=NULL, correction=TRUE){
  ## if unweighted
  if(is.null(weights)) weights <- rep(1,nrow(YY))
  
  rho.hat <- getRho(YY, xx, id, weights, correction)
  mat <- diag(length(id))
  pos <- get.pos(id)
  mat[pos] <- rho.hat
  WW <- diag(1/sqrt(weights))
  Sigma.hat <- WW %*% mat %*% WW
  return(Sigma.hat)
}


######################
## Useful functions ##
######################

## split a data matrix into a list by columns
my.split <- function(mat){
  split(mat, rep(1:ncol(mat), each = nrow(mat)))
}

## get positions of the common correlation paramter (rho) in the structure matrix (Sigma)
get.pos <- function(id){
  pos.list <- lapply(id,function(z) which(id %in% z))
  col.pos <- rep(1:length(pos.list), unlist(lapply(pos.list,length)))
  row.pos <- unlist(pos.list)
  matched.pos <- (col.pos-1)*length(id)+row.pos
  diag.pos <- (0:(length(id)-1))*length(id)+(1:length(id))
  rho.pos <- setdiff(matched.pos, diag.pos)
  rho.pos
}

## calculate the common correlation (i.e. rho)
getRho <- function(YY, xx, id, weights, correction=TRUE, return.all=FALSE){
  ## NOTE: YY can be a MATRIX!!!
  
  ## regress out weights and get residuals from WLS
  res.lm <- lm(YY~xx, weights=weights)$residuals #residuals from WLS
  res.lm <- diag(sqrt(weights))%*%res.lm #reweight residuals
  ## subset data for only correlated data
  tab <- table(id)
  id.selected <- names(tab)[tab>1]
  res.lm1 <- as.matrix(res.lm[id%in%id.selected,])
  id1 <- id[id%in%id.selected]
  ## order residuals in blocks
  res.lm1.o <- as.matrix(res.lm1[order(id1),])
  id1.o <- id1[order(id1)]
  
  ## getRho.each()
  getRho.each <- function(res.lm1.o.each, id1.o){
    res.list <- split(res.lm1.o.each, id1.o)
    SS1 <- sum(unlist(res.list)^2)
    SS2 <- sum(sapply(res.list, function(z) sum(z%o%z)))
    n.l <- sapply(res.list, length)
    r <- (SS2-SS1) / (sum(n.l*(n.l-1))/sum(n.l)) / SS1
    ## rho.hat.each
    if(correction){ r*(1 + (1-r^2)/2/(length(unique(id1.o))-3)) }
    else r
  }
  
  ## average over individal rhos
  rho.hat.all <- sapply(my.split(res.lm1.o), getRho.each, id1.o=id1.o)
  rho.hat <- mean(rho.hat.all, na.rm = TRUE)
  
  ## output
  if(return.all) return(rho.hat.all)
  else return(rho.hat)
}




