#' @export
#'
#'
#'
#'
#'
asymp_test <- function(Y,X,Z=NULL,method=method, fast=TRUE){

  if (is.null(Z)){
    n <- length(Y)
    F_Y_X <- ebvcdf(X, Y)
    cor <- matrix(NA,n,n)
    F_Y <- sapply(1:n, function(i){ecdf(Y)(Y[i])})
    F_X <- sapply(1:n, function(i){ecdf(X)(X[i])})
    for (i in 1:n){
      for (j in 1:n){
        if (F_X[i]==1|F_X[i]==0|F_Y[j]==1|F_Y[j]==0){
          cor[i,j] <- 0
        }
        else{
          cor[i,j] <- (((F_Y_X(X[i],Y[j])) - (F_X[i]*F_Y[j]))^2)/((F_X[i]*(1-F_X[i]))*(F_Y[j]*(1-F_Y[j])))
        }
        #cor[i,j] <- (((F_Y_X(X[i],Y[j])) - (F_X[i]*F_Y[j]))^2)
      }
    }
    return((n^{-2})*sum(apply(cor,1,sum)))
  }

  else{
    Z <- rep(1,length(Y))

    n <- length(Y)
    V <- CDF(X,Z,method=method, fast=fast) # F(X/U)
    #V <- NWP(X,Z,0.28)
    W <- CDF(Y,Z,method=method,fast=fast) # F(Y/U)
    #W <- NWP(Y,Z,0.28)

    F_V_W <- ebvcdf(V$ccdf, W$ccdf)
    #F_V_W <- ebvcdf(V, W)

    cor <- matrix(NA,n,n)

    F_V <- sapply(1:n, function(i){ecdf(V$ccdf)(V$ccdf[i])})
    F_W <- sapply(1:n, function(i){ecdf(W$ccdf)(W$ccdf[i])})
    for (i in 1:n){
      for (j in 1:n){
        cor[i,j] <- ((F_V_W(V$ccdf[i],W$ccdf[j])) - (F_V[i]*F_W[j]))^2
      }
    }
    return((n^{-2})*sum(apply(cor,1,sum)))
  }

}
