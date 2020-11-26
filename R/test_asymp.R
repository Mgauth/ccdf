#' CDF function
#'
#' @param Y
#' @param X
#'
#'
#' @export



test_asymp <- function(Y,X,Z=NULL){

  y <- sort(unique(Y))
  if (is.null(Z)){
    #modelmat <- model.matrix(Y~X)
    modelmat <- cbind(rep(1,length(Y)),X)
  }
  else{
    modelmat <- model.matrix(Y~X+Z)
  }
  beta <- rep(NA,length(y))
  indi_pi <- matrix(NA,length(Y),length(y))

  Phi <- (1/length(Y))*(t(modelmat)%*%modelmat)
  H <- (solve(Phi)%*%t(modelmat))

  for (i in 1:length(y)){ # on fait varier le seuil
    indi_Y <- 1*(Y<=y[i])
    indi_pi[,i] <- indi_Y
    reg <- lm(indi_Y ~ 1 + modelmat[,-1])
    beta[i] <- reg$coefficients[2]
  }


  prop <- colMeans(indi_pi)
  Sigma <- sapply(1:(length(y)-1), function(i){sapply(1:(length(y)-1), function(j){
    if (i<=j){
      prop[i]-(prop[j]*prop[i])
    }
    else{
      prop[j]-(prop[j]*prop[i])
    }
  })})


  H_square <- sum(H[2,]^2)

  Sigma <- (1/(length(Y)-1))*(H_square*Sigma)
  decomp <- eigen(Sigma)
  A <- matrix(0,length(y)-1,length(y)-1)
  diag(A) <- decomp$values
  z <- sqrt(length(y)-1)*(beta[-length(y)])
  STAT <- sum(t(z)*z)

  pval <- CompQuadForm::davies(q=STAT, lambda=diag(A), lim = 15000, acc = 5e-04)$Qq

  return(data.frame(raw_pval=pval,Stat=STAT))

}
