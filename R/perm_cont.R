#' permutation Z continuous
#'
#' @param X
#'
#'
#' @export

perm_cont <- function(Y,X,Z){
  prob <- matrix(NA,length(Z),length(Z))
  new_prob <- matrix(NA,length(Z),length(Z))
  reg <- lm(X~Z)
  X_star <- rep(NA,length(X))
  n <- length(Z)

  for (z in 1:length(Z)){
    for (l in 1:length(Z)){
      prob[z,l] <- abs(predict(reg,newdata=data.frame(Z=Z[z]))-predict(reg,newdata=data.frame(Z=Z[l])))/
          (sum(abs(predict(reg,newdata=data.frame(Z=Z[z]))-predict(reg,newdata=data.frame(Z=Z)))))

    }
    new_prob[z,-z] <- (1/prob[z,-z])/sum(1/prob[z,-z])
    new_prob[z,z] <- 0
  }

  for (i in 1:length(X)){
    X_star[i] <- sample(X[-i], size=1, prob = prob[i,-i])
  }

  return(X_star)
}
