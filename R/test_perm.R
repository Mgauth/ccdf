#' Permutation test using the test statistic
#'
#' @param Y
#' @param X
#' @param Z
#' @param n_perm
#' @param parallel
#' @param n_cpus
#'
#'
#' @export



test_perm <- function(Y, X, Z=NULL, n_perm=100, parallel = TRUE, n_cpus = NULL){
  
  if(parallel){
    if(is.null(n_cpus)){
      n_cpus <- parallel::detectCores() - 1
    }
    cl <- parallel::makeCluster(n_cpus)
    doParallel::registerDoParallel(cl)
  }else{
    cl <- NULL
  }
  
  y <- sort(unique(Y))
  
  if (is.null(Z)){
    modelmat <- model.matrix(Y~X)
  }
  else{
    modelmat <- model.matrix(Y~X+Z)
  }
  
  beta <- rep(NA,length(y))
  for (i in 1:length(y)){
    indi_Y <- 1*(Y<=y[i])
    reg <- lm(indi_Y ~ 1 + modelmat[,-1])
    beta[i] <- reg$coefficients[2]
  }
  
  z <- sqrt(length(Y))*(beta[-length(y)])
  STAT_obs <- sum(t(z)*z)
  
  if (is.null(Z)){
    
    results <- foreach(i = 1:n_perm, .combine = 'c') %dopar% {
      
      X_star <- sample(X)
      modelmat_perm <- model.matrix(Y~X_star)
      beta_perm<- rep(NA,length(y))
      
      for (i in 1:length(y)){
        indi_Y <- 1*(Y<=y[i])
        reg <- lm(indi_Y ~ 1 + modelmat_perm[,-1])
        beta_perm[i] <- reg$coefficients[2]
      }
      
      z <- sqrt(length(Y))*(beta_perm[-length(y)])
      STAT_perm <- sum(t(z)*z)
      
    }
  }
  
  else{
    
    results <- foreach(i = 1:n_perm, .combine = 'c') %dopar% {
      
      sample_X <- function(X,Z,z){
        X_sample <- rep(NA,length(Z))
        for (i in 1:length(z)){
          X_sample[which(Z==z[i])] <- sample(X[which(Z==z[i])])
        }
        return(X_sample)
      }
      
      X_star <- switch(class(Z),
                       "factor" = sample_X(X,Z,unique(Z)),
                       "numeric" = perm_cont(Y,X,Z))
      
      modelmat_perm <- model.matrix(Y~X_star+Z)
      
      for (i in 1:length(y)){
        indi_Y <- 1*(Y<=y[i])
        reg <- lm(indi_Y ~ 1 + modelmat_perm[,-1])
        beta_perm[i] <- reg$coefficients[2]
      }
      
      z <- sqrt(length(Y))*(beta_perm[-length(y)])
      STAT_perm <- sum(t(z)*z)
    }
  }
  
  if(parallel){
    parallel::stopCluster(cl)
  }
  
  pval <- (sum(1*(results>=STAT_obs))+1)/(n_perm+1)
  
  
  return(data.frame(raw_pval=pval))
  
}


