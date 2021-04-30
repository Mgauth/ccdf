#' Permutation test using the test statistic
#'
#' @param Y
#' @param X
#' @param Z
#' @param n_perm
#' @param parallel
#' @param n_cpus
#' @param space_y a logical flag indicating whether the y thresholds are spaced. 
#'When \code{space_y} is \code{TRUE}, a regular sequence between the minimum and 
#'the maximum of the observations is used. Default is \code{FALSE}.
#'
#' @param number_y a integer indicating the number of y thresholds (and therefore
#' the number of regressions) to perform the test. Default is \code{NULL}.
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#'
#' @export

test_perm <- function(Y, X, Z = NULL, n_perm = 100, parallel = FALSE, n_cpus = NULL, space_y = FALSE, number_y = 5, log = FALSE){

  if(parallel){
    if(is.null(n_cpus)){
      n_cpus <- parallel::detectCores() - 1
    }
    cl <- parallel::makeCluster(n_cpus)
    doParallel::registerDoParallel(cl)
  }
  else{
    n_cpus <- 1
  }
  
  Y <- as.numeric(Y)
  
  if (space_y){
    if (log){
      y <- exp(seq(log(min(Y)),log(max(Y)),length.out=number_y))
    }
    else{    
      y <- seq(min(unique(Y)),max(unique(Y)),length.out=number_y)
    }
  }
  
  else{
    y <- sort(unique(Y))
  }
  
  
  if (is.null(Z)){ 
    colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
    modelmat <- model.matrix(~.,data=X)
  }
  
  # with covariates Z
  else{
    colnames(X) <- sapply(1:ncol(X), function(i){paste0('X',i)})
    colnames(Z) <- sapply(1:ncol(Z), function(i){paste0('Z',i)})
    modelmat <- model.matrix(~.,data=cbind(X,Z))
  }
  
  ind_X <- which(substring(colnames(modelmat),1,1)=="X")
  
  beta <- matrix(NA,(length(y)-1),length(ind_X))
  indi_pi <- matrix(NA,length(Y),(length(y)-1))
  
  for (i in 1:(length(y)-1)){ # on fait varier le seuil
    indi_Y <- 1*(Y<=y[i])
    indi_pi[,i] <- indi_Y
    reg <- lm(indi_Y ~ as.matrix(modelmat[,-1]))
    beta[i,] <- reg$coefficients[ind_X]
  }
  
  beta <- as.vector(beta)
  
  z <- sqrt(length(Y))*(beta[-length(y)])
  STAT_obs <- sum(t(z)*z)
  
  if (parallel){
    if (is.null(Z)){
      
      results <- foreach(i = 1:n_perm, .combine = 'c') %dopar% {
        
        X_star <- data.frame(X=X[sample(1:nrow(X)),])
        #colnames(X_star) <- sapply(1:ncol(X), function(i){paste0('X',i)})
        modelmat_perm <- model.matrix(~.,data=X_star)
        beta_perm <- matrix(NA,(length(y)-1),length(ind_X))
        indi_pi <- matrix(NA,length(Y),(length(y)-1))
        
        for (i in 1:(length(y)-1)){ # on fait varier le seuil
          indi_Y <- 1*(Y<=y[i])
          indi_pi[,i] <- indi_Y
          reg <- lm(indi_Y ~ as.matrix(modelmat_perm[,-1]))
          beta_perm[i,] <- reg$coefficients[ind_X]
        }
        
        beta_perm <- as.vector(beta_perm)
        z <- sqrt(length(Y))*(beta_perm[-length(y)])
        STAT_perm <- sum(t(z)*z)
        STAT_perm
      }
    }
    
    else{
      
      if ((ncol(X)!=1)|ncol(Z)!=1){
        stop("X and Z must be univariate. The permutation test will soon be improved to handle multivariate variable of interest and multivariate covariate.")
      }
      
      results <- foreach(i = 1:n_perm, .combine = 'c') %dopar% {
        
        sample_X <- function(X,Z,z){
          X_sample <- rep(NA,length(Z))
          for (i in 1:length(z)){
            X_sample[which(Z==z[i])] <- sample(X[which(Z==z[i])])
          }
          return(X_sample)
        }
        
        X_star <- switch(class(Z[,1]),
                         "factor" = sample_X(X[,1],as.numeric(Z[,1]),unique(as.numeric(Z[,1]))), # vérifier
                         "numeric" = perm_cont(Y,as.numeric(levels(X[,1]))[X[,1]],as.numeric(Z[,1])))
        if (is.factor(X[,1])){
          X_star <- data.frame(X=as.factor(X_star))
        }
        
        X_star <- data.frame(X=X_star)
        #colnames(X_star) <- sapply(1:ncol(X), function(i){paste0('X',i)})
        modelmat_perm <- model.matrix(~.,data=cbind(X_star,Z))
        beta_perm <- matrix(NA,(length(y)-1),length(ind_X))
        indi_pi <- matrix(NA,length(Y),(length(y)-1))
        
        for (i in 1:(length(y)-1)){
          indi_Y <- 1*(Y<=y[i])
          indi_pi[,i] <- indi_Y
          reg <- lm(indi_Y ~ as.matrix(modelmat_perm[,-1]))
          beta_perm[i,] <- reg$coefficients[ind_X]
        }
        
        beta_perm <- as.vector(beta_perm)
        z <- sqrt(length(Y))*(beta_perm[-length(y)])
        STAT_perm <- sum(t(z)*z)
        STAT_perm
      }
    }
    parallel::stopCluster(cl)
    
    score <- sum(1*(results>=STAT_obs))
    pval <- (sum(1*(results>=STAT_obs))+1)/(n_perm+1)
  }
  
  else{
    if (is.null(Z)){
      
      STAT_perm <- rep(NA,n_perm)
      
      for (k in 1:n_perm){
        
        X_star <- data.frame(X=X[sample(1:nrow(X)),])
        #colnames(X_star) <- sapply(1:ncol(X), function(i){paste0('X',i)})
        modelmat_perm <- model.matrix(~.,data=X_star)
        beta_perm <- matrix(NA,(length(y)-1),length(ind_X))
        indi_pi <- matrix(NA,length(Y),(length(y)-1))
        
        for (i in 1:(length(y)-1)){ # on fait varier le seuil
          indi_Y <- 1*(Y<=y[i])
          indi_pi[,i] <- indi_Y
          reg <- lm(indi_Y ~ as.matrix(modelmat_perm[,-1]))
          beta_perm[i,] <- reg$coefficients[ind_X]
        }
        
        beta_perm <- as.vector(beta_perm)
        z <- sqrt(length(Y))*(beta_perm[-length(y)])
        STAT_perm[k] <- sum(t(z)*z)
        
      }
    }
    
    else{
      
      STAT_perm <- rep(NA,n_perm)
      
      for(k in 1:n_perm){
        
        sample_X <- function(X,Z,z){
          X_sample <- rep(NA,length(Z))
          for (i in 1:length(z)){
            X_sample[which(Z==z[i])] <- sample(X[which(Z==z[i])])
          }
          return(X_sample)
        }
        
        X_star <- switch(class(Z[,1]),
                         "factor" = sample_X(X[,1],as.numeric(Z[,1]),unique(as.numeric(Z[,1]))), # vérifier
                         "numeric" = perm_cont(Y,as.numeric(levels(X[,1]))[X[,1]],as.numeric(Z[,1])))
        
        if (is.factor(X[,1])){
          X_star <- data.frame(X=as.factor(X_star))
        }

        #colnames(X_star) <- sapply(1:ncol(X), function(i){paste0('X',i)})
        modelmat_perm <- model.matrix(~.,data=cbind(X_star,Z))
        beta_perm <- matrix(NA,(length(y)-1),length(ind_X))
        indi_pi <- matrix(NA,length(Y),(length(y)-1))
        
        for (i in 1:(length(y)-1)){
          indi_Y <- 1*(Y<=y[i])
          indi_pi[,i] <- indi_Y
          reg <- lm(indi_Y ~ as.matrix(modelmat_perm[,-1]))
          beta_perm[i,] <- reg$coefficients[ind_X]
        }
        
        beta_perm <- as.vector(beta_perm)
        z <- sqrt(length(Y))*(beta_perm[-length(y)])
        STAT_perm[k] <- sum(t(z)*z)
      }
    }
    
    score <- sum(1*(STAT_perm>=STAT_obs))
    pval <- (sum(1*(STAT_perm>=STAT_obs))+1)/(n_perm+1)
  }

  return(data.frame(score=score,raw_pval=pval))
  
}


