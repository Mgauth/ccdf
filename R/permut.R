#' permut
#'
#' @param Y
#'
#' @import doParallel, parallel, foreach
#' @export
#'
#' @examples
#' Y <- c(rep(0,40),rnorm(n = 160))
#' X <- as.factor(rbinom(n=200, size = 1, prob = 0.5))
#' Z <- rnorm(n = 100)
#' permut(Y,X,n_perm=100)

permut <- function(Y, X, Z = NULL, distance = "L2", n_perm, method="logistic regression",
                     parallel = TRUE, n_cpus = NULL, fast=TRUE){

  if(parallel){
    if(is.null(n_cpus)){
      n_cpus <- parallel::detectCores() - 1
    }
    cl <- parallel::makeCluster(n_cpus)
    doParallel::registerDoParallel(cl)
  }else{
    cl <- NULL
  }


  if (is.null(Z)){

    res <- CDF(Y, X, method=method, fast=fast)
    w_init <- weights_ccdf(Y,X)

    if (distance=="L2"){
      init_dist <- sqrt(sum(w_init*(res$cdf-res$ccdf)^2))
    }
    if (distance=="L1"){
      init_dist <- sum(w_init*abs(res$cdf-res$ccdf))
    }
    if (distance=="L_sup"){
      init_dist <- max(w_init*abs(res$cdf- res$ccdf))
    }

    results <- foreach(i = 1:n_perm, .combine = 'c') %dopar% {
      X_star <- sample(X)
      res_perm <- CDF(Y=Y, X=X_star, method=method, fast = fast)
      w_perm <- weights_ccdf(Y,X_star)
      switch(distance,
             "L2" = sqrt(sum(w_perm*(res$cdf-res_perm$ccdf)^2)),
             "L1" = sum(w_perm*abs(res$cdf-res_perm$ccdf)),
             "L_sup" = max(w_perm*abs(res$cdf- res_perm$ccdf)))
    }

  }

  else{

    res_init <-  CDF(Y, X, Z, method=method, fast=fast)
    w_init <- weights_ccdf(Y,X,Z)

    init_dist <- sqrt(sum(w_init*(res_init$ccdf_nox-res_init$ccdf_x)^2))

    if (distance=="L2"){
      init_dist <- sqrt(sum(w_init*(res_init$ccdf_nox-res_init$ccdf_x)^2))
    }
    if (distance=="L1"){
      init_dist <- sum(w_init*abs(res_init$ccdf_nox-res_init$ccdf_x))
    }
    if (distance=="L_sup"){
      init_dist <- max(w_init*abs(res_init$ccdf_nox-res_init$ccdf_x))
    }

    if (is.factor(Z)){
      Z_ind <- unlist(lapply(unique(Z),function(z){which(Z==z)}))
    }

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

      res <- CDF(Y,X_star,Z, method=method, fast=fast)
      w_perm <- weights_ccdf(Y,X_star,Z)

      switch(distance,
             L2 = sqrt(sum(w_perm*(res_init$ccdf_nox-res$ccdf_x)^2)),
             L1 = sum(w_perm*abs(res_init$ccdf_nox-res$ccdf_x)),
             L_sup = max(w_perm*abs(res_init$ccdf_nox-res$ccdf_x)))

      #AD = length(Y)*sum(((res$ccdf_nox-res$ccdf_x)^2)/(res$ccdf_nox*(1-res$ccdf_nox)))
    }
  }

  pval <- (sum(1*(results>=init_dist))+1)/(n_perm+1)

  if(parallel){
    parallel::stopCluster(cl)
  }
  return(list(pval=pval))

}
