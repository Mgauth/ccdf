#' CDF function
#'
#' @param exprmat
#'
#'
#' @export
#' @examples
#' X <- as.factor(rbinom(n=100, size = 1, prob = 0.5))
#' Y <- matrix(rnorm(n = 1000,0,1),nrow=10,ncol=10)
#' Y <- ((X==1)*rnorm(n = 50,0,1)) + ((X==0)*rnorm(n = 50,2,1))
#' Y <- replicate(10,Y)
#' Z <- rnorm(n = 100)
#' res1 <- ccdf_testing(Y,X,Z,method="logistic regression",test="asymptotic",n_cpus=15)
#'
ccdf_testing <- function(exprmat = NULL,
                         variables2test = NULL,
                         covariates = NULL,
                         method = "logistic regression",
                         distance = "L2",
                         test = c("asymptotic","permutations","dist_permutations"),
                         n_perm = 100,
                         n_perm_start = 100,
                         n_perm_end = 1000,
                         parallel = TRUE,
                         n_cpus = NULL,
                         fast = FALSE,
                         adaptive=FALSE){
  
  if(parallel){
    if(is.null(n_cpus)){
      n_cpus <- parallel::detectCores() - 1
    }
    cl <- parallel::makeCluster(n_cpus)
    doParallel::registerDoParallel(cl)
  }
  else{
    cl <- 1
  }
  
  # if(is.null(exprmat)){
  #   stop("'exprmat' should be specified")
  # }
  #
  # stopifnot(!is.null(variables2test))
  
  # checking for 0 variance genes
  # v_g <- matrixStats::rowVars(exprmat)
  # if(sum(v_g==0) > 0){
  #   warning("Removing ", sum(v_g==0), " genes with 0 variance from ",
  #           "the testing procedure.\n",
  #           "  Those genes should probably have been removed ",
  #           "beforehand...")
  #   y <- y[v_g>0, ]
  # }
  if (test=="dist_permutations"){
    
    if (adaptive==TRUE){ # rajouter verif n_perm
      
      step_perm <- c(n_perm_start, n_perm_start*4, n_perm_end)
      
      print(paste("Computing", step_perm[1], "permutations..."))
      
      res <- pbapply::pbsapply(1:nrow(exprmat), FUN=function(i){permut(
        Y = exprmat[i,],
        X = variables2test,
        Z = covariates,
        distance=distance,
        n_perm = step_perm[1],
        method = method,
        parallel = parallel,
        n_cpus = n_cpus,
        fast = fast)},cl=1)
      res <- as.vector(unlist(res))
      nb_perm <- rep(n_perm_start+1,nrow(exprmat))
      nb <- round((res$pval*(n_perm_start+1))-1)
      #res <- as.vector(unlist(res$pval))
      
      for (k in 2:length(step_perm)){
        
        index <- which(df$adj_pval<(0.2/(k-1)))
        
        print(paste("Computing", step_perm[k], "permutations..."))
        
        res_perm <- pbapply::pbsapply(1:nrow(exprmat[index,]), FUN=function(i){permut(
          Y = exprmat[index,][i,],
          X = variables2test,
          Z = covariates,
          distance=distance,
          n_perm = step_perm[k],
          method = method,
          parallel = parallel,
          n_cpus = n_cpus,
          fast = fast)},cl=1)
        res_perm <- as.vector(unlist(res_perm))
        nb_perm[index] <- rep(step_perm[k]+1,nrow(exprmat[index,]))
        nb <- nb[index] + round((res_perm*(step_perm[k]+1))-1)
        #res_perm <- as.vector(unlist(res_perm$pval))
        
      }
      
      pval <- (nb+1)/nb_perm
      df <- data.frame(raw_pval = pval,
                       adj_pval = p.adjust(pval, method = "BH"))
    }
    
    else{
      
      print(paste("Computing", n_perm, "permutations..."))
      
      res <- pbapply::pbsapply(1:nrow(exprmat), FUN=function(i){permut(
        Y = exprmat[i,],
        X = variables2test,
        Z = covariates,
        distance=distance,
        n_perm = n_perm,
        method = method,
        parallel = parallel,
        n_cpus = n_cpus,
        fast = fast)},cl=1)
      res <- as.vector(unlist(res))
      
      df <- data.frame(raw_pval = res,
                       adj_pval = p.adjust(res, method = "BH"))
      
    }
    
  }
  
  
  
  else if (test=="permutations"){
    
    
    print(paste("Computing", n_perm, "permutations..."))
    
    res <- do.call("rbind",pbapply::pblapply(1:nrow(exprmat), FUN=function(i){
      test_perm(Y = exprmat[i,],
                X = variables2test,
                Z = covariates,
                n_perm = n_perm,
                parallel = parallel,
                n_cpus = cl)},cl=1))
    
    #res <- as.vector(unlist(res))
    
    df <- data.frame(raw_pval = res,
                     adj_pval = p.adjust(res, method = "BH"))
    
    
  }
  
  else if (test=="asymptotic"){
    Y <- exprmat
    X <- variables2test
    Z <- covariates
    res <- do.call("rbind",pbapply::pblapply(1:nrow(Y), function(i){test_asymp(Y[i,], X, Z)}, cl=n_cpus))
    df <- data.frame(raw_pval = res$raw_pval, adj_pval = p.adjust(res$raw_pval, method = "BH"), test_statistic = res$Stat)
  }
  
  return(df)
  
}
