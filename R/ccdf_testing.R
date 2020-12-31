#' Complex dis
#'
#' @param exprmat
#'
#'
#' @export
#' @examples
#' X <- as.factor(rbinom(n=100, size = 1, prob = 0.5))
#' Y <- replicate(10, ((X==1)*rnorm(n = 50,0,1)) + ((X==0)*rnorm(n = 50,2,1)))
#' Y <- t(Y)
#' Z <- rnorm(n = 100)
#' res1 <- ccdf_testing(Y,X,test="permutations",n_cpus=16,adaptive=TRUE)
#'
ccdf_testing <- function(exprmat = NULL,
                         variable2test = NULL,
                         covariate = NULL,
                         method = c("logistic regression","linear regression"),
                         distance = c("L2","L1","L_sup"),
                         test = c("asymptotic","permutations","dist_permutations"),
                         n_perm = 100,
                         n_perm_adaptative = c(100,150,250,500),
                         threshold = c(0.1,0.05,0.01),
                         parallel = TRUE,
                         n_cpus = NULL,
                         fast = TRUE,
                         adaptive=FALSE){
  
  # check
  
  stopifnot(is.matrix(exprmat) | is.data.frame(exprmat))
  stopifnot(is.factor(variable2test) | is.numeric(variable2test))
  stopifnot(is.null(dim(variable2test)))
  stopifnot(is.factor(covariate) | is.numeric(covariate) | is.null(covariate))
  stopifnot(is.null(dim(covariate)))
  stopifnot(is.logical(parallel))
  stopifnot(is.logical(fast))
  stopifnot(is.logical(adaptive))
  stopifnot(is.numeric(n_perm))
  
  if (is.data.frame(exprmat)){
    genes_names <- rownames(exprmat)
    exprmat <- as.matrix(exprmat)
  }
  else{
    genes_names <- c(1:nrow(exprmat))
  }
  
  if (sum(is.na(exprmat)) > 1) {
    warning("'y' contains", sum(is.na(y)), "NA values. ",
            "\nCurrently they are ignored in the computations but ",
            "you should think carefully about where do those NA/NaN ",
            "come from...")
    exprmat <- exprmat[complete.cases(exprmat),]
  }
  
  
  # checking for 0 variance genes
  v_g <- matrixStats::rowVars(exprmat)
  if(sum(v_g==0) > 0){
    warning("Removing ", sum(v_g==0), " genes with 0 variance from ",
            "the testing procedure.\n",
            "  Those genes should probably have been removed ",
            "beforehand...")
    exprmat <- exprmat[v_g>0, ]
  }

  
  if (length(method) > 1) {
    method <- method[1]
  }
  stopifnot(method %in% c("logistic regression","linear regression"))
  
  if (length(test) > 1) {
    test <- test[1]
  }
  stopifnot(test %in% c("asymptotic","permutations","dist_permutations"))
  
  if (length(distance) > 1) {
    distance <- distance[1]
  }
  stopifnot(distance %in% c("L2","L1","L_sup"))
  
  if ((length(n_perm_adaptative)!=(length(threshold)+1))){
    warning("length of threshold + 1 must be equal to length of n_perm_adaptive. \n",
            "Consider using the default parameters.")
  }
  
  if (which_test == "permutation") {
    N_possible_perms <- factorial(ncol(exprmat))
    if (n_perm > N_possible_perms) {
      warning("The number of permutations requested 'n_perm' is ",
              n_perm, "which is larger than the total number of ",
              "existing permutations ", N_possible_perms,
              ". Try a lower number for 'n_perm' (currently ",
              "running with 'nperm=", N_possible_perms, "').")
      n_perm <- N_possible_perms
    }
  }
  
  
  
  # parallel
  
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
  
  # test

  if (test=="dist_permutations"){
    
    if (adaptive==TRUE){
      
      print(paste("Computing", n_perm_adaptative[1], "permutations..."))
      
      res <- pbapply::pbsapply(1:nrow(exprmat), FUN=function(i){permut(
        Y = exprmat[i,],
        X = variables2test,
        Z = covariates,
        n_perm = n_perm_adaptative[1],
        parallel = TRUE,
        n_cpus = n_cpus)$score},cl=1)
      perm <- rep(n_perm_adaptative[1],nrow(exprmat))
      
      for (k in 1:length(threshold)){
        
        index <- which((res/(perm+1))<threshold[k])
        
        if (length(index)==0){break}
        
        else{
          
          print(paste("Computing", sum(n_perm_adaptative[1:(k+1)]), "permutations..."))
          
          res_perm <- pbapply::pbsapply(1:nrow(exprmat[index,]), FUN=function(i){permut(
            Y = exprmat[index,][i,],
            X = variables2test,
            Z = covariates,
            n_perm = n_perm_adaptative[k+1],
            parallel = parallel,
            n_cpus = n_cpus)$score},cl=1)
          res[index] <- res[index] + res_perm
          perm <- perm[index] + rep(n_perm_adaptative[k+1],nrow(exprmat[index,]))
        }
        
        
      }
      
      df <- data.frame(raw_pval = res/(perm+1),
                       adj_pval = p.adjust(res/(perm+1), method = "BH"))
    }
    
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
        fast = fast)$pval},cl=1)

      
      df <- data.frame(raw_pval = res,
                       adj_pval = p.adjust(res, method = "BH"))
      
    }
    
  }
  
  
  
  else if (test=="permutations"){
    
    if (adaptive==TRUE){
      
      print(paste("Computing", n_perm_adaptative[1], "permutations..."))
      
      res <- pbapply::pbsapply(1:nrow(exprmat), FUN=function(i){test_perm(
        Y = exprmat[i,],
        X = variables2test,
        Z = covariates,
        n_perm = n_perm_adaptative[1],
        parallel = TRUE,
        n_cpus = n_cpus)$score},cl=1)
      perm <- rep(n_perm_adaptative[1],nrow(exprmat))

      for (k in 1:length(threshold)){

        index <- which((res/(perm+1))<threshold[k])
        
        if (length(index)==0){break}
        
        else{
          
          print(paste("Computing", sum(n_perm_adaptative[1:(k+1)]), "permutations..."))
          
          res_perm <- pbapply::pbsapply(1:nrow(exprmat[index,]), FUN=function(i){test_perm(
            Y = exprmat[index,][i,],
            X = variables2test,
            Z = covariates,
            n_perm = n_perm_adaptative[k+1],
            parallel = parallel,
            n_cpus = n_cpus)$score},cl=1)
          res[index] <- res[index] + res_perm
          perm <- perm[index] + rep(n_perm_adaptative[k+1],nrow(exprmat[index,]))
        }
      
        
      }
      
      df <- data.frame(raw_pval = res/(perm+1),
                       adj_pval = p.adjust(res/(perm+1), method = "BH"))
    }
    
    else{
      
      print(paste("Computing", n_perm, "permutations..."))
      
      res <- do.call("rbind",pbapply::pblapply(1:nrow(exprmat), FUN=function(i){
        test_perm(Y = exprmat[i,],
                  X = variables2test,
                  Z = covariates,
                  n_perm = n_perm,
                  parallel = parallel,
                  n_cpus = cl)},cl=1))
      
      #res <- as.vector(unlist(res))
      
      df <- data.frame(raw_pval = res$raw_pval,
                       adj_pval = p.adjust(res$raw_pval, method = "BH"))
      
    }
    
    
  }
  
  else if (test=="asymptotic"){
    Y <- exprmat
    X <- variables2test
    Z <- covariates
    res <- do.call("rbind",pbapply::pblapply(1:nrow(Y), function(i){test_asymp(Y[i,], X, Z)}, cl=n_cpus))
    df <- data.frame(raw_pval = res$raw_pval, adj_pval = p.adjust(res$raw_pval, method = "BH"), test_statistic = res$Stat)
  }
  
  rownames(df) <- genes_names
  
  return(df)
  
}
