#' Complex hypothesis testing using (un)conditional independence test
#'
#'@param exprmat a numeric matrix or data frame of size \code{G x n} containing the
#'preprocessed expressions from \code{n} samples (or cells) for \code{G}
#'genes. Default is \code{NULL}.
#'
#'@param covariate 
#'\code{covariate} must be a numeric or factor vector of size \code{n}
#'containing the model covariates for \code{n} samples.
#'
#'@param variable2test
#'\code{variable2test} must be a numeric or factor vector of size \code{n}
#' containing the variable to be tested (the condition).
#'
#'@param which_test a character string indicating which method to use to
#'compute the test, either \code{'asymptotic'}, \code{'permutations'} or 
#'\code{'dist_permutations'}.
#'Default is \code{'asymptotic'}.
#'
#'@param n_perm the number of permutations. Default is \code{100}.
#'
#'@param adaptive a logical flag indicating whether adaptive permutations
#'should be performed. Default is \code{FALSE}.
#'
#'@param n_perm_adaptive a vector of the increasing numbers of 
#'adaptive permutations when \code{adaptive} is \code{TRUE}. 
#'\code{length(n_perm_adaptive)} should be equal to \code{length(thresholds)+1}. 
#'Default is \code{c(0.1,0.05,0.01)}.
#'
#'@param thresholds a vector of the decreasing thresholds to compute
#'adaptive permutations when \code{adaptive} is \code{TRUE}. 
#'\code{length(thresholds)} should be equal to \code{length(n_perm_adaptive)-1}. 
#'Default is \code{c(100,150,250,500)}.
#'
#'@param distance a character string indicating which distance to use to
#'compute the test, either \code{'L2'}, \code{'L1'} or 
#'\code{'L_sup'}, when \code{method} is \code{'dist_permutations'}, 
#'Default is \code{'L2'}.
#'
#'@param fast a logical flag indicating whether the fast version of
#'the logistic regression should be performed when \code{method} 
#'is \code{'dist_permutations'}
#'Default is \code{TRUE}.
#'
#'@param parallel a logical flag indicating whether parallel computation
#'should be enabled. Default is \code{TRUE}.
#'
#'@param n_cpus an integer indicating the number of cores to be used when
#'\code{parallel} is \code{TRUE}.
#'Default is \code{parallel::detectCores() - 1}.
#'
#'
#'@return A list with the following elements:\itemize{
#'   \item \code{which_test}: a character string carrying forward the value of
#'   the '\code{which_test}' argument indicating which test was performed (either
#'   'asymptotic','permutations','dist_permutations').
#'   \item \code{n_perm}: an integer carrying forward the value of the
#'   '\code{n_perm}' argument or '\code{n_perm_adaptive}' indicating the number of permutations performed
#'   (\code{NA} if asymptotic test was performed).
#'   \item \code{pval}: computed p-values. A data frame with one raw for
#'   each gene, and with 2 columns: the first one '\code{raw_pval}' contains
#'   the raw p-values, the second one '\code{adj_pval}' contains the FDR adjusted p-values
#'   using Benjamini-Hochberg correction.
#' }
#'
#' @export
#' 
#' @examples
#' X <- as.factor(rbinom(n=100, size = 1, prob = 0.5))
#' Y <- replicate(10, ((X==1)*rnorm(n = 50,0,1)) + ((X==0)*rnorm(n = 50,2,1)))
#' Y <- t(Y)
#' Z <- rnorm(n = 100)
#' res <- ccdf_testing(exprmat=Y, variable2test=X, test="asymptotic") # asymptotic test
#' res2 <- ccdf_testing(exprmat=Y, variable2test=X, test="permutations", adaptive=TRUE) # adaptive permutation test


ccdf_testing <- function(exprmat = NULL,
                         variable2test = NULL,
                         covariate = NULL,
                         method = c("logistic regression","linear regression"),
                         distance = c("L2","L1","L_sup"),
                         test = c("asymptotic","permutations","dist_permutations"),
                         n_perm = 100,
                         n_perm_adaptative = c(100,150,250,500),
                         thresholds = c(0.1,0.05,0.01),
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
  
  if ((length(n_perm_adaptative)!=(length(thresholds)+1))){
    warning("length of thresholds + 1 must be equal to length of n_perm_adaptive. \n",
            "Consider using the default parameters.")
  }
  
  if (which_test == "permutations"){
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
      
      for (k in 1:length(thresholds)){
        
        index <- which((res/(perm+1))<thresholds[k])
        
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
      
      for (k in 1:length(thresholds)){
        
        index <- which((res/(perm+1))<thresholds[k])
        
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
  
  if (adaptive == TRUE){
    n_perm <- n_perm_adaptative
  }
  
  if (test == "asymptotic"){
    n_perm <- NA
  }
  
  return(list(which_test = method,
              n_perm = n_perm, 
              pvals = df))
  
  
}
