#' Asymptotic test
#'
#' @param Y
#' @param X
#'
#' @import RcppNumerical
#' @import CompQuadForm
#'
#' @export
#' 
#'


test_asymp <- function(Y, X, Z = NULL, space_y = FALSE, number_y = length(unique(Y)), log = FALSE, keep_zeros = TRUE){
  
  Y <- as.numeric(Y)
  
  if (space_y){
    if (keep_zeros){
      if (log){
        y <- exp(seq(log(ifelse(length(which(Y==0))==0,min(Y),min(Y[-which(Y==0)]))),log(max(Y[-which.max(Y)])),length.out=number_y))
      }
      else{    
        y <- seq(min(Y),max(Y[-which.max(Y)]),length.out=number_y)
      }
    }
    else{
      if (log){
        y <- exp(seq(log(min(Y[-which(Y==0)])),log(max(Y[-which.max(Y)])),length.out=number_y))
      }
      else{    
        y <- seq(min(Y[-which(Y==0)]),max(max(Y[-which.max(Y)])),length.out=number_y)
      }
    }
  }
  else{
    if (keep_zeros){
      y <- sort(unique(Y[-which.max(Y)]))
    }
    else{
      y <- sort(unique(Y[-which(Y==0)]))
    }
  }
  
  # no covariates Z
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
  #nb_fact <- colSums(modelmat[,ind_X])
  
  beta <- matrix(NA,(length(y)-1),length(ind_X))
  indi_pi <- matrix(NA,length(Y),(length(y)-1))
  
  Phi <- (1/length(Y))*(t(as.matrix(modelmat))%*%as.matrix(modelmat))
  H <- (solve(Phi)%*%t(as.matrix(modelmat))) # ginv
  H <- H[ind_X,]
  
  for (i in 1:(length(y)-1)){ # on fait varier le seuil
    indi_Y <- 1*(Y<=y[i])
    indi_pi[,i] <- indi_Y
    #beta[i,] <- (1/length(indi_Y))*rowSums(sapply(1:length(Y),function(i){H[i]*indi_pi[i,]}))
    reg <- lm(indi_Y ~ as.matrix(modelmat[,-1]))
    beta[i,] <- reg$coefficients[ind_X]
  }
  
  beta <- as.vector(beta)
  prop <- colMeans(indi_pi)

  if (is.null(dim(H))){
    H_square <- sum(H^2)
    Sigma <- sapply(1:(length(y)-1), function(i){sapply(1:((length(y)-1)*length(ind_X)), function(j){
      if (i<=j){
        (prop[i]-(prop[j]*prop[i]))
      }
      else{
        (prop[j]-(prop[j]*prop[i]))
      }
    })})
    Sigma <- (1/length(Y))*(H_square*Sigma)
  }
  else{
    # temp_Sigma <-  lapply(1:ncol(H), function(k){sapply(1:nrow(H), function(s){sapply(1:nrow(H), function(r){H[s,k]*H[r,k]})})})
    # sum_temp_Sigma <- temp_Sigma[[1]]
    # for (i in 2:ncol(H)){
    #   sum_temp_Sigma <- sum_temp_Sigma + temp_Sigma[[i]]
    # }
    # 
    # ind_sig <- rep(1:(length(y)-1),length(ind_X))
    # Sigma <- sapply(1:((length(y)-1)*length(ind_X)), function(i){sapply(1:((length(y)-1)*length(ind_X)), function(j){
    #   if (i<=j){
    #     sum_temp_Sigma[floor(i/(length(y))+1),floor(j/(length(y))+1)]*(prop[ind_sig[i]]-(prop[ind_sig[j]]*prop[ind_sig[i]]))
    #   }
    #   else{
    #     sum_temp_Sigma[floor(i/(length(y))+1),floor(j/(length(y))+1)]*(prop[ind_sig[j]]-(prop[ind_sig[j]]*prop[ind_sig[i]]))
    #   }
    # })})
    # Sigma <- (1/length(Y))*Sigma
    
    temp_Sigma <-  lapply(1:ncol(H), function(k){sapply(1:nrow(H), function(s){sapply(1:nrow(H), function(r){H[s,k]*H[r,k]})})})
    sum_temp_Sigma <- temp_Sigma[[1]]
    for (i in 2:ncol(H)){
      sum_temp_Sigma <- sum_temp_Sigma + temp_Sigma[[i]]
    }
    
    ind_sig <- rep(1:(length(y)-1),length(ind_X))
    
    Sigma <- matrix(NA,((length(y)-1)*length(ind_X)),((length(y)-1)*length(ind_X)))
    for (i in 1:((length(y)-1)*length(ind_X))){
      for (j in 1:((length(y)-1)*length(ind_X))){
        if (i<=j){
          Sigma[i,j] <- sum_temp_Sigma[floor(i/(length(y))+1),floor(j/(length(y))+1)]*(prop[ind_sig[i]]-(prop[ind_sig[j]]*prop[ind_sig[i]]))
        }
        else{
          Sigma[i,j] <- sum_temp_Sigma[floor(i/(length(y))+1),floor(j/(length(y))+1)]*(prop[ind_sig[j]]-(prop[ind_sig[j]]*prop[ind_sig[i]]))
        }
      }
    }
    Sigma <- (1/length(Y))*Sigma
  }

  
  # if(length(ind_X)==1){
  #   H_square <- sum(H[ind_X,]^2)
  #   Sigma <- (1/length(Y))*(H_square*Sigma)
  # }
  

  #H_square <- h_i^2
  #Sigma <-  lapply(1:length(H_square),function(i){(H_square[i]*Sigma)})
  #Sigma_sum <- matrix(0,length(ind_X)*(length(y)-1),length(ind_X)*(length(y)-1))
  #for (i in 1:length(H_square)){
  #  Sigma_sum <- Sigma_sum + Sigma[[i]]
  #}
  #Sigma_sum <- (1/length(H_square))*Sigma_sum

  
  #Sigma <- sapply(1:size_X,function(i){(1/length(Y))*(H_square[i]*Sigma[(1+(length(y)*i)-length(y)):(length(y)*i),(1+(length(y)*i)-length(y)):(length(y)*i)])}) 
  
  #decomp <- lapply(1:length(ind_X),function(i){eigen(Sigma[[i]])}) 
  decomp <- eigen(Sigma)
  A <- matrix(0,(length(ind_X)*(length(y)-1)),(length(ind_X)*(length(y)-1)))
  diag(A) <- decomp$values
  z <- (sqrt(length(Y)))*beta
  STAT <- sum(t(z)*z)
  
  param <- list(lim=15000,acc= 5e-04)
  
  pval <- CompQuadForm::davies(q=STAT, lambda=diag(A), lim = param$lim, acc = param$acc)$Qq
  
  
  times <- 2
  while ((pval>1)&(times<11)){
    pval <- CompQuadForm::davies(q=STAT, lambda=diag(A), lim = times*param$lim, acc = param$acc)$Qq
    times <- times*2
  }
  
  if (pval>1){pval<-1}
  
  times <- 0.1
  while ((pval>1)&(times<5e-08)){
    pval <- CompQuadForm::davies(q=STAT, lambda=diag(A), lim = param$lim, acc = times*param$acc)$Qq
    times <- 0.1*times
  }
  
  return(data.frame(raw_pval=pval,Stat=STAT))
  
}


