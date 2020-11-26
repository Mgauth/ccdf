#' CDF function
#'
#' @param X
#'
#'
#' @export


CDF <- function(Y,X,Z=NULL,method="logistic regression", fast=FALSE){

  if (class(Z)=="NULL"){
    n_Y <- length(Y)
    temp_order <- sort(Y,index.return=TRUE)$ix
    y <- sort(unique(Y))
    y_sort <- sort(Y)
    x_sort <- X[temp_order]
    modelmat <- model.matrix(Y~X)

    cdf <- list()
    ccdf <- list()

    for (i in 1:length(y)){
      w <- which(Y==y[i])
      indi_Y <- 1*(Y<=y[i])

      # unCDF
      cdf[[i]] <- rep(sum(indi_Y)/n_Y, length(w))

      if (length(unique(indi_Y))==1){
        ccdf[[i]] <- rep(1,length(w))
      }

      else{

        if (method=="RF"){
          # CCDF
          rf <- randomForest(data.frame(X=X, row.names = NULL), indi_Y,ntree=100)
          ccdf[[i]] <- predict(rf)[w]
        }

        else if (method=="CART tree"){
          # CDF
          tree <- rpart(indi_Y~X)
          ccdf[[i]] <- predict(tree, newdata = new_data)
        }

        else if (method=="logistic regression"){
          # CDF
          if(fast){
            #fast
            glm_coef <- RcppNumerical::fastLR(x=modelmat, y=indi_Y,
                                              eps_f = 1e-08, eps_g = 1e-08)$coefficients
          }else{
            #safe
            glm_coef <- glm(indi_Y ~ 1 + X, family = binomial(link = "logit"))$coefficients
          }
          if (is.null(dim(modelmat[w,]))){
            exp_predlin <- exp(-sum(t(glm_coef*t(modelmat[w, ]))))
          }
          else{
            exp_predlin <- exp(-apply(t(glm_coef*t(modelmat[w, ])),1,sum))
          }
          ccdf[[i]] <- 1/(1+exp_predlin) # exp_predlin/(1+exp_predlin)
        }

        else if (method=="linear regression"){
          reg <- lm(indi_Y ~ 1 + modelmat[,2])
          ccdf[[i]] <- predict(reg)[w]
        }
      }
    }

    ccdf <- unlist(ccdf, use.names = FALSE)
    cdf <- unlist(cdf, use.names = FALSE)
    return(list(cdf=cdf, ccdf=ccdf, y=y_sort, x=x_sort))
  }

  else{

    n_Y <- length(Y)
    temp_order <- sort(Y,index.return=TRUE)$ix
    y <- sort(unique(Y))
    y_sort <- sort(Y)
    x_sort <- X[temp_order]
    z_sort <- Z[temp_order]
    modelmat <- model.matrix(Y~X+Z)

    ccdf_x <- list()
    ccdf_nox <- list()
    cdf <- list()

    for (i in 1:length(y)){
      w <- which(Y==y[i])
      indi_Y <- 1*(Y<=y[i])
      new_data <- data.frame(X[w],Z[w])
      names(new_data) <- c("X","Z")

      cdf[[i]] <- rep(sum(indi_Y)/n_Y, length(w))

      if (method=="RF"){
        # CDF
        rf_x <- randomForest(data.frame(X=X, Z=Z, row.names = NULL), indi_Y,ntree=50)
        ccdf_x[[i]] <- predict(rf_x)[w]
        # CCDF
        new_data <- data.frame(Z[w])
        names(new_data) <- c("Z")
        rf_nox <- randomForest(data.frame(Z=Z, row.names = NULL), indi_Y,ntree=50)
        ccdf_nox[[i]] <-  predict(rf_nox)[w]

      }

      if (method=="CART tree"){
        # CDF
        tree_x <- rpart(indi_Y~X+Z)
        ccdf_x[[i]] <- predict(tree_x,newdata = new_data)
        # CCDF
        new_data <- data.frame(Z[w])
        names(new_data) <- c("Z")
        tree_nox <- rpart(indi_Y~Z)
        ccdf_nox[[i]] <- predict(tree_nox,newdata = new_data)

      }

      else if (method=="logistic regression"){
        if (fast){
          glm_coef_x <- RcppNumerical::fastLR(x=modelmat, y=indi_Y)$coefficients
          glm_coef_nox <- RcppNumerical::fastLR(x=modelmat[,-2], y=indi_Y)$coefficients

        }
        else{
          glm_coef_x <- glm.fit(x=modelmat,y=indi_Y, family = binomial())$coefficients
          glm_coef_nox <- glm.fit(x=modelmat[,-2],y=indi_Y, family = binomial())$coefficients
        }
        if (is.null(dim(modelmat[w,]))){
          exp_predlin_x <- exp(-sum(t(glm_coef_x*t(modelmat[w, ]))))
          exp_predlin_nox <- exp(-sum(t(glm_coef_nox*t(modelmat[w,-2]))))
        }
        else{
          exp_predlin_x <- exp(-apply(t(glm_coef_x*t(modelmat[w, ])),1,sum))
          exp_predlin_nox <- exp(-apply(t(glm_coef_nox*t(modelmat[w,-2])),1,sum))
        }
        ccdf_x[[i]] <- 1/(1+exp_predlin_x) # exp_predlin_x/(1+exp_predlin_x)
        ccdf_nox[[i]] <- 1/(1+exp_predlin_nox) # exp_predlin_nox/(1+exp_predlin_nox)

      }
    }

    ccdf_x <- unlist(ccdf_x, use.names = FALSE)
    ccdf_nox <- unlist(ccdf_nox, use.names = FALSE)
    cdf <- unlist(cdf, use.names = FALSE)

    return(list(cdf=cdf, ccdf_nox=ccdf_nox, ccdf_x=ccdf_x, y=y_sort, x_sort=x_sort, z=z_sort))
  }
}


CDF_boot <- function(Y,X,Z=NULL,method="logistic regression", fast=TRUE){

  B <- 1000

  if (!is.null(Z)){
    Z_boot <- Z[boot_index]
  }

  if (class(Z)=="NULL"){

    cdf <- empirical_cdf(Y, ubounds=Y)$CDF

    Y_ordonne <- sort(Y, index.return=TRUE)
    Y <- Y_ordonne$x
    if (class(X)=="factor") X_prime <- factor(rep(X[1], length(X)), levels = levels(X))
    else X_prime <- rep(NA, length(X))


    for (k in 1:length(Y)){
      X_prime[k] <- X[Y_ordonne$ix[k]]
    }
    y <- sort(unique(Y))
    X <- X_prime

    modelmat <- model.matrix(Y~X)

    x <- list()
    ccdf <- list()
    MAT <- list()

    list_OOB <- list()
    list_boot <- list()
    glm_coef_pred <- list()


    # if (method=="RF"){
    #   # CCDF
    #   rf <- randomForest(data.frame(X=X, row.names = NULL), indi_Y, ntree = 50)
    # }


    for (b in 1:B){


      boot_index <- unique(sample(length(Y),length(Y) ,replace = TRUE))
      Y_boot <- Y[boot_index]
      X_boot <- X[boot_index]
      y_oob <- sort(unique(Y[-boot_index]))
      y_boot <- sort(unique(Y[boot_index]))
      list_boot[[b]] <- boot_index
      glm_coef_mat <- matrix(NA,length(y_oob),2)

      list_OOB[[b]] <- sort(unique(Y[-boot_index]))
      pred_OOB <- matrix(NA,length(list_OOB[[b]]),length(list_OOB[[b]]))


      for (i in 1:length(y_oob)){


        indi_Y <- 1*(Y_boot<=y_boot[i])
        #x[[i]] <- X_boot[w]


        if (method=="CART tree"){
          # CDF
          tree <- rpart(indi_Y~X_boot)
          pred_OOB[b,] <- predict(tree, newdata = X[list_OOB[[b]]])
          #ccdf[[i]] <- predict(tree, newdata = X[list_OOB[[b]]])
        }

        else if (method=="logistic regression"){
          # CDF
          if(fast){
            #fast
            #glm_coef <-RcppNumerical::fastLR(x=modelmat[sort(boot_index),], y=indi_Y,
            #eps_f = 1e-08, eps_g = 1e-08)$coefficients
            glm_coef <- RcppNumerical::fastLR(x=modelmat[sort(boot_index),], y=indi_Y)$coefficients
            glm_coef_mat[i,] <- glm_coef

          }else{
            #safe
            glm_coef <- glm(indi_Y ~ 1 + X_boot, family = binomial(link = "logit"))$coefficients
          }

          exp_predlin <- exp(rowSums(glm_coef*modelmat[-list_boot[[b]],]))
          pred_OOB[,i] <- exp_predlin/(1+exp_predlin)
        }
      }
      glm_coef_pred[[b]] <- glm_coef_mat
      MAT[[b]] <- pred_OOB
    }


    pred <- rep(NA,length(y))
    for (i in 1:length(y)){
      pred_courant <- NULL
      for (b in 1:B){
        if (is.element(y[i], list_OOB[[b]])==TRUE){
          w <- which(list_OOB[[b]]==y[i])
          pred_courant <- c(pred_courant,MAT[[b]][w,w])
        }
      }
      pred[i] <- mean(pred_courant)
    }

    ####

    ccdf <- switch(method,
                   "RF" = rf$predicted[i],
                   "CART" = 0, # ???
                   "logistic regression" = pred)

    x <- unlist(x)
    #ccdf <- unlist(ccdf)
    #cdf <- unlist(cdf)
    return(list(y=y, x=x, cdf=cdf, ccdf=ccdf, coef=glm_coef_pred, list_boot=list_boot, list_OOB=list_OOB))
  }

  else{

    y <- sort(unique(Y))
    x <- list()
    z <- list()
    ccdf_x <- list()
    ccdf_nox <- list()

    for (i in 1:length(y)){
      w <- which(Y==y[i])
      indi_Y <- 1*(Y<=y[i])
      x[[i]] <- X[w]
      z[[i]] <- Z[w]
      new_data <- data.frame(X[w],Z[w])
      names(new_data) <- c("X","Z")

      if (method=="RF"){
        # CDF
        rf_x <- randomForest(data.frame(X=X, Z=Z, row.names = NULL), indi_Y,ntree=50)
        ccdf_x[[i]] <- predict(rf_x)[w]
        # CCDF
        new_data <- data.frame(Z[w])
        names(new_data) <- c("Z")
        rf_nox <- randomForest(data.frame(Z=Z, row.names = NULL), indi_Y,ntree=50)
        ccdf_nox[[i]] <-  predict(rf_nox)[w]

      }

      if (method=="CART tree"){
        # CDF
        tree_x <- rpart(indi_Y~X+Z)
        ccdf_x[[i]] <- predict(tree_x,newdata = new_data)
        # CCDF
        new_data <- data.frame(Z[w])
        names(new_data) <- c("Z")
        tree_nox <- rpart(indi_Y~Z)
        ccdf_nox[[i]] <- predict(tree_nox,newdata = new_data)

      }

      else if (method=="logistic regression"){
        # CDF
        reg_x <- glm(indi_Y ~ X+Z, family = binomial(logit))
        ccdf_x[[i]] <- predict(reg_x,type="response",newdata = new_data)
        # CCDF
        new_data <- data.frame(Z[w])
        names(new_data) <- c("Z")
        reg_nox <- glm(indi_Y ~ Z, family = binomial(logit))
        ccdf_nox[[i]] <- predict(reg_nox,type="response",newdata = new_data)


        #####


        # else if (method=="logistic regression"){
        #   # CDF
        #   if(fast){
        #     #fast
        #     glm_coef <- RcppNumerical::fastLR(x=modelmat, y=indi_Y,
        #                                       eps_f = 1e-08, eps_g = 1e-08)$coefficients
        #   }else{
        #     #safe
        #     glm_coef <- glm(indi_Y ~ 1 + X, family = binomial(link = "logit"))$coefficients
        #   }
        #   exp_predlin <- exp(sum(glm_coef*modelmat[w, ]))
        #   ccdf[[i]] <- exp_predlin/(1+exp_predlin)
        # }


      }
    }

    x <- unlist(x)
    z <- unlist(z)
    ccdf_x <- unlist(ccdf_x)
    ccdf_nox <- unlist(ccdf_nox)
    return(list(y=y, x=x, z=z, ccdf_nox=ccdf_nox, ccdf_x=ccdf_x))
  }
}


reg_estimate <- function(Y,X,Z=NULL,method="RF", fast=TRUE, perm=FALSE){

  if (class(Z)=="NULL"){
    mse <- NULL
    if (perm==FALSE){
      rf <- randomForest(data.frame(X=X, row.names = NULL), Y, ntree = 500)
      mse <- rf$mse[500]
      mean((Y-rf$predicted)^2)
    }

    else{
      rf <- randomForest(data.frame(X=X, row.names = NULL), Y, ntree = 50)
      for (b in 1:100){
        pred <- predict(rf,newdata = data.frame(X=sample(X), row.names = NULL))
        mse[b] <- mean((Y-pred)^2)

      }

    }

  }
}

