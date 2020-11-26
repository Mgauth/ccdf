#' CDF_permut function
#'
#' @param X
#'
#'
#' @export
CDF_permut <- function(Y,X,Z=NULL,method="logistic regression", fast=TRUE, coef=CDF$coeff, list_boot=list_boot, list_OOB=list_OOB){

  B <- 500

  if (!is.null(Z)){
    Z_boot <- Z[boot_index]
  }

  if (class(Z)=="NULL"){


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

    # if (method=="RF"){
    #   # CCDF
    #   rf <- randomForest(data.frame(X=X, row.names = NULL), indi_Y, ntree = 50)
    # }


    for (b in 1:B){

      pred_OOB <- matrix(NA,length(list_OOB[[b]]),length(list_OOB[[b]]))
      y_oob <- sort(unique(Y[-list_boot[[b]]]))

      for (i in 1:length(y_oob)){

        if (method=="CART tree"){
          # CDF
          tree <- rpart(indi_Y~X_boot)
          pred_OOB[b,] <- predict(tree, newdata = X[list_OOB[[b]]])
          #ccdf[[i]] <- predict(tree, newdata = X[list_OOB[[b]]])
        }

        else if (method=="logistic regression"){
          glm_coef <- coef
          exp_predlin <- exp(rowSums(glm_coef[[b]][i,]*modelmat[-list_boot[[b]],]))
          pred_OOB[,i] <- exp_predlin/(1+exp_predlin)
        }
      }
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
    return(list(y=y, x=x, ccdf=ccdf))
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
