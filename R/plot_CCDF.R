#' CDF function
#'
#' @param X
#'
#' @import RcppNumerical
#' @importFrom  randomForest randomForest
#' @import rpart
#' 
#' @export
#' 
#' 


plot_CCDF <- function(Y,X,Z=NULL,method="linear regression", fast=TRUE, space_y=FALSE, number_y=length(Y)){
  
  stopifnot(is.data.frame(Y))
  stopifnot(is.data.frame(X))
  stopifnot(is.data.frame(Z) | is.null(Z))
  stopifnot(is.logical(fast))
  
  genes_names <- colnames(Y)
  
  if (sum(is.na(Y)) > 1) {
    warning("'y' contains", sum(is.na(y)), "NA values. ",
            "\nCurrently they are ignored in the computations but ",
            "you should think carefully about where do those NA/NaN ",
            "come from...")
    Y <- Y[complete.cases(Y)]
  }
  
  
  if (length(method) > 1) {
    method <- method[1]
  }
  stopifnot(method %in% c("linear regression","logistic regression", "CART", "RF"))
  
  
  if (space_y){
    if (is.null(number_y)){
      warning("Missing argument", number_y, ". No spacing is used.")
      space_y <- FALSE
    }
  }
  
  
  res <- CCDF(as.numeric(Y[,1]),X,Z,method,fast,space_y,number_y)
  
  
  if (class(Z)=="NULL"){
    
    df_plot <- data.frame("y" = res$y, "x" = res$x, "cdf" = res$cdf, "ccdf" = res$ccdf)
    
    if (is.factor(X[,1])){
      l_X <- length(unique(X[,1]))
      ggplot() + ggtitle(colnames(Y)) +
        geom_step(data = df_plot,aes(x = y, y = cdf, color = viridis(n=(l_X+1))[1]), size = 1) +
        geom_step(data = df_plot, aes(x = y, y = ccdf, color = as.factor(x)), size = 1) +
        scale_color_manual(name = "", labels = c("CDF", paste0("CCDF X=", unique(x)), "CCDF2"),
                           values = c(viridis(n=(l_X+1))[1],viridis(n=(l_X+1))[-1])) + xlab("y") +
        ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    }
    else{
      ggplot() + ggtitle(colnames(Y)) +
        geom_step(data = df_plot,aes(x = y, y = cdf, color = viridis(n=2)[1]), size = 1) +
        geom_point(data = df_plot, aes(x = y, y = ccdf, color =  viridis(n=2)[2]), size = 1) +
        scale_color_manual(name = "", labels = c("CDF", "CCDF"),
                           values = c(viridis(n=3))) + xlab("y") +
        ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    }
  }
  
  
  else{
    
    df_plot <- data.frame("y" = res$y, "x" = res$x, "z" = res$z, "cdf" = res$cdf, "ccdf_nox" = res$ccdf_nox,  "ccdf_x" = res$ccdf_x)
    
    if (is.factor(X[,1])){
      l_X <- length(unique(X[,1]))
      ggplot() + ggtitle(colnames(Y)) +
        geom_step(data = df_plot, aes(x = y, y = cdf, color = viridis(n=(l_X+2))[1]), size = 1) +
        geom_point(data = df_plot, aes(x = y, y = ccdf_x, color = as.factor(x)), shape=16, size = 2) +
        geom_point(data = df_plot, aes(x = y, y = ccdf_nox, color="black"), shape=2, size = 2) +
        scale_color_manual(name = "", 
                           labels = c("CDF", paste0("CCDF_X X=", unique(x)), "CCDF_noX"),
                           values = c(viridis(n=(l_X+2))[-(l_X+2)],"gold"),
                           guide = guide_legend(override.aes = list(linetype = c(1,0,0,0),
                                                                    shape = c(NA,16,16,2)))) +
        xlab("y") +
        ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))

      
      
    }
    else{
      ggplot() + ggtitle(colnames(Y)) +
        geom_step(data = df_plot, aes(x = y, y = cdf, color = viridis(n=3)[1]), size = 1) +
        geom_point(data = df_plot, aes(x = y, y = ccdf_x, color = viridis(n=3)[2]), shape=16, size = 2) +
        geom_point(data = df_plot, aes(x = y, y = ccdf_nox, color=viridis(n=3)[3]), shape=2, size = 2) +
        scale_color_manual(name = "", 
                           labels = c("CDF", "CCDF_X", "CCDF_noX"),
                           values = c(viridis(n=3)),
                           guide = guide_legend(override.aes = list(linetype = c(1,0,0),
                                                                    shape = c(NA,16,2)))) +
        xlab("y") +
        ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    }
    
    plot(res$y,res$ccdf_nox)
    #return(list(cdf=cdf, ccdf_nox=ccdf_nox, ccdf_x=ccdf_x, y=y_sort, x=x_sort, z=z_sort))
  }
}

