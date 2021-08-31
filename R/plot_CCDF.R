#' Function to plot the CCDF according to the type of X et Z
#'
#'@param Y a numeric vector of size \code{n} containing the
#'preprocessed expressions from \code{n} samples (or cells).
#'
#'@param X a numeric or factor vector of size \code{n}
#'containing the variable to be tested (the condition to be tested). 
#' 
#'@param Z a numeric or factor vector of size \code{n}
#'containing the covariate. Multiple variables are not allowed.
#'
#'@param method a character string indicating which method to use to
#'compute the CCDF, either \code{'linear regression'}, \code{'logistic regression'}
#' and  \code{'permutations'} or \code{'RF'} for Random Forests.
#'Default is \code{'linear regression'} since it is the method used in the test.
#'
#'@param fast a logical flag indicating whether the fast implementation of
#'logistic regression should be used. Only if \code{'dist_permutations'} is specified.
#'Default is \code{TRUE}.
#'
#'@param space_y a logical flag indicating whether the y thresholds are spaced. 
#'When \code{space_y} is \code{TRUE}, a regular sequence between the minimum and 
#'the maximum of the observations is used. Default is \code{FALSE}.
#'
#'@param number_y an integer value indicating the number of y thresholds (and therefore
#'the number of regressions) to perform the test. Default is \code{length(Y)}.
#' 
#'@param log a logical flag indicating whether the y thresholds are spaced in logarithmic scale. 
#'When \code{log} is \code{TRUE}, a regular sequence between the minimum and
#'the maximum in the logarithmic scale of the observations is used. If the observations
#'are sampled from a count distribution, \code{log} should be \code{TRUE}. Default is \code{FALSE}.
#'
#'@param plot_X a logical flag indicating whether the y thresholds are spaced in logarithmic scale. 
#'
#'@import cowplot
#'
#'@export

plot_CCDF <- function(Y,X,Z=NULL,method="linear regression",fast=TRUE,space_y=FALSE,number_y=length(Y),plot_X=TRUE, plot_Z=FALSE){
  
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
  
  
  if (class(Z[,1])=="NULL"){
    
    df_plot <- data.frame("y" = res$y, "x" = res$x, "cdf" = res$cdf, "ccdf" = res$ccdf)
    
    if (is.factor(X[,1])){
      levels(df_plot$x) <- unique(X[,1])
      df_plot$x <- ordered(df_plot$x, levels = levels(df_plot$x))
      l_X <- length(unique(X[,1]))
      ggplot() + ggtitle(colnames(Y)) +
        geom_step(data = df_plot,aes(x = y, y = cdf, color = viridis(n=(l_X+1))[1]), size = 0.5, linetype="dotted") +
        geom_step(data = df_plot, aes(x = y, y = ccdf, color = as.factor(x)), size = 0.5) +
        scale_color_manual(name = "", labels=c("CDF", paste0("CCDF X=", levels(df_plot$x)[ordered(levels(df_plot$x))])),
                           values = viridis(n=4),
                           guide = guide_legend(override.aes = list(linetype = c("dotted",rep("solid",l_X))))) + xlab("gene expression") +
        ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    }
    else{
      ggplot() + ggtitle(colnames(Y)) +
        geom_step(data = df_plot,aes(x = y, y = cdf, color = viridis(n=2)[1]), size = 0.5, linetype="dotted") +
        geom_point(data = df_plot, aes(x = y, y = ccdf, color =  viridis(n=2)[2]), size = 0.5) +
        scale_color_manual(name = "", labels = c("CDF", "CCDF"),
                           values = c(viridis(n=3)),
                           guide = guide_legend(override.aes = list(linetype = c("dotted","blank"),
                                                                    shape = c(NA,16)))) + xlab("gene expression") +
        ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    }
  }
  
  
  
  else{
    
    res <- CCDF(as.numeric(Y[,1]),X,Z,method,fast,space_y,number_y)
    
    df_plot <- data.frame("y" = res$y, "x" = res$x, "z" = res$z, "cdf" = res$cdf, "ccdf_nox" = res$ccdf_nox,  "ccdf_x" = res$ccdf_x, "shape1"= "20", "shape2" = "3")
    
    
    if (plot_X){
      if (is.factor(X[,1])){
        levels(df_plot$x) <- unique(X[,1])
        df_plot$x <- ordered(df_plot$x, levels = levels(df_plot$x))
        l_X <- length(unique(X[,1]))
        p_cdf <- ggplot() + ggtitle(colnames(Y)) + ####
          geom_step(data = df_plot, aes(x = y, y = cdf, color = viridis(n=(l_X+2))[1]), size = 0.5, linetype="dotted") +
          geom_point(data = df_plot, aes(x = y, y = ccdf_x, color = as.factor(x), shape = shape1), size = 2) +
          geom_point(data = df_plot, aes(x = y, y = ccdf_nox, color = as.factor(x), shape = shape2), size = 2) +
          scale_color_manual(name = "", labels = c("CDF", levels(df_plot$x)[ordered(levels(df_plot$x))]),
                             values = c(viridis(n=(l_X+1))),
                             guide = guide_legend(override.aes = list(linetype = c("dotted",rep("blank",l_X)),
                                                                      shape = c(NA,rep(16,l_X))))) +
          scale_shape_manual(name="CCDF type", labels = c("CCDF_X", "CCDF_noX"),
                             values = c(20,3),
                             guide = guide_legend(override.aes = list(shape = c(20,3)))) +
          xlab("y") +
          ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
        
        # Z
        levels(df_plot$x) <- unique(X[,1])
        df_plot$x <- ordered(df_plot$x, levels = levels(df_plot$x))
        l_X <- length(unique(X[,1]))
        l_Z <- length(unique(Z[,1]))
        p_facet <- ggplot() + ggtitle(colnames(Y)) +
          geom_point(data = df_plot, aes(x = y, y = ccdf_x, color = viridis(n=3)[2]), size = 0.5) +
          geom_point(data = df_plot, aes(x = y, y = ccdf_nox, color = "gold"), size = 0.5) +
          scale_color_manual(name = "CCDF type", label =  c("CCDF_X", "CCDF_noX"),
                             values = c(viridis(n=3)[2], "gold"),
                             guide = guide_legend(override.aes = list(linetype = rep("blank",2),
                                                                      shape = rep(16,2)))) +
          xlab("gene expression") +
          ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(.~as.factor(df_plot$z), nrow=3)
        
        ggdraw() +
          draw_plot(p_cdf, 0, .75, 0.75, .25) +
          draw_plot(p_facet, 0, 0, 1, .75)
        
      }
      
      else{
        ggplot() + ggtitle(colnames(Y)) +
          geom_step(data = df_plot, aes(x = y, y = cdf,  color=viridis(n=3)[2]), size = 0.5, linetype="dotted") +
          geom_point(data = df_plot, aes(x = y, y = ccdf_x, color=viridis(n=3)[1]), shape=16, size = 0.5) +
          geom_point(data = df_plot, aes(x = y, y = ccdf_nox, color='gold'), shape=2, size = 0.5) +
          scale_color_manual(name = "", 
                             labels = c("CDF", "CCDF_X", "CCDF_noX"),
                             values = c(viridis(n=3)),
                             guide = guide_legend(override.aes = list(linetype = c("dotted","blank","blank"),
                                                                      shape = c(NA,16,2)))) +
          xlab("gene expression") +
          ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
      }
    }
    else{
      
      if (plot_Z){
        res_X <- CCDF(as.numeric(Y[,1]),X,Z=NULL,method,fast,space_y,number_y)
        df_X <- data.frame("y" = res_X$y, "x" = res_X$x, "cdf" = res_X$cdf, "ccdf" = res_X$ccdf)
        if (is.factor(Z[,1])){
          # X
          if (is.factor(X[,1])){
            levels(df_X$x) <- c(levels(X[,1]), "CDF")
            df_X$x <- ordered(df_X$x, levels = levels(df_X$x))
            l_X <- length(unique(X[,1]))
            
            p_cdf <- ggplot() + ggtitle(colnames(Y)) +
              geom_step(data = df_X, aes(x = y, y = cdf, color = "CDF"), size = 0.5, linetype="dotted") +
              geom_step(data = df_X, aes(x = y, y = ccdf, color = factor(x, exclude = "CDF")), size = 0.5) +
              scale_color_manual(name = "", limits = levels(df_X$x), #labels=c("CDF", paste0("CCDF X=", levels(df_X$x)[ordered(levels(df_X$x))])),
                                 values = c(viridis(n=4)[-4],"#CC0066"),
                                 guide = guide_legend(override.aes = list(linetype = c(rep("solid",l_X),"dotted")))) + xlab("gene expression") +
              ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
          }
          else{
            p_cdf <- ggplot() + ggtitle(colnames(Y)) +
              geom_step(data = df_X,aes(x = y, y = cdf, color = viridis(n=2)[1]), size = 0.5, linetype="dotted") +
              geom_point(data = df_X, aes(x = y, y = ccdf, color =  viridis(n=2)[2]), size = 0.5) +
              scale_color_manual(name = "", labels = c("CDF", "CCDF"),
                                 values = c(viridis(n=3)),
                                 guide = guide_legend(override.aes = list(linetype = c("dotted","blank"),
                                                                          shape = c(NA,16)))) + xlab("gene expression") +
              ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
          }
          
          # Z
          #cat_X <- unique(X[,1])  #paste0("CCDF X=", levels(df_X$x)[ordered(levels(df_X$x))])
          #cat_X <- c(cat_X, "Marginal on Z")
          levels(df_plot$x) <- c(levels(X[,1]), "Marginal on Z")
          df_plot$x <- ordered(df_plot$x, levels = levels(df_plot$x))
          l_X <- length(unique(X[,1]))
          l_Z <- length(unique(Z[,1]))
          
          p_facet <- ggplot() + ggtitle(colnames(Y)) +
            geom_step(data = df_plot, aes(x = y, y = ccdf_x, color = factor(x, exclude = "Marginal on Z"))) +
            geom_step(data = df_plot, aes(x = y, y = ccdf_nox, color = "Marginal on Z"), size = 0.5, linetype = "dotted") +
            scale_color_manual(name = "CCDF", 
                               values = c(viridis(n=4)[-4],"gold"), limits = levels(df_plot$x),
                               guide = guide_legend(override.aes = list(linetype = c("solid","solid","solid","dotted")))) +
            xlab("gene expression") +
            ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(.~as.factor(df_plot$z), nrow=3, scales = "free_y")
          
          ggdraw() +
            draw_plot(p_cdf, 0, 0.75, 0.75, 0.25) +
            draw_plot(p_facet, 0, 0, 1, 0.75) +
            draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.75), size = 15)
          
        }
        
        else{
          ggplot() + ggtitle(colnames(Y)) +
            geom_step(data = df_plot, aes(x = y, y = cdf,  color=viridis(n=3)[2]), size = 0.5, linetype="dotted") +
            geom_point(data = df_plot, aes(x = y, y = ccdf_x, color=viridis(n=3)[1]), shape=16, size = 0.5) +
            geom_point(data = df_plot, aes(x = y, y = ccdf_nox, color='gold'), shape=2, size = 0.5) +
            scale_color_manual(name = "", 
                               labels = c("CDF", "CCDF_X", "CCDF_noX"),
                               values = c(viridis(n=3)),
                               guide = guide_legend(override.aes = list(linetype = c("dotted","blank","blank"),
                                                                        shape = c(NA,16,2)))) +
            xlab("gene expression") +
            ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
        }
      }
      
      else{
        ggplot() + ggtitle(colnames(Y)) +
          geom_step(data = df_plot, aes(x = y, y = cdf,  color=viridis(n=3)[2]), size = 0.5, linetype="dotted") +
          geom_point(data = df_plot, aes(x = y, y = ccdf_x, color=viridis(n=3)[1]), shape=16, size = 0.5) +
          geom_point(data = df_plot, aes(x = y, y = ccdf_nox, color='gold'), shape=2, size = 0.5) +
          scale_color_manual(name = "", 
                             labels = c("CDF", "CCDF_X", "CCDF_noX"),
                             values = c(viridis(n=3)),
                             guide = guide_legend(override.aes = list(linetype = c("dotted","blank","blank"),
                                                                      shape = c(NA,16,2)))) +
          xlab("gene expression") +
          ylab("value") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
      }
    }
  }
}
