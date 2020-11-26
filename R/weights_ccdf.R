#' weights
#'
#' @param X
#'
#'
#' @export


weights_ccdf <- function(Y,X,Z=NULL){

  find_unique <- function(couple, unique_couple){
    index_list <- list()
    for (j in 1:nrow(unique_couple)){
      index <- NULL
      for (i in 1:nrow(couple)){
        if (identical(as.numeric(couple[i,]),as.numeric(unique_couple[j,]))){
          index <- c(index,i)
          index_list[[j]] <- index
        }
      }
    }
    return(index_list)
  }

  temp_order <- sort(Y,index.return=TRUE)$ix

  if (is.null(Z)){
    couple <- data.frame(Y=Y[temp_order],X=X[temp_order])
  }
  else{
    couple <- data.frame(Y=Y[temp_order],X=X[temp_order],Z=Z[temp_order])
  }

  unique_couple <- unique(couple)
  index_list <- find_unique(couple,unique_couple)

  w <- rep(NA, length(index_list))
  for (i in 1:length(index_list)){
    if (length(index_list[[i]])==1){
      w[index_list[[i]]] <- 1
    }
    else{
      w[index_list[[i]]] <- rep(1/length(index_list[[i]]),length(index_list[[i]]))
    }
  }
  return(w)
}
