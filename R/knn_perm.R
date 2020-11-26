knn_perm <- function(X,Y,Z,k=5,B){
  N_i <- matrix(NA,length(Z),k)
  distance <- rep(NA,length(Z))
  res_init <- CDF(X,Y,Z,method="logistic regression")
  statistic_init <- max(abs(res_init$ccdf_nox-res_init$ccdf_x))
  statistic_b <- rep(NA,B)

  # knn
  for (z in 1:length(Z)){
    distance <- abs(Z-Z[z])
    N_i[z,] <- sort(distance,index.return=TRUE)$ix[1:k]
  }

  # permutations
  for (b in 1:B){
    U <- c()
    x_star <- c()
    N_i <- N_i[,sample(ncol(N_i))]
    rand_perm <- sample(1:length(Z))
    for (i in 1:length(Z)){
      j <- N_i[rand_perm[i],1] # N_i(0)
      m <- 0
      while ((j %in% U) & (m<(k-1))){
        m <- m+1
        j <- N_i[rand_perm[i],m]
      }
      x_star[rand_perm[i]] <- X[j]
      U <- c(U,j)
    }
    res <- CDF(x_star,Y,Z,method="logistic regression")
    statistic_b[b] <- max(abs(res$ccdf_nox-res$ccdf_x))
  }
  pval <- (1/(B+1))*(sum(statistic_b>statistic_init)+1)
  return(X_star)
}

###

knn_perm <- function(X,Y,Z,k=5){
  N_i <- matrix(NA,length(Z),k)
  distance <- rep(NA,length(Z))

  # knn
  for (z in 1:length(Z)){
    distance <- abs(Z-Z[z])
    N_i[z,] <- sort(distance,index.return=TRUE)$ix[1:k]
  }

  U <- c()
  x_star <- c()
  N_i <- N_i[,sample(ncol(N_i))]
  rand_perm <- sample(1:length(Z))
  for (i in 1:length(Z)){
    j <- N_i[rand_perm[i],1] # N_i(0)
    m <- 0
    while ((j %in% U) & (m<(k-1))){
      m <- m+1
      j <- N_i[rand_perm[i],m]
    }
    x_star[rand_perm[i]] <- X[j]
    U <- c(U,j)
  }

  return(X_star)
}

