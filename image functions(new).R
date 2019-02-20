library(RSpectra)
# NJW #
image.njw <- function(data,b,sigma_p,x) {
  ## produce W ##
  W <- matrix(0,ncol=nrow(data)*ncol(data),nrow = nrow(data)*ncol(data))
  #### do by row #####
  ## small box b*b (assume b is odd) ###
  sigma_l <- (b-1)/2
  for (i in 1:nrow(data)) {
    for (j in 1:ncol(data)) {
      # put similities values in small boxes of (i,j) to the W 
      for (k in (i-(b-1)/2):(i+(b-1)/2)) {
        if (k>0 & k<=nrow(data)) {
          for (l in (j-(b-1)/2):(j+(b-1)/2)) {
            if (l>0 & l<=ncol(data)) {
              d1 <- max(abs(i-k),abs(j-l))^2/(2*sigma_l^2)
              d2 <- (data[k,l]-data[i,j])^2/(2*sigma_p^2) 
              W[nrow(data)*(l-1)+k,nrow(data)*(j-1)+i] <- exp(-d1-d2)
            }
          }
        }
      }
    }
  }
  #########
  ## cluster the image by W ##
  A <- W
  D <- 1/sqrt(rowSums(A))
  A1 <- t(t(A*D)*D)
  U <- eigs(A1, k = x)$vectors
  # Kmeans on U
  labelsout <- kmeans(U, centers = x)$cluster
  cluster <- as.integer(as.factor(labelsout))
  result <- matrix(cluster,ncol = ncol(data))   ## convert the pixels from a vector to the image ###
  plot((as.cimg((result))))  ## plot result ##
}

#landmarks##
## LSC ##
#############################
### B size of big box B*B ###
image.lsc <- function(data,B,b,sigma_p,x) {
  landmarks_col = ceiling(ncol(data)/b) - 1
  landmarks_row = ceiling(nrow(data)/b) - 1
  n_l = landmarks_col*landmarks_row
  d_col = b
  d_row = b
  ### compute the similarity matrix W ##
  W <- matrix(0,nrow = nrow(data)*ncol(data),ncol = n_l)
  #### do by colw #####
  small_l <- (b-1)/2
  big_l <- (B-1)/2          
  sigma_l <- (B-1)/2
  for (i in 1:landmarks_row) {
    for (j in 1:landmarks_col) {
      n <- 0
      total <- 0
      for (k in (i*d_row-(b-1)/2):(i*d_row+(b-1)/2)) {
        if (k>0 & k<=nrow(data)) {
          for (l in (j*d_col-(b-1)/2):(j*d_col+(b-1)/2)) {
            if (l>0 & l<=ncol(data)) {
              n <- n+1
              total <- total+data[k,l]
            }
          }
        }
      }
      ave_small <- total/n     # average in small box
      for (k in (i*d_row-(B-1)/2):(i*d_row+(B-1)/2)) {
        if (k>0 & k<=nrow(data)) {
          for (l in (j*d_col-(B-1)/2):(j*d_col+(B-1)/2)) {
            if (l>0 & l<=ncol(data)) {
              d1 <- max(abs(i*d_row-k),abs(j*d_col-l))^2/2/sigma_l^2
              if (sigma_p==0) d2 <- 0
              else d2 <- (data[k,l]-ave_small)^2/2/sigma_p^2 
              W[nrow(data)*(l-1)+k,landmarks_row*(j-1)+i] <- exp(-d1-d2)
            }
          }
        } 
      }
    }
  }
  ###################################
  ## cluster the pixels by W ###
  A <- W
  A[(rowSums(A) == 0),] <- .Machine$double.eps
  A1  <- A / rowSums(A)  # Calculate A1
  A2 <- t(t(A1) / sqrt(rowSums(t(A1))))  # Calculate A2
  U <- svds(A2, k = x,nu = x, nv = 0)$u
  labelsout <- kmeans(U, centers = x)$cluster
  cluster <- as.integer(as.factor(labelsout))
  result <- matrix(cluster,ncol = ncol(data))   # convert the pixels as vextor to an image matrix
  plot((as.cimg((result))))  # plot the results
}

plot((as.cimg(data)))


