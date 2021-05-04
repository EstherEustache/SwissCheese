#' @title Linear imputation of Swiss cheese nonresponse
#'
#' @description Impute missing values by means of a deterministic approach using a matrix of imputation probabilities.
#'
#'
#' @param X a matrix with NA values. The rows correspond to the units.
#' @param d a vector containing the sampling weights of the units. If NULL (default), all sampling weights are equal to 1.
#' @param k the number of neighbors considered to impute each nonrespondent. If NULL (default), the smallest \code{k}
#' allowing a solution to the calibration equations will be chosen.
#' @param tol a tolerance parameter. Default value is 1e-3.
#' @param max_iter the maximum number of iterations to consider convergence.
#'
#'
#' @return Returns a list including:
#'
#' @return \code{X_new} the imputed matrix of \code{X}.
#' @return \code{k} the number of neighbors taken into account.
#'
#'
#' @details
#' First, the \code{k} nearest neighbors of each nonrespondent are chosen using the Euclidean distance in the function \code{\link{indKnn}}.
#' Next, imputation probabilities are computed for the nearest neighbors of each nonrespondent.
#' The imputation probabilities satisfy some calibration constraints for all variables simultaneously.
#' The matrix of imputation probabilities is computed using the function \code{\link{calibrateKnn}}.
#' Then, each missing value is imputed by the sum of the multiplication of the imputation probabilities and the corresponding observed values.
#'
#'
#' @seealso \code{\link{indKnn}} and \code{\link{calibrateKnn}}
#'
#' @examples
#' Xr  <- rbind(c(0.1,0.3,0.4,0.1), c(0.1,0.3,0.2,0.1), c(0.1,0.2,0.3,0.1),
#'              c(0.2,0.3,0.2,0.3), c(0.1,0.1,0.2,0.1))
#' Xm  <- rbind(c(NA,0.1,NA,0.1), c(0.1,NA,0.2,NA))
#' X   <- rbind(Xr,Xm)
#' deterministicImput(X)
#'
#' @export
#'
deterministicImput <- function(X, d = NULL, k = NULL, tol = 1e-3, max_iter = 2000){
  ##------------------
  ##  Initialization
  ##------------------

  #--- Standardized data
  X_init <- X
  for(i in 1:ncol(X)){
    X[!is.na(X[,i]),i] <- as.vector(scale(X[!is.na(X[,i]),i]))
  }

  N  <- nrow(X)
  J  <- ncol(X)
  S  <- 1:nrow(X)
  r  <- (!is.na(X))*1                                   # r: matrix of responds

  Sr <- which(apply(r, 1, function(x) !any(0 %in% x)))  # Sr: responding units
  Sm <- which(apply(r, 1, function(x) any(0 %in% x)))   # Sm: nonresponding units
  nr <- length(Sr)                                      # nr: number of respondents
  nm <- length(Sm)                                      # nm: number of nonrespondents

  R  <- r[Sm,]                      # r: matrix of responds among the nonrespondents

  Xr <- X[Sr,]
  Xm <- X[Sm,]

  dr <- d[Sr,]
  dm <- d[Sm,]

  #-------Variable with at least 2 respondents
  TEST <- which(colSums(R) != 0)
  Xr   <- Xr[,TEST]
  Xm   <- Xm[,TEST]
  R    <- R[,TEST]
  J    <- ncol(X)

  #------Respondents ordering it term of distance with nonresponding units in rows
  knn   <- indKnn(Xr, Xm)

  #------k
  k_init <- k
  if(is.null(k)) { k  <- 2 }
  knn_k <- knn[,1:k]


  ##----------
  ##  Method
  ##----------

  #-----Imputation probabilities
  cvg <- 0
  while(cvg == 0){
    psi <- calibrateKnn(Xr, Xm, dr, dm, knn = knn_k, tol, max_iter)
    if(is.null(psi)){
      k     <- k+1
      if(k > nr){ stop('The algorithm does not converge. k is bigger than the number of respondents.') }
      knn_k <- knn[,1:k]
    } else {
      cvg <- 1
    }
  }

  #------Imputed matrix
  Xm_init     <- as.matrix(X_init[Sm,])
  Xm_init[!R] <- 0
  Xr_init     <- as.matrix(X_init[Sr,])

  for (i in 1:nm) {
    Xm_init[i,] <- R[i,]*Xm_init[i,] + (1-R[i,])*colSums(psi[i,]*Xr_init[knn_k[i,],])
  }
  X_init[Sm,] <- Xm_init

  return(list(X_new = X_init, k = k))
}




