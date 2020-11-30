#' @title Linear imputation of Swiss cheese nonresponse
#'
#' @description Impute missing values using the linear imputation.
#'
#'
#' @param X a matrix with NA values. Rows correspond to units.
#' @param d a vector containing the sampling weights of units. If NULL (default), all sampling weights are equal to 1.
#' @param k the number of neighbors considered to impute each nonrespondent. If NULL (default), the smaller k
#' that makes possible to satisfying calibration equations will be chosen.
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
#' @export
#'
#' @seealso \link{indKnn} and \link{calibrateKnn}
#'
#' @examples
#'
#'
linearImput <- function(X, d = NULL, k = NULL, tol = 1e-3, max_iter = 2000){
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

  Sr <- which(apply(r, 1, function(x) !any(0 %in% x)))  # Sr: respondent units
  Sm <- which(apply(r, 1, function(x) any(0 %in% x)))   # Sm: nonrespondent units
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

  #------Respondents ordering it term of distance with nonrespondent units in rows
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

  #if((!is.null(k_init)) && (k>k_init)){ warning('K must have been increased.') }

  return(list(X_new = X_init, k = k))
}




