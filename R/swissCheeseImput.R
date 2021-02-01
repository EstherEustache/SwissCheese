#' @title Donor imputation of Swiss cheese nonresponse using the cube method
#'
#' @description Impute missing values that have not a particular pattern in the dataset by using a donor imputation.
#' It extends the balanced k-nearest neighbor imputation (Hasler and Tille, 2016) to the treatment of the Swiss cheese nonresponse.
#'
#'
#'
#' @param X a matrix with NA values. The rows correspond to the units.
#' @param d a vector containing the sampling weights of units. If NULL (default), all sampling weights are equal to 1.
#' @param k the number of neighbors considered to impute each nonrespondent. If NULL (default), the smaller \code{k}
#' that makes possible to satisfy calibration equations will be chosen.
#' @param tol a tolerance parameter. Default value is 1e-3.
#' @param max_iter the maximum number of iterations to consider convergence.
#'
#'
#'
#' @details
#' First, the \code{k} nearest neighbors of each nonrespondent are chosen in terms of Euclidean distance.
#' Next, imputation probabilities from nearest neighbors of each nonrespondent  are computed satisfying some calibration constraints for all variables simultaneously
#' and to sum to one for each nonrespondent. The method of calibration is based on four requirements (see function \code{\link{calibrateKnn}} and article ... on Arxiv
#' to a complete description of the calibration constraint).
#' Then, the function \code{\link[StratifiedSampling::stratifiedcube]{StratifiedSampling::stratifiedcube}} from the package \code{StratifiedSampling} is used to compute the final imputation matrix to choose the donors.
#' This function uses the cube method (see (Deville and Tille, 2004)).
#'
#'
#'
#' @references
#' Deville, J. C. and Tille, Y. (2004). Efficient balanced sampling: the cube method. Biometrika, 91(4), 893-912.
#' \emph{Biometrika}, 91(4), 893-912
#'
#'
#'
#' @return Returns a list including:
#' @return \code{X_new} the imputed matrix of \code{X}.
#' @return \code{k} the number of neighbors taken into account.
#'
#'
#'
#' @author Esther Eustache \email{esther.eustache@@unine.ch}
#'
#'
#'
#' @seealso \code{\link[StratifiedSampling::stratifiedcube]{StratifiedSampling::stratifiedcube}}, \code{\link{indKnn}}, \code{\link{calibrateKnn}}, \code{\link{linearImput}}
#'
#' @examples
#' Xr  <- rbind(c(0.1,0.3,0.4,0.1), c(0.1,0.3,0.2,0.1), c(0.1,0.2,0.3,0.1),
#'              c(0.2,0.3,0.2,0.3), c(0.1,0.1,0.2,0.1))
#' Xm  <- rbind(c(NA,0.1,NA,0.1), c(0.1,NA,0.2,NA))
#' X <- rbind(Xr,Xm)
#' swissCheeseImput(X)
#'
#' @export
#'
swissCheeseImput <- function(X, d = NULL, k = NULL, tol = 1e-3, max_iter = 2000)
{
  ##------------------
  ##  Initialization
  ##------------------

  #--- Standardized data
  X_init <- X
  for(i in 1:ncol(X)){
    X[!is.na(X[,i]),i] <- as.vector(scale(X[!is.na(X[,i]),i]))
  }

  N  <- nrow(X)
  S  <- 1:nrow(X)
  r  <- (!is.na(X))*1                                   # r: matrix of responds

  Sr <- which(apply(r, 1, function(x) !any(0 %in% x)))  # Sr: respondent units
  Sm <- which(apply(r, 1, function(x) any(0 %in% x)))   # Sm: nonrespondent units
  nr <- length(Sr)                                      # nr: number of respondents
  nm <- length(Sm)                                      # nm: number of nonrespondents

  R  <- r[Sm,]                      # r: matrix of responds among the nonrespondents

  Xr <- as.matrix(X[Sr,])
  Xm <- as.matrix(X[Sm,])

  if(is.null(d)){ d <- rep(1, N) }
  dr <- d[Sr]
  dm <- d[Sm]

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
  if (is.null(k)) { k <- 2 } else { k <- min(k,nr) }


  ##-----------------------------------
  ##  First condition for convergence
  ##-----------------------------------
  cvg       <- 0
  donor_max <- rep(0,J)
  tot_m     <- colSums(na.omit(Xm))

  while (cvg == 0) {
    #----Test if k is not to big
    if(k > nr){ stop('The algorithm does not converge. K is bigger than the number of respondents.') }

    knn_k <- knn[,1:k]

    #----Test of the condition
    for(j in 1:J){
      R0_j <- which(R[,j] == 0)
      # donor_j: contains all the possible donors for each variable
      if(sum(R[,j]) != nm){
        donor_j      <- unique(as.vector(knn_k[R0_j,]))
      } else {
        donor_j      <- unique(as.vector(knn_k))
      }
      # select the higgest value for each variable
      donor_max[j] <- max(dr[donor_j]*Xr[donor_j,j])
    }
    if (any(colSums(R)*donor_max < tot_m)) {
      cvg <- 0
      k   <- k+1
    }else{
      cvg <- 1
    }
  }

  ##----------------------------
  ##  Imputation probabilities
  ##----------------------------
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
    print(k)
  }

  #------Imputed matrix
  Xm_init     <- as.matrix(X_init[Sm,])
  Xm_init[!R] <- 0
  Xr_init     <- as.matrix(X_init[Sr,])

  psi_cube <- imputMatrixKnn(X_init, psi, knn_k)

  for(i in 1:nm){
    Xm_init[i,] <- R[i,]*Xm_init[i,] + (1-R[i,])*Xr_init[knn_k[i,][psi_cube[i,] > (1-1e-6)],]
  }

  X_init[Sm,] <- Xm_init

  return(list(X_new = X_init, k = k))
}
