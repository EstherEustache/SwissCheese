#' @title calib2
#'
#' @description Compute imputation probabilities matrix satisfying calibration equations.
#' Only the \code{k}-nearest respondent units of each nonrespondent have non-zero imputation probabilities.
#'
#'
#'
#' @param Xr a matrix without NA values. The rows correspond to respondent units.
#' @param Xm a matrix with at least one NA values on each of its rows. The rows correspond to nonrespondent units.
#' @param dr a vector containing the sampling weights of the respondent units in \code{Xr}. If NULL (default), all sampling weights are equal to 1.
#' @param dm a vector containing the sampling weights of the nonrespondent units in \code{Xm}. If NULL (default), all sampling weights are equal to 1.
#' @param knn a matrix that is returned by the function \code{\link{indKnn}}.
#' @param tol a tolerance parameter. Default value is 1e-3.
#' @param max_iter the maximum number of iterations to consider convergence.
#'
#'
#'
#'
#' @details
#' The main idea is to have imputation probabilities satisfying some calibration constraints for all variables simultaneously
#' and summing to one for each nonrespondent. The method is based on four requirements:
#' \enumerate
#'       \item The imputed values is selected among the completely observed units: a donor imputation method is used.
#'       \item Only one donor is selected per nonrespondent: all missing values of a unit should be imputed by the same donor.
#'       \item The donors is selected among the \code{k} nearest neighbors of units with missing values.
#'       \item If the observed values of the nonrespondents were imputed, the total estimator of each variable will remain unchanged.
#' \enumerate
#' See the package's vignette to a complete description of the calibration constraints.
#'
#'
#' @return a matrix with the same length as \code{knn} containing imputation probabilities for each neighbors.
#'
#'
#'
#' @author Esther Eustache \email{esther.eustache@@unine.ch}
#'
#'
#'
#' @seealso \link{indKnn}
#'
#' @examples #A faire
#'
#' @export

calibrateKnn <- function(Xr, Xm, dr = NULL, dm = NULL, knn, tol = 1e-2, max_iter = 2000)
{
  ##------------------
  ##  Initialization
  ##------------------

  nr <- nrow(Xr)
  nm <- nrow(Xm)
  R  <- !is.na(Xm)
  Xm[is.na(Xm)] <- 0
  J  <- ncol(X)

  #------sampling weights
  if (is.null(dr)) { dr <- rep(1, nrow(Xr)) }
  if (is.null(dm)) { dm <- rep(1, nrow(Xm)) }

  #------Sums of variables of nonrespondents
  tot_m   <- colSums(dm * Xm * R)

  #------k
  k     <- ncol(knn)
  knn_k <- knn

  ##--------------------------------------
  ##  Matrix of imputation probabilities
  ##--------------------------------------

  psi   <- matrix(rep(1/k, nm*k), c(nm,k))
  w_j   <- rep(0, nr)
  w_tmp <- rep(0, nr)

  #----Tolerance
  error  <- rep(5, J) #condition of the while loop
  it_2   <- 0

  ###----Main loop
  while (max(error) > tol){
    error_p <- error
    it      <- 1
    while(it < 10){
      #----Calibration and update of the psi weigths
      for (j in 1:J) {
        ##-----Only Sm units that responded to variable j
        R1_j  <- which(R[,j]==1)
        knn_j <- knn_k[R1_j,]
        psi_j <- psi[R1_j,]
        dm_j  <- dm[R1_j]

        for(i in 1:nr){
          w_j[i] <- sum((psi_j*dm_j)[knn_j == i])
        }

        ##-----update psi_j for Sm units that have responded to variable j
        lambda <- 0
        c      <- 0
        div    <- 1
        V      <- var(Xr[,j])

        while((div > 1e-10) & (c < 100)){
          c      <- c+1
          lambda <- lambda - (sum(w_j*Xr[,j]*exp(lambda*Xr[,j]))-tot_m[j])/sum(w_j*Xr[,j]*Xr[,j]*exp(lambda*Xr[,j]))
          div    <- (sum(w_j*Xr[,j]*exp(Xr[,j]*lambda))-tot_m[j])/V
        }

        if(c < 100){
          g <- exp(lambda*Xr[,j])
        } else {
          return(NULL)
        }

        g[is.nan(g)] <- rep(0, sum(is.nan(g)))

        for (i in 1:nr) {
          psi_j[knn_j == i] <- psi_j[knn_j == i]*g[i]
        }
        psi[R1_j,] <- psi_j
        psi        <- t(apply(psi, 1, function(x) x/sum(x)))  #----normalization
      }
      it <- it+1
    }

    #----Normalization
    psi <- t(apply(psi, 1, function(x) x/sum(x)))

    for(j in 1:J){
      #----Consider only respondents of variable j among Sm
      R1_j     <- which(R[,j]==1)
      knn_tmp  <- knn_k[R1_j,]
      psi_tmp  <- psi[R1_j,]
      dm_tmp   <- dm[R1_j]

      for(i in 1:nr){
        w_tmp[i] <- sum((psi_tmp*dm_tmp)[knn_tmp == i])
      }
      error[j] <- abs((Xr[,j] %*% w_tmp - tot_m[j]) / tot_m[j])
    }
    it_2 <- it_2+1
    if(((max(error_p)-max(error))<1e-8) | (it_2 == max_iter)){
      return(NULL)
    }
  }

  return(psi)
}
