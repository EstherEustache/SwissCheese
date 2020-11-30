#' @title Donor imputation using the cube method
#'
#' @description Impute missing values with one of the nearest respondents using the cube method.
#'
#'
#'
#' @param Xr a matrix without NA values. The rows correspond to respondent units.
#' @param Xm a matrix with at least one NA values on each of its rows. The rows correspond to nonrespondent units.
#' @param knn a matrix containing in rows the nonrespondent units and in columns the decreasing rank
#' of the \code{K}-Nearest Neighbour among the respondents units.
#' @param psi a matrix with the same dimension of \code{KNN} contains imputation probabilities of the
#' \code{K}-Nearest Neighbor of the nonrespondent units.
#'
#'
#'
#'
#'
#' @references
#' Deville, J. C. and Tille, Y. (2004). Efficient balanced sampling: the cube method. Biometrika, 91(4), 893-912.
#' \emph{Biometrika}, 91(4), 893-912
#'
#' Chauvet, G. and Tille, Y. (2006). A fast algorithm for balanced sampling.
#' \emph{Computational Statistics}, 21(1), 53-62.
#'
#'
#'
#'
#'
#' @return the dataframe of nonrespondents with NA values imputed.
#'
#'
#' @author Esther Eustache, \email{esther.eustache@@unine.ch}
#'
#'
#'
#' @seealso \link{calibrateKnn} and \link[BalancedSampling:flightphase]{BalancedSampling::flightphase}
#'
#' @export
#'
#' @examples #A faire
#'
cubeImput <- function(Xr, Xm, knn, psi){
      if ((nrow(psi)!=nrow(knn))&(ncol(psi)!=ncol(knn))) { stop('PSI and knn matrix must have the same dimension.') }

      ##------------------
      ##  Initialization
      ##------------------
      nm  <- nrow(Xm)
      k   <- nrow(knn)
      J   <- ncol(Xr)
      EPS <- 1e-7
      R   <- !is.na(Xm)
      Xm[!R] <- 0

      ##---------------
      ##  Cube method
      ##---------------
      X_pik      <- as.vector(psi)
      X_strat    <- Xr[as.vector(knn),]
      num_strat  <- rep(1:nm, each=k)
      XX         <- matrix(X_pik, ncol=J, nrow=nm*k, byrow=FALSE) * R[num_strat,] * X_strat

      pik_cube   <- balancedStratification(X = XX, pik = X_pik, stratum = num_strat, EPS)


      ##------------------
      ##  Imputed matrix
      ##------------------
      psi_cube <- matrix(pik_cube, nrow=k, ncol=nm, byrow=FALSE)
      Xm_imput <- R*Xm + (1-R)*Xr[knn[psi_cube > (1-EPS)],]

      return(Xm_imput)
}
