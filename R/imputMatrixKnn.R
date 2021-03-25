#' @title Compute an imputation matrix
#'
#' @description Compute an imputation matrix using the cube method and a matrix of imputation probabilities.
#' Only the \code{k}-nearest respondent of each nonrespondent have non-zero imputation probabilities.
#'
#'
#'
#' @param Xr a matrix without NA values. The rows correspond to the responding units.
#' @param psi a matrix of imputation probabilities  that is returned by the function \code{\link{calibrateKnn}}.
#' @param knn a matrix that is returned by the function \code{\link{indKnn}}.
#'
#'
#'
#' @details
#' The function \code{\link[StratifiedSampling::stratifiedcube]{StratifiedSampling::stratifiedcube}} from the package \code{StratifiedSampling} is used to compute the imputation matrix to choose the donors.
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
#' @return The imputation matrix with the same size as matrix \code{\link{knn}} composed of 0s and 1s specifying which respondent will imput each nonrespondent.
#'
#'
#'
#' @author Esther Eustache \email{esther.eustache@@unine.ch}
#'
#'
#'
#' @seealso \code{\link[StratifiedSampling::stratifiedcube]{StratifiedSampling::stratifiedcube}}, \code{\link{calibrateKnn}}
#'
#'
#' @examples
#' Xr  <- rbind(c(0.1,0.3,0.4,0.1), c(0.1,0.3,0.2,0.1), c(0.1,0.2,0.3,0.1),
#'              c(0.2,0.3,0.2,0.3), c(0.1,0.1,0.2,0.1), c(0.3,0.1,0.1,0.1))
#' Xm  <- rbind(c(NA,0.1,NA,0.1), c(0.1,NA,0.2,NA))
#' X   <- rbind(Xr,Xm)
#' knn <- indKnn(Xr,Xm)[,1:5]
#' psi <- calibrateKnn(Xr, Xm, knn = knn)
#' imputMatrixKnn(X, psi, knn)
#'
#' @export
#'
imputMatrixKnn <- function(X, psi, knn)
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

  Sr <- which(apply(r, 1, function(x) !any(0 %in% x)))  # Sr: responding units
  Sm <- which(apply(r, 1, function(x) any(0 %in% x)))   # Sm: nonresponding units
  nr <- length(Sr)                                      # nr: number of respondents
  nm <- length(Sm)                                      # nm: number of nonrespondents

  R  <- r[Sm,]                      # r: matrix of responds among the nonrespondents

  Xr <- as.matrix(X[Sr,])
  Xm <- as.matrix(X[Sm,])


  #-------Variable with at least 2 respondents
  TEST <- which(colSums(R) != 0)
  Xr   <- Xr[,TEST]
  Xm   <- Xm[,TEST]
  R    <- R[,TEST]
  J    <- ncol(X)

  k    <- ncol(knn)


  ##--------------------
  ##  Imputation matrix
  ##--------------------

  #-------Cube method
  Xm[!R]       <- 0
  X_psi        <- as.vector(t(psi))
  X_strat      <- Xr[as.vector(t(knn)),]
  num_strat    <- rep(1:nm, each=k)
  XX           <- matrix(X_psi, ncol=J, nrow=k*nm, byrow=FALSE) * X_strat * R[num_strat,]
  psi_cube     <- StratifiedSampling::stratifiedcube(cbind(X_psi,XX),as.matrix(num_strat),X_psi)
  imput_matrix <- matrix(psi_cube, byrow = TRUE, nrow = nrow(knn))

  return(imput_matrix)
}
