#' @title Donor imputation using the cube method
#'
#' @description Impute missing values with one of the nearest respondents using the cube method.
#'
#'
#'
#' @param X a dataframe with some NA values.
#'
#' @param KNN a matrix containing in rows the nonrespondent units and in columns the decreasing rank
#' of the \code{K}-Nearest Neighbour among the respondents units.
#'
#' @param PSI a matrix with the same dimension of \code{KNN} contains imputation probabilities of the
#' \code{K}-Nearest Neighbor of the nonrespondent units.
#'
#'
#'
#' @details (cube method)
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
#' @return Returns a list including:
#'
#' @return \code{X_imput} the new dataframe with NA values imputed.
#'
#' @return \code{PSI_new} the imputation matrix with the same length of \code{KNN} which contains imputation probabilities equal to 0 or 1.
#'
#'
#'
#' @author Esther Eustache, \email{esther.eustache@@unine.ch}
#'
#'
#'
#' @seealso \link{calibrate.knn} and \link[BalancedSampling:flightphase]{BalancedSampling::flightphase}
#'
#' @export
#'
#' @examples
#' df         <- matrix(c(1,2,1,3,2,
#'                        NA,3,4,5,2,
#'                        1,3,NA,NA,3,
#'                        3,2,5,7,8,
#'                        5,3,4,2,2,
#'                        6,NA,2,2,1,
#'                        4,2,4,4,3,
#'                        2,5,6,3,7,
#'                        1,4,2,3,1,
#'                        NA,NA,NA,NA,2), nrow = 10, ncol = 5, byrow = TRUE)
#' res        <- indKnn(X = df,  k=5)
#' KNN        <- res$KNN
#' res1       <- calibrateKnn(X = df, d = NULL, KNN = KNN, PSI.init = NULL, k = NULL, max.iter = 1000, tolerance = 0.01)
#' k_opt      <- res3$k_opt
#' PSI        <- res3$PSI
#' res_imput  <- cubeImput(X = df, PSI = PSI, KNN = KNN[,1:k.opt])
#' df_new     <- res_imput$X_imput # df with NA values imputed
#' df_new
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
