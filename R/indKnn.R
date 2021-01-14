#' @title K-nearest neighbors of nonrespondent units
#'
#' @description Find the k-nearest neighbors of nonrespondent units among respondent units, in terms of Euclidean distance.
#' The function \code{\link[FNN:knn]{knn}} from the package \code{KNN} is used.
#'
#'
#' @param Xr a matrix without NA values. The rows correspond to respondent units.
#' @param Xm a matrix with at least one NA value on each of its rows. The rows correspond to nonrespondent units.
#'
#'
#' @return a matrix containing in columns the nonrespondent units and in rows the nearest respondent units for each nonrespondent in decreasing rank.
#'
#'
#'
#' @author Esther Eustache \email{esther.eustache@@unine.ch}
#'
#'
#'
#' @seealso \code{\link[FNN:knn]{FNN::knn}}
#'
#'
#'
#' @examples
#' Xr  <- rbind(c(0.1,0.3,0.4,0.1), c(0.1,0.3,0.2,0.1), c(0.1,0.2,0.3,0.1), c(0.2,0.3,0.2,0.3), c(0.1,0.1,0.2,0.1))
#' Xm  <- rbind(c(NA,0.1,NA,0.1), c(0.1,NA,0.2,NA))
#' knn <- indKnn(Xr,Xm)
#'
#'
#' @export
#'
indKnn <- function(Xr, Xm){
  ##------------------
  ##  Initialization
  ##------------------
  nr <- nrow(Xr)
  nm <- nrow(Xm)
  if (nm < 2) { stop('Xm must contains at least 2 units.') }
  if (nr < 2) { stop('Xr must contains at least 2 units.') }
  R <- !is.na(Xm)

  ##--------------------------
  ##  compute the neighbors
  ##--------------------------

  ##---Standardize matrices and replace NA by 0.
  X_std               <- apply(rbind(Xr, Xm), MARGIN=2, function(x) x/sqrt(1/length(stats::na.omit(x))*sum((stats::na.omit(x)-mean(stats::na.omit(x)))^2)))
  X_std[is.na(X_std)] <- 0
  Xr_std              <- X_std[1:nr,]
  Xm_std              <- X_std[(nr+1):(nm+nr),]

  ##---Compute the Euclidean distance between respondents and non respondents
  knn       <- matrix(rep(0, nm*nr), nrow = nm, ncol = nr)
  for (i in 1:nm){
    res     <- FNN::knn(Xr_std*R[i,], Xm_std[i,]*R[i,], c(1:nr), k = nr)
    knn[i,] <- attr(res, "nn.index")[1,]
  }

  return(knn)
}

