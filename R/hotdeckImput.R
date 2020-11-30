#' @title Random hot deck imputation of Swiss cheese nonresponse
#'
#' @description Impute each nonrespondent unit by selecting randomly a donor among its nearest respondents in term of Euclidean distance
#' (see (Andridge and Little, 2004)).
#'
#'
#'
#' @param X a matrix with NA values. Rows correspond to units.
#' @param k the number of neighbors considered to impute each nonrespondent. If NULL (default), the minimum between
#' 30 will be chosen.
#'
#'
#' @return the new dataframe with NA values imputed.
#' @export
#'
#'
#'
#' @author Esther Eustache \email{esther.eustache@@unine.ch}
#'
#'
#'
#' @seealso \link{indKnn}
#'
#'
#'
#' @examples
#'
#' @export

hotDeckImput <- function(X, k=NULL){
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

  Xr     <- X[Sr,]
  Xm     <- X[Sm,]

  #-------Variable with at least 2 respondents
  TEST <- which(colSums(R) != 0)
  Xr   <- Xr[,TEST]
  Xm   <- Xm[,TEST]
  R    <- R[,TEST]
  J    <- ncol(X)


  ##----------
  ##  Method
  ##----------

  if(is.null(k)){
    hotdeck  <- sample(1:nr, nm)
  } else {
    knn      <- indKnn(Xr, Xm)
    hotdeck  <- apply(knn[,1:k], 1, function(x) sample(x,1))
  }

  #-----Imputed matrix
  Xm_init     <- as.matrix(X_init[Sm,])
  Xm_init[!R] <- 0
  Xr_init     <- as.matrix(X_init[Sr,])

  for (i in 1:nm) {
    Xm_init[i,] <- R[i,]*Xm_init[i,] + (1-R[i,])*Xr_init[hotdeck[i],]
  }
  X_init[Sm,] <- Xm_init

  return(X_init)
}
