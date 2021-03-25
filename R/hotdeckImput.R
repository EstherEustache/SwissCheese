#' @title Random hot deck imputation of Swiss cheese nonresponse
#'
#' @description Impute each nonresponding unit by selecting randomly a donor among its nearest respondents in terms of Euclidean distance.
#'
#'
#'
#' @param X a matrix with NA values. The rows correspond to the units.
#' @param k the number of neighbors considered to impute each nonrespondent. If NULL (default), all the respondents will be considered.
#'
#'
#' @return The new dataframe with NA values imputed.
#'
#'
#'
#' @author Esther Eustache \email{esther.eustache@@unine.ch}
#'
#'
#' @details First, the \code{k} nearest neighbors of each nonrespondent are chosen in terms of Euclidean distance using the function \code{\link{indKnn}}.
#' Next, a donor is selected randomly among the neighbors of each nonrespondent in order to impute the missing values.
#'
#'
#' @seealso \code{\link{indKnn}}
#'
#' @references
#' Rebecca R. Andridge, Roderick J. A. Little (2010). A Review of Hot Deck Imputation for Survey Non-response. International Statistical Review, Vol. 78, Issue 1, pp. 40-64.
#'
#'
#' @examples
#' Xr  <- rbind(c(0.1,0.3,0.4,0.1), c(0.1,0.3,0.2,0.1), c(0.1,0.2,0.3,0.1),
#'              c(0.2,0.3,0.2,0.3), c(0.1,0.1,0.2,0.1))
#' Xm  <- rbind(c(NA,0.1,NA,0.1), c(0.1,NA,0.2,NA))
#' X <- rbind(Xr,Xm)
#' hotDeckImput(X)
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

  Sr <- which(apply(r, 1, function(x) !any(0 %in% x)))  # Sr: responding units
  Sm <- which(apply(r, 1, function(x) any(0 %in% x)))   # Sm: nonresponding units
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
