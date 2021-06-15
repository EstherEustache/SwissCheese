#' @title Donor imputation of Swiss cheese nonresponse using the cube method
#'
#' @description Impute missing values that have no particular pattern by using a donor imputation method.
#' It extends the balanced \code{k}-nearest neighbor imputation (Hasler and Tille, 2016) to the treatment of the Swiss cheese nonresponse.
#'
#'
#'
#' @param X a matrix with NA values. The rows correspond to the units.
#' @param d a vector containing the sampling weights of the units. If NULL (default), all sampling weights are equal to 1.
#' @param w
#' @param k
#'
#'
#' @return the imputed matrix of \code{X}.
#'
#'
#'
#' @author Esther Eustache \email{esther.eustache@@unine.ch}
#'
#'
#'
#' @seealso \code{\link[StratifiedSampling:stratifiedcube]{StratifiedSampling::stratifiedcube}}, \code{\link{indKnn}}, \code{\link{calibrateKnn}}
#'
#' @examples
#' Xr  <- rbind(c(0.1,0.3,0.4,0.1), c(0.1,0.3,0.2,0.1), c(0.1,0.2,0.3,0.1),
#'              c(0.2,0.3,0.2,0.3), c(0.1,0.1,0.2,0.1))
#' Xm  <- rbind(c(NA,0.1,NA,0.1), c(0.1,NA,0.2,NA))
#' X   <- rbind(Xr,Xm)
#' lpSwissCheeseImput(X)
#'
#' @export
#'
lpSwissCheeseImput <- function(X, d = NULL, w = NULL, k = 5)
{
  ##------------------
  ##  Initialization
  ##------------------
  eps <- 1e-4

  #--- Standardized data
  X_init <- X
  var_0 <- which(!apply(X, 2, function(x) sd(x[!is.na(x)]) == 0))
  X0 <- X[,var_0]
  for(i in 1:length(var_0)){
    X0[!is.na(X0[,i]),i] <- as.vector(scale(X0[!is.na(X0[,i]),i]))
  }
  X[,var_0] <- X0

  N  <- nrow(X)
  r  <- (!is.na(X))*1                                   # r: matrix of responds

  Sr <- which(apply(r, 1, function(x) !any(0 %in% x)))  # Sr: responding units
  Sm <- which(apply(r, 1, function(x) any(0 %in% x)))   # Sm: nonresponding units
  nr <- length(Sr)                                      # nr: number of respondents
  nm <- length(Sm)                                      # nm: number of nonrespondents

  Rr  <- r[Sm,]                      # r: matrix of responds among the nonrespondents
  X[is.na(X)] <- 0
  Xr <- as.matrix(X[Sr,])
  Xm <- as.matrix(X[Sm,])

  if(is.null(d)){ d <- rep(1, N) }
  dr <- d[Sr]
  dm <- d[Sm]

  #-------Variable with at least 2 respondents
  TEST <- which(colSums(Rr) != 0)
  Xr   <- Xr[,TEST]
  Xm   <- Xm[,TEST]
  R    <- Rr[,TEST]
  J    <- ncol(Xr)

  ##---Compute the Euclidean distance between respondents and non respondents
  if(is.null(w)){ w <- rep(1, J) }
  D <- rep(0, nr*nm)
  for(i in 1:nm){
    D[((i-1)*nr+1):(i*nr)] <- (sqrt(rowSums(t(t(Xr)*R[i,]*w-Xm[i,]*R[i,]*w)^2)))/sum(R[i,])
  }

  ##----------------------------
  ##  Imputation probabilities
  ##----------------------------

  ##---Constraints of calibration
  C  <- matrix(rep(0, nr*nm*ncol(Xm)), nrow = ncol(Xm))
  for(j in 1:(ncol(Xm))){
    C1 <- matrix(rep(0, nr*nm), ncol = nm, nrow = nr)
    for(v in 1:nm){
      for(u in 1:nr){
        C1[u,v] <- dm[v]*R[v,j]*Xr[u,j]
      }
    }
    C[j,] <- as.vector(C1)  # matrix of constraints
  }
  B1 <- colSums(dm*R*Xm) # matrix of the right hand side of the constraint equation

  ##---constraints of summing to 1
  I <- matrix(rep(0,nr*nm*nm), nrow = nm)
  for(v in 1:nm){
    I[v,((v-1)*nr+1):(v*nr)] <- 1
  }
  B2 <- rep(1,nm)

  ##--- concatenation of the constraints
  A <- rbind(C, I)
  B <- c(B1, B2)

  ##---Lpsolve
  dir <- rep("=", length(B1)+length(B2))


  lp_res <- lpSolve::lp (direction = "min",
                         objective.in = D,
                         const.mat = A,
                         const.dir = dir,
                         const.rhs = B,
                         transpose.constraints = TRUE)$solution

  if(!is.null(k)){
    p  <- 1/k
    II <- NULL

    while(any(lp_res>(p+eps))){
      ##---constraints of being lower than p
      larger <- which(lp_res>(p+eps)) #proba larger than p

      I_tmp  <- matrix(rep(0, length(larger)*nm*nr), nrow = length(larger), ncol = nm*nr)
      for(i in 1:length(larger)){
        I_tmp[i,larger[i]] <- 1
      }
      II <- rbind(II, I_tmp)

      B3 <- rep(p,nrow(II))

      ##--- concatenation of the constraints
      A <- rbind(C, I, II)
      B <- c(B1, B2, B3)

      ##---Lpsolve
      dir <- c(rep("=", length(B1)+length(B2)), rep("<", length(B3)))

      lp_res <- lpSolve::lp (direction = "min",
                             objective.in = D,
                             const.mat = A,
                             const.dir = dir,
                             const.rhs = B,
                             transpose.constraints = TRUE)$solution
    }
  }


  if(!any((apply(A, 1, function(x) sum(x*lp_res))-B)>1e-6)){
    cat('Calibration equation are exactly satisfied.')
  }else{
    cat('Calibration equation are not satisfied.')
  }

  ##----------------------
  ##  Imputation matrix
  ##----------------------
  Xm_init      <- as.matrix(X_init[Sm,])
  Xm_init[!Rr] <- 0
  Xr_init      <- as.matrix(X_init[Sr,])

  M_prob  <- matrix(lp_res, ncol = nm, byrow = F)

  M_X <- matrix(rep(0,nr*nm*nm), ncol = nm)
  for(v in 1:nm){
    M_X[((v-1)*nr+1):(v*nr),v] <- lp_res[((v-1)*nr+1):(v*nr)]
  }
  M_imput <- matrix(StratifiedSampling::stratifiedcube(X = as.matrix(rep(0,length(lp_res))), strata = rep(1:nm,each = nr), pik = lp_res), byrow = FALSE, ncol = nm)

  for(i in 1:nm){
    Xm_init[i,] <- Rr[i,]*Xm_init[i,] + (1-Rr[i,])*Xr_init[M_imput[,i] > (1-1e-6),]
  }
  X_init[Sm,] <- Xm_init

  return(X_init)
}


