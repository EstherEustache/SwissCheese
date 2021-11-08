#' @title Balanced donor imputation method handling Swiss cheese nonresponse
#'
#' @description It imputes missing values by using a donor imputation method.
#' It extends the balanced \code{k}-nearest neighbor imputation (Hasler and Tille, 2016) to the treatment of the Swiss cheese nonresponse.
#' A linear programming is used to compute the imputation probability matrix while satisfying some constraints.
#'
#'
#' @param X a matrix with NA values. The rows correspond to the units.
#' @param d a vector containing the sampling weights of the units. If NULL (default), all sampling weights are equal to 1.
#' @param w a vector containing the weights of the variables.
#' Each variable is multiplied by its weight in order to give more (if the weight is greater than 1) or less importance (if the weight is less than 1) when calculating the distance between the units.
#' If NULL (default), all weights are equal to 1.
#' @param k_min minimum number of neighbors with a non-zero imputation probability for each nonrespondent.
#' @param k_max maximum number of neighbors with a non-zero imputation probability for each nonrespondent.
#' @param rand if TRUE, the imputation will be random. If FALSE, the imputation will be determinist.
#' @param test if TRUE, the function checks if the constraints are met.
#'
#' @return the imputed matrix of \code{X}.
#'
#' @author Esther Eustache \email{esther.eustache@@unine.ch}
#'
#'
#'
#' @seealso \code{\link[StratifiedSampling:stratifiedcube]{StratifiedSampling::stratifiedcube}},
#' \code{\link[Rglpk:Rglpk_solve_LP]{Rglpk::Rglpk_solve_LP}}, \code{\link[slam:simple_triplet_matrix]{slam::simple_triplet_matrix}}
#'
#' @examples
#' Xr  <- rbind(c(0.1,0.3,0.4,0.1), c(0.1,0.3,0.2,0.1), c(0.1,0.2,0.3,0.1),
#'              c(0.2,0.3,0.2,0.3), c(0.1,0.1,0.2,0.1))
#' Xm  <- rbind(c(NA,0.1,NA,0.1), c(0.1,NA,0.2,NA))
#' X   <- rbind(Xr,Xm)
#' lpSparseSwissCheese(X)
#'
#' @export
#'
lpSparseSwissCheese <- function(X, d = NULL, w = NULL, k_min = 5, k_max = NULL, rand = TRUE, test = FALSE)
{
  ##------------------
  ##  Initialization
  ##------------------
  eps <- 1e-3

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
  if(is.null(k_max)){ k_max <- nr }
  D   <- rep(0, k_max*nm)
  ind <- rep(0, k_max*nm)
  # D_all <- matrix(rep(0, nm*nr), nrow = nr, ncol = nm)
  for(i in 1:nm){
    d_tmp                          <- (sqrt(rowSums(t(t(Xr)*R[i,]*w-Xm[i,]*R[i,]*w)^2)))/sum(R[i,])
    ind[((i-1)*k_max+1):(i*k_max)] <- which(order(d_tmp) %in% (1:k_max))
    D[((i-1)*k_max+1):(i*k_max)]   <- d_tmp[which(order(d_tmp) %in% (1:k_max))]
    # D_all[,i]                      <- d_tmp
  }


  ##-----------------------------------
  ##  First condition for convergence
  ##-----------------------------------
  tot_m     <- colSums(Xm)

  for(j in 1:J){
    R1_j      <- which(R[,j] == 1)
    sum_max_j <- 0
    for(v in R1_j){
      ind_tmp   <- ind[((v-1)*k_max+1):(v*k_max)]
      sum_max_j <- sum_max_j+max(Xr[ind_tmp,j])
    }
    if(sum_max_j<tot_m[j]){ stop('k_max must be bigger.') }
  }


  ##----------------------------
  ##  Imputation probabilities
  ##----------------------------

  C_sparse <- NULL
  for(j in 1:ncol(Xm)){
    for(v in 1:nm){
      ind_tmp <- ind[((v-1)*k_max+1):(v*k_max)]
      C_v   <- rep(0,k_max)
      non_0 <- rep(0,k_max)
      for(u in 1:k_max){
        u_ind  <- ind_tmp[u]
        C_v[u] <- dm[v]*R[v,j]*Xr[u_ind,j]
      }
      non_0    <- abs(C_v)>1e-8
      if(any(non_0)){
        C_sparse <- rbind(C_sparse, cbind(rep(j, sum(non_0)), (((v-1)*k_max+1):(v*k_max))[non_0], C_v[non_0]))
      }
    }
  }
  B1 <- colSums(dm*R*Xm)

  ##---constraints of summing to 1
  move     <- max(C_sparse[,1])
  I_sparse <- cbind(rep(1:nm, each = k_max)+move, 1:(nm*k_max), rep(1,nm*k_max))

  B2 <- rep(1,nm)

  ##--- concatenation of the constraints
  A_tmp    <- rbind(C_sparse, I_sparse)
  A_sparse <- slam::simple_triplet_matrix(i = A_tmp[,1], j = A_tmp[,2], v = A_tmp[,3], nrow = max(A_tmp[,1]), ncol = max(A_tmp[,2]), dimnames = NULL)
  B        <- c(B1, B2)

  ##--The linear program
  dir_sparse <- rep("==", length(B1)+length(B2))

  if(!is.null(k_min)){
    p      <- 1/k_min
    bounds <- list(upper = list(ind = 1:max(A_tmp[,2]), val = rep(p,max(A_tmp[,2]))))

    lp_res <- Rglpk::Rglpk_solve_LP(obj = D,
                                    mat = A_sparse,
                                    dir = dir_sparse,
                                    rhs = B,
                                    max = FALSE,
                                    bounds)$solution
  }else{
    lp_res <- Rglpk::Rglpk_solve_LP(obj = D,
                                    mat = A_sparse,
                                    dir = dir_sparse,
                                    rhs = B,
                                    max = FALSE)$solution
  }


  if(all(lp_res < 1e-6)){ stop('No feasible solution found for the linear program.') }

  ##----------------------
  ##  Imputation matrix
  ##----------------------

  prob_imp     <- rep(0,nr*nm)
  for(v in 1:nm){
    prob_tmp                                 <- rep(0,nr)
    prob_tmp[ind[((v-1)*k_max+1):(v*k_max)]] <- lp_res[((v-1)*k_max+1):(v*k_max)]
    prob_imp[((v-1)*nr+1):(v*nr)]            <- prob_tmp
  }




  # ###############################################################################
  if(test){
    ##--------CHECK CONSTRAINTS--------##

    ##---Constraints of calibration
    C  <- matrix(rep(0, nr*nm*ncol(Xm)), nrow = ncol(Xm))
    for(j in 1:(ncol(Xm))){
      C1 <- matrix(rep(0, nm*nr), ncol = nm, nrow = nr)
      for(v in 1:nm){
        ind_tmp <- ind[((v-1)*k_max+1):(v*k_max)]
        for(u in 1:k_max){
          u_ind   <- ind_tmp[u]
          C1[u_ind,v] <- dm[v]*R[v,j]*Xr[u_ind,j]
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

    if(all(((A%*%prob_imp)-B)<1e-6)){ cat('Constraints of the program are met') }
  }
  ####################################################################################

  Xm_init      <- as.matrix(X_init[Sm,])
  Xm_init[!Rr] <- 0
  Xr_init      <- as.matrix(X_init[Sr,])

  if(rand){
    M_imput <- matrix(StratifiedSampling::stratifiedcube(X = as.matrix(rep(0,length(prob_imp))), strata = rep(1:nm,each = nr), pik = prob_imp), byrow = FALSE, ncol = nm)

    for(i in 1:nm){
      Xm_init[i,] <- Rr[i,]*Xm_init[i,] + (1-Rr[i,])*Xr_init[M_imput[,i] > (1-1e-6),]
    }
  }else{
    M_prob   <- matrix(prob_imp, byrow = FALSE, nrow = nr)

    for(i in 1:nm){
      Xm_init[i,] <- Rr[i,]*Xm_init[i,] + (1-Rr[i,])*colSums(M_prob[,i]*Xr_init[,])
    }
  }

  return(X_imp = X_init)
}


