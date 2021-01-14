#' @noRd
balancedStratification=function(X, pik, stratum, EPS=1e-6)
{
  # cleaning and random permutation of the strata
  stratum  <- sampling::cleanstrata(stratum)
  perm     <- sample(max(stratum),max(stratum))
  stratum1 <- stratum
  for (i in 1:max(stratum)) stratum1[stratum==i] <- perm[i]
  stratum  <- stratum1
  # balancing in each stratum
  pikt  <- pik
  XX    <- cbind(pik,X)
  for(i in 1:max(stratum))  pikt[stratum==i] <- BalancedSampling::flightphase(pikt[stratum==i], as.matrix(XX[stratum==i,]))
  # balancing between the strata
  A     <- (1/pik)*X
  for(i in 2:max(stratum))
  {
    T       <- (stratum %in% 1:i) & (pikt>EPS) & (pikt < 1-EPS)
    ss      <- sampling::cleanstrata(stratum[T])
    m       <- matrix(0, length(stratum[T]), length(unique(stratum[T])))
    for (i in 1:length(ss)) m[i, ss[i]] = 1
    At      <- cbind(m,A[T,])
    pikt[T] <- BalancedSampling::flightphase(pikt[T], as.matrix(pikt[T]*At))
  }
  for(i in 1:max(stratum)){
    if(sum(pikt[stratum == i]) < (1-1e-5) ){
      cat('\nProblem with the flight phase.\n')
      stop('break')
    }
  }
  # landing phase
  p  <- ncol(A)
  for(i in p:1){
    T       <- (pikt>EPS) & (pikt < 1-EPS)
    At      <- cbind(sampling::disjunctive(stratum[T]),A[T,1:i])
    pikt[T] <- BalancedSampling::flightphase(pikt[T], as.matrix(pikt[T]*At))
  }
  T       <-  (pikt>EPS) & (pikt < 1-EPS)
  At      <- cbind(sampling::disjunctive(stratum[T]))
  pikt[T] <- BalancedSampling::flightphase(pikt[T], as.matrix(pikt[T]*At))
  for(i in 1:max(stratum)){
    if(sum(pik[stratum == i]) < (1-1e-6) ){
      stop('\nProblem with the landing phase.\n')
      stop('break')
    }
  }
  # return the sample
  return(pikt)
}




