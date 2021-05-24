#' Performs one iteration of DFP for the Bayesian Lasso
#'
#' This function performs one step of the DFP with the Bayesian Lasso prior.
#' The inputs are described in Algorithm 1 and the Supplementary Materials of
#' Guhaniyogi & Gutierrez.
#'
#' @param P     Current partition to be used.
#' @param XX    Sufficient statistic X'X.
#' @param Xy    Sufficient statistic X'y.
#' @param yy    Sufficient statistic y'y.
#' @param sN    Number of Observations.
#' @param hb    Point estimate of the coefficients.
#' @param hs    Point estimates of the variance of the error.
#' @param hl    Point estimate for the global shrinking parameter.
#' @param ht    Point estimate for the local shrinking parameters.
#' @param hlsh  Hyper-prior for the lambda parameter shape.
#' @param hlsc  Hyper-prior for the lambda parameter rate.
#' @param nmcmc Number of samples of the parameters.
#'
#' @return A list containing samples for the model parameters.
#' \describe{
#'   \item{sb}{A matrix with samples for coefficients, 1 sample per row.}
#'   \item{st}{A matrix with samples for local shrinkage parameters, 1 sample
#'    per row.}
#'   \item{ss}{Samples for sigma.}
#'   \item{sl}{Samples for the global shrinkage parameter.}
#' }
#'
#' @author Rene Gutierrez Marquez

#' @export

dfpBayLasSte <- function(P,
                         XX,
                         Xy,
                         yy,
                         sN,
                         hb,
                         hs,
                         hl,
                         ht,
                         hlsc,
                         hlsh,
                         nmcmc){

  # Gets the number of coefficients
  p <- length(hb)

  # Batch Samples
  ## Sample for beta
  sb <- matrix(0, nmcmc, p)
  ## Sample for tau
  st <- matrix(0, nmcmc, p)
  ## Sample for sigma
  ss <- vector('numeric', nmcmc)
  ## Sample for lambda
  sl <- vector('numeric', nmcmc)

  # Samples Beta and Tau
  ## Gets the Number Of Components in the Regression Coefficients
  parDiv <- length(P)
  ## Centers the Beta Values (to mimic coordinate ascent)
  cb <- hb
  for(j in 1:parDiv){
    # Component to work with
    ind     <- P[[j]]
    cb[ind] <- solve(XX[ind, ind], Xy[ind] - XX[ind, -ind] %*% cb[-ind])
  }
  ## Initial Values
  b <- cb
  t <- ht
  ## Samples in Parallel
  aux <- foreach::foreach(j         = 1:parDiv,
                          .combine  = cbind,
                          .packages = c("statmod","mvtnorm")) %dopar% {
    # Component to work with
    ind    <- P[[j]]
    # Number of elements in the partition
    parCoe <- length(ind)
    for(k in 1:nmcmc){
      # Samples Tau
      ## Mean of the component tau2
      tauMea    <- sqrt(hl * hs / b^2)
      ## Shape of component tau2
      tauSha    <- hl
      ## Samples
      t         <- 1 / statmod::rinvgauss(parCoe, tauMea, tauSha)
      ## Saves the Sample
      st[k,ind] <- t

      # Samples Beta
      if( parCoe == 1 )
      {
        indVar <- solve(XX[ind, ind] + 1/t) * hs
      } else {
        indVar <- solve(XX[ind, ind] + diag(1/t)) * hs
      }
      # Mean of the component beta
      # indMea <- indVar %*% (Xy[ind] - XX[ind,-ind] %*% hb[-ind]) / hs
      indMea <- indVar %*% (Xy[ind] - XX[ind,-ind] %*% cb[-ind]) / hs
      # Samples the component beta
      b      <- mvtnorm::rmvnorm(1, indMea, indVar)
      # Samples the component beta
      sb[k,ind] <- b
    }
    # Output Parallel Computation
    rbind(matrix(sb[,ind],nmcmc,parCoe),matrix(st[,ind],nmcmc,parCoe))
  }
  ## Saves the Samples
  ### Samples the component beta
  sb[,unlist(P)] <- aux[1:nmcmc,]
  ### Saves Samples for tau2
  st[,unlist(P)] <- aux[(nmcmc+1):(2*nmcmc),]
  ### Point Estimate for tau2
  ht <- colMeans(st)
  ### Point Estimate for beta
  hb <- colMeans(sb)

  # Samples Lambda
  ## Scale lambda2
  lamSca <- sum(ht)/2 + hlsc
  # Shape lambda2
  p      <- length(hb)
  lamSha <- hlsh + p
  # Samples lambda
  sl <-  1/rgamma(nmcmc, lamSha, lamSca)
  # Point Estimate of lambda2
  hl <- mean(sl)

  # Samples Sigma2
  ## Scale sigma2
  sigSca <- matrix(( yy - 2 * t(hb) %*% Xy + t(hb) %*% XX %*% hb  + sum(hb^2/ht) )/2)
  ## Shape sigma2
  p      <- length(hb)
  sigSha <- (sN + p - 1) / 2
  ## Samples sigma2
  ss <- 1/rgamma(nmcmc, sigSha, sigSca)
  hs <- mean(ss)

  return(list(sb = sb, st = st, ss = ss, sl = sl))
}
