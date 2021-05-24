#' Performs one iteration of DFP for the Horseshoe Lasso
#'
#' This function performs one step of the DFP with the Horseshoe Lasso prior.
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
#' @param hl    Point estimate for the local shrinking parameter.
#' @param ht    Point estimate for the global shrinking parameters.
#' @param hv    Point estimate for the auxilary local shrinking hyper-parameter.
#' @param he    Point estimate for the auxilary global shrinking
#' hyper-parameters.
#' @param nmcmc Number of samples of the parameters.
#'
#' @return A list containing samples for the model parameters.
#' \describe{
#'   \item{sb}{A matrix with samples for coefficients, 1 sample per row.}
#'   \item{sl}{A matrix with samples for local shrinkage parameters, 1 sample
#'    per row.}
#'   \item{sv}{A matrix with samples for auxilary local shrinkage
#'   hyper-parameters, 1 sample per row.}
#'   \item{ss}{Samples for sigma.}
#'   \item{st}{Samples for the global shrinkage parameter.}
#'   \item{se}{Samples for the auxilary global shrinkage hyper-parameter.}
#' }
#'
#' @author Rene Gutierrez Marquez

#' @export

dfpHorShoeSte <- function(P,
                         XX,
                         Xy,
                         yy,
                         sN,
                         hb,
                         hs,
                         hl,
                         ht,
                         hv,
                         he,
                         nmcmc){

  # Gets the number of coefficients
  p <- length(hb)

  # Batch Samples
  ## Sample for beta
  sb <- matrix(0, nmcmc, p)
  ## Sample for lambda
  sl <- matrix(0, nmcmc, p)
  ## Sample for v
  sv <- matrix(0, nmcmc, p)
  ## Sample for sigma
  ss <- vector('numeric', nmcmc)
  ## Sample for tau2
  st <- vector('numeric', nmcmc)
  ## Sample for xi
  se <- vector('numeric', nmcmc)

  # Samples Beta and Lambda
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
  l <- hl
  ## Samples in Parallel
  aux <- foreach::foreach(j         = 1:parDiv,
                          .combine  = cbind,
                          .packages = c("statmod","mvtnorm")) %dopar% {
    # Component to work with
    ind    <- P[[j]]
    bb     <- b[ind]
    ll     <- l[ind]
    # Number of elements in the partition
    parCoe <- length(ind)
    for(k in 1:nmcmc){
      shal   = 1                            # Computes the Shape Parameter
      scal   = 1/hv[ind]+bb^2/(2*ht*hs)     # Computes the Rate Parameter
      ll     = 1/rgamma(parCoe, shal, scal) # Samples lambda
      ## Saves the Sample
      sl[k,ind] <- ll

      # Samples Beta
      if( parCoe == 1 )
      {
        indVar <- solve(XX[ind, ind] + 1/ll/ht) * hs
      } else {
        indVar <- solve(XX[ind, ind] + diag(1/ll)/ht) * hs
      }
      # Mean of the component beta
      indMea <- indVar %*% (Xy[ind] - XX[ind,-ind] %*% cb[-ind]) / hs
      # Samples the component beta
      bb     <- mvtnorm::rmvnorm(1, indMea, indVar)
      # Samples the component beta
      sb[k,ind] <- bb
    }
    # Output Parallel Computation
    rbind(matrix(sb[,ind],nmcmc,parCoe),matrix(st[,ind],nmcmc,parCoe))
  }
  ## Saves the Samples
  ### Samples the component beta
  sb[,unlist(P)] <- aux[1:nmcmc,]
  ### Saves Samples for Lambda2
  sl[,unlist(P)] <- aux[(nmcmc+1):(2*nmcmc),]
  ### Point Estimate for Lambda2
  hl <- colMeans(sl)
  ### Point Estimate for beta
  hb <- colMeans(sb)

  ### Samples nu
  sv = 1/rgamma(sim*p,1,rep(1+1/hl,sim))              # Samples nu from an Inverse Gamma
  sv = matrix(sv, nrow = sim, ncol = p, byrow = TRUE) # Samples to Matrix Form
  # Point Estimate of nu
  hv <- mean(sv)

  ### Samples Sigma, tau and xi
  for(a in 1:sim)
  {
    ### Samples Sigma2
    shas = (n*i+p)/2                                                                        # Sigma2 Shape
    scas = matrix(( yy - 2 * t(hb) %*% Xy + t(hb) %*% XX %*% hb  + sum(hb * hb / hl)/t )/2) # Sigma2 Scale
    s    = 1/rgamma(1, shas, scas)                                                          # Sigma2 Sample from an Inverse Gamma

    ### Samples tau
    shat = (p+1)/2                # tau shape
    scat = 1/e + sum(hb^2/hl)/(2*s) # tau scale
    t    = 1/rgamma(1,shat,scat)  # Samples tau from an Inverse Gamma

    ### Samples xi
    e = 1/rgamma(1,1,1+t) # Samples xi from an Inverse Gamma

    ### Saves the Samples
    ss[a]  = s  # Saves the sigma sample
    st[a]  = a  # Saves the tau sample
    se[a]  = e  # Saves the xi sample
  }
  ### Point Estimates
  hs = mean(ss) # sigma2 Point Estimate
  ht = mean(st) # tau Point Estimate
  he = mean(se) # Point Estimate xi

  return(list(sb = sb, st = st, ss = ss, sl = sl, hv, he))
}
