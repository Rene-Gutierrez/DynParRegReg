#' Performs one iteration of DFP for the Bayesian Lasso
#'
#' This function performs one step of the DFP with the Bayesian Lasso prior.
#' The inputs are described in Algorithm 1 and the Supplementary Materials of
#' Guhaniyogi & Gutierrez.
#'
#' @param XX    Sufficient statistic X'X.
#' @param Xy    Sufficient statistic X'y.
#' @param yy    Sufficient statistic y'y.
#' @param sN    Number of Observations.
#' @param c     Spike variance parameter.
#' @param hb    Point estimate of the coefficients.
#' @param hs    Point estimates of the variance of the error.
#' @param hl    Point estimate for the global shrinking parameter.
#' @param ht    Point estimate for the local shrinking parameters.
#' @param hd    Point estimate for the sparsity parameter.
#' @param hg    Point estimate for the indicator functions.
#' @param hlsh  Hyper-prior for the lambda parameter shape.
#' @param hlsc  Hyper-prior for the lambda parameter rate.
#' @param hds1  Hyper-prior for the sparsity parameter.
#' @param hds2  Hyper-prior for the sparsity parameter.
#' @param nmcmc Number of samples of the parameters.
#'
#' @return A list containing samples for the model parameters.
#' \describe{
#'   \item{sb}{A matrix with samples for coefficients, 1 sample per row.}
#'   \item{st}{A matrix with samples for local shrinkage parameters, 1 sample
#'    per row.}
#'   \item{ss}{Samples for sigma.}
#'   \item{sl}{Samples for the global shrinkage parameter.}
#'   \item{sd}{Samples for the sparsity parameter.}
#'   \item{sg}{Samples for the Lasso indicators.}
#' }
#'
#' @author Rene Gutierrez Marquez

#' @export

dfpSpiLasSte <- function(XX,
                         Xy,
                         yy,
                         sN,
                         c,
                         hb,
                         ht,
                         hl,
                         hs,
                         hd,
                         hg,
                         hlsh  = 1,
                         hlsc  = 1,
                         hds1  = 1,
                         hds2  = 1,
                         nmcmc){
  # Gets the number of coefficients
  p <- length(hb)

  # Batch Samples
  ## Sample for beta
  sb <- matrix(0, nmcmc, p)
  ## Sample for tau
  st <- matrix(0, nmcmc, p)
  ## Sample for lambda
  sl <- vector('numeric', nmcmc)
  ## Sample for sigma
  ss <- vector('numeric', nmcmc)
  ## Sample for gamma
  sg <- matrix(0, nmcmc, p)
  ## Sample for Theta
  sd <- vector('numeric', nmcmc)

  # Samples Beta and Tau
  ## Checks for relevant and irrelevant coefficients
  ### Current Selected Parameters
  sel <- seq(1,p)[hg > 0.5]
  ### Current Deselected Parameters
  des <- seq(1,p)[hg <= 0.5]
  ### Number of deselected coefficients
  nd  <- length(des)
  ### Number of selected coefficients
  ns  <- length(sel)
  ## Works with the Deselected Parameters
  ### Checks if there are deselected parameters
  if(nd > 0)
  {
    # Samples tau for the deselected parameters
    st[, des] <- matrix(rexp(nd * nmcmc, hl / 2), nrow = nmcmc, ncol = nd)
    # Samples beta for deselected coefficients
    for(j in 1:nd)
    {
      # Variance of Deselected Parameter
      desVar       <- hs / (XX[des[j], des[j]] + 1 / c)
      # Mean of the Beta Deselected Block
      desMea       <- desVar * (Xy[des[j]] - XX[des[j],-des[j]] %*% as.matrix(hb[-des[j]])) / hs
      # Samples the Beta Deselected Block
      sb[, des[j]] <- rnorm(nmcmc, desMea, sqrt(desVar))
      # Updates the estimate for Beta Deselected Block
      hb[des[j]]   <- mean(sb[,des[j]])
    }
  }
  ## Works with the Selected Parameters
  ### Initial Value
  b <- hb[sel]
  ### Checks if there are Selected parameters
  if(ns > 0)
  {
    for(s in 1:nmcmc)
    {
      # Samples tau
      ## Mean of Selected Block tau
      tauMea <- sqrt(hl * hs / b^2)
      ## Mean of Selected Block tau
      t      <- 1/statmod::rinvgauss(ns, tauMea, hl)
      # Samples beta
      ## Variance Selected Parameters
      if(ns == 1)
      {
        # Precision Matrix Selected Block
        selVar <- solve( XX[sel,sel] + 1/t)
      } else {
        # Precision Matrix Selected Block
        selVar <- solve( XX[sel,sel] + diag(1/t))
      }
      ## Mean Selected Parameters
      selMea <- selVar %*% (Xy[sel] - XX[sel,des] %*% as.matrix(hb[des]))
      ## Samples the Beta Selected Block
      b      <- mvtnorm::rmvnorm(1, selMea, selVar)
      ### Saves the Samples
      sb[s, sel] <- b
      st[s, sel] <- t
    }
  }
  ## Sets the mean for tau and beta
  ht <- colMeans(st) # Point Estimate tau2
  hb <- colMeans(sb) # Point Estimate beta

  ## Shape lambda2
  shal <- hlsh + p
  ## Scale lambda2
  scal <- sum(ht) / 2 + hlsc

  ## Shape sigma2
  shas <- (sN + p - 1) / 2
  ## Scale sigma2
  scas <- matrix((yy - 2 * t(hb) %*% Xy + t(hb) %*% XX %*% hb + sum(hb * hb / ht) ) / 2)

  # Samples Lambda, Gamma Sigma and Theta
  for(k in 1:nmcmc)
  {
    # Samples Lambda
    ## Shape lambda2
    shal <- hlsh + p
    ## Scale lambda2
    scal <- sum(ht) / 2 + hlsc
    ## Samples Lambda
    l    <- 1 / rgamma(n     = 1,
                       shape = shal,
                       rate  = scal)

    # Samples Sigma2
    ## Shape sigma2
    shas <- (sN + p - 1) / 2
    ## Scale sigma2
    scas <- matrix((yy - 2 * t(hb) %*% Xy + t(hb) %*% XX %*% hb + sum(hb * hb / ht) ) / 2)
    ## Samples
    s    <- 1 / rgamma(n     = 1,
                       shape = shas,
                       rate  = scas)

    # Samples theta
    ## Shape 1 theta
    sh1d <- hds1 + sum(hg)
    ## Shape 2 theta
    sh2d <- hds2 + p - sum(hg)
    ## Samples
    d    <- rbeta(1, sh1d, sh2d)

    # Saves the samples
    ## Sample for sigma2
    ss[k]  <- s
    ## Sample for lambda2
    sl[k]  <- l
    ## Sample for theta
    sd[k]  <- d
  }
  ## Point Estimate
  ### Point Estimate lambda2
  hl <- mean(sl)
  ### Point Estimate
  hs <- mean(ss)
  ### Point Estimate theta
  hd <- mean(sd)

  # Samples each gamma
  for(j in 1:p)
  {
    # Computes the probability parameter
    ## Proportional to the Probability of being selected
    psel <- hd * exp(-hb[j]^2 / (2 * hs * ht[j])) / sqrt(hs * ht[j])
    ## Proportional to the Probability of being deselected
    pdes <- (1 - hd) * exp(-hb[j]^2 / (2 * hs * c)) / sqrt(hs * c)
    ## Samples gamma
    sg[, j] <- rbinom(nmcmc, 1, psel / (psel + pdes))
  }
  # Point Estimate Gamma
  hg <- round(colMeans(sg))

  return(list(sb = sb, st = st, ss = ss, sl = sl, sd = sd, sg = sg))
}
