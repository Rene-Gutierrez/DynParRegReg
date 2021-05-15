#' Statistics Creator
#'
#' This function creates the performance statistics used in Gutierrez &
#' Guhaniyogi.
#'
#' @param y Responses
#' @param X Predictors
#' @param b Samples for the coefficients.
#' @param s Samples for the variance of the error.
#'
#' @return A list containing selected parameter samples, sufficient statistics,
#' and parameter partition.
#' \describe{
#'   \item{preCreCov}{Prediction 95 Credible Interval Coverage.}
#'   \item{preSquErr}{Mean Prediction Square Error.}
#'   \item{preCreLen}{Prediction 95 Credible Interval Length.}
#'   \item{preIntSco}{Prediction 95 Credible Interval Score.}
#' }
#'
#' @author Rene Gutierrez Marquez

#' @export

pre.sta = function(y, X, b, s)
{
  #  Gets the Dimensions
  ## Number of Observations
  n     <- length(y)
  ## Number of Simulations
  nmcmc <- length(s)
  # Computes Prediction Samples for the Test
  ## Computes the Prediction Mean
  preMea <- X %*% t(b)
  ## Starts the Prediction Samples
  pre    <- matrix(0, n, nmcmc)
  for(i in 1:n)
  {
    # For every observation we compute in the Prediction Set
    ## Samples predictions
    pre[i,] <- rnorm(nmcmc, preMea[i,], sqrt(s) )
  }

  # Prediction Statistics
  ## Prediction Mean
  preMea <- rowMeans(pre)
  ## Credibility Interval
  # Lower Quantile
  preLowQua <- apply(t(pre), 2, quantile, 0.025 )
  # Upper Quantile
  preUppQua <- apply(t(pre), 2, quantile, 0.975 )
  ## Credibility Interval length
  preCreLen <- mean(preUppQua - preLowQua)/sqrt(var(y))
  ## Coverage of Credibility Interval
  lg <- as.numeric(preLowQua < y)
  ug <- as.numeric(y < preUppQua)
  preCreCov <- mean(lg * ug)
  ## Square Error
  preSquErr <- mean((y - preMea)^2)/var(y)
  ## Interval Score
  preIntSco <- (preUppQua - preLowQua) + (40) * (preLowQua - y) * (1-lg) + (40) * (y - preUppQua) * (1-ug)
  preIntSco <- mean(preIntSco)

  ### Returns Results
  return( c(preCreCov, preSquErr, preCreLen, preIntSco))
}
