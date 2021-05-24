#' Performs DFP for the Bayesian Lasso Prior
#'
#' This function runs DFP using the Bayesian Lasso prior on data batches
#' previously stored in the indicated paths. It saves the performance statistics
#' indicated in the paper as well as the current point estimates and sufficient
#' statistics used for the following batch. It returns samples for selected
#' parameters for monitoring, as well as the updated sufficient statistics.
#' Notice that in this case, the partition is implied in the spike indicator
#' functions.
#'
#' @param staBat Starting batch number.
#' @param endBat End batch number. As default is set as the starting batch.
#' @param XDatPat Path to the batch predictors. A default is provided as a
#' suggestion.
#' @param yDatPat Path to the batch responses.
#' @param staPat  Path to the statistics of previous batches. If no file has
#' been created a new one is created with the path provided.
#' @param curParPat Path to the current parameter point estimates, partition and
#' sufficient statistics.
#' @param savPat Path to the file to save the parameters requested.
#' @param suffix Suffix to be added to the file paths.
#' @param cores  Number of cores to be used in the parallel computations. Set as
#' default as half the number of cores detected.
#' @param staFla A boolean indicating if the statistics will be saved.
#' @param savCurPar A boolean indicating if the updated parameters will be
#' saved.
#' @param nmcmc Number of samples of the parameters.
#' @param hlsh  Hyper-prior for the lambda parameter shape. Default is 1.
#' @param hlsc  Hyper-prior for the lambda parameter rate. Default is 1.
#' @param hds1  Hyper-prior for the theta (sparsity) parameter. Default is 1.
#' @param hds2  Hyper-prior for the theta (sparsity) parameter. Default is 1.
#' @param progress A boolean indicating if a progress bar is desired. Default is
#'  TRUE.
#' @param staPar Starting value for the parameters. Used only for the first
#' batch.
#' @param savBol Boolean indicating if the selected parameters will be saved.
#' @param savCoe Vector indicating the index of the coefficients to be saved.
#' This will also saved the associated parameters to each coefficient.
#' @param datPer Vector containing the order in which the data is to be
#' processed.
#' @param c Spike parameter value.
#' @param iniPar Initial partition.
#'
#' @return A list containing selected parameter samples, sufficient statistics,
#' and parameter partition.
#' \describe{
#'   \item{sb}{A matrix with samples for the selected coefficients, 1 sample per
#'   row.}
#'   \item{st}{A matrix with samples for the selected local shrinkage
#'   parameters associated to the selected coefficients, 1 sample per row.}
#'   \item{ss}{Samples for sigma.}
#'   \item{sl}{Samples for the global shrinkage parameter.}
#'   \item{sd}{Samples for the sparsity parameter.}
#'   \item{sg}{Samples for the Lasso indicators.}
#'   \item{XX}{Updated sufficient statistic X'X.}
#'   \item{Xy}{Updated sufficient statistic X'y.}
#'   \item{yy}{Updated sufficient statistic y'y.}
#'   \item{sN}{Updated sufficient statistic for the number of observations
#'   processed so far.}
#' }
#'
#' @author Rene Gutierrez Marquez

#' @export

dfpSpiLas <- function(staBat    = 1,
                      endBat    = staBat,
                      XDatPat   = "./dat/datX-complete-",
                      yDatPat   = "./dat/daty-complete-",
                      staPat    = "./out/sta.RData",
                      curParPat = "./out/curPar.RData",
                      savPat    = "./out/savPar.RData",
                      suffix    = "complete-dfp-BayLas",
                      cores     = parallel::detectCores() / 2,
                      staFla    = TRUE,
                      savCurPar = FALSE,
                      nmcmc     = 100,
                      hlsh      = 1,
                      hlsc      = 1,
                      hds1      = 1,
                      hds2      = 1,
                      progress  = TRUE,
                      staPar,
                      savBol    = TRUE,
                      savCoe,
                      datPer    = seq(1:endBat),
                      c,
                      iniPar){

  # Parallel Set-Up
  print(paste0("Number of Cores = ", cores))
  doParallel::registerDoParallel(cores)

  # Progress Bar Initialization
  if(progress){
    pb <- txtProgressBar(min     = 0,
                         max     = 1,
                         initial = 0,
                         style   = 3,
                         width   = 72)
  }

  # Saved Parameters Initialization and Retrieve
  ## Checks if the Parameter Saving Is requested
  if(savBol){
    # Checks if a saved parameter Exists
    if(file.exists(savPat)){
      # Retrieves the Saved Parameters
      savPar <- get(load(file = savPat))
    } else {
      # Initializes a Saved Parameter List
      savPar <- list()
    }
  }

  # Stats Initialization
  ## Checks if it is the first Batch
  if(staBat == 1){
    # Initializes the Parameters
    sta <- list()
  } else {
    # Recalls the Stats if Available
    if(file.exists(staPat)){
      sta <- get(load(file = staPat))
    } else {
      sta <- list()
    }
  }

  # Sufficient Statistics Initialization
  ## Checks if it is the first Batch
  if(staBat == 1){
    # Initializes the Sufficient Statistics
    XX <- 0
    Xy <- 0
    yy <- 0
    sN <- 0
  } else {
    # Gets the Sufficient Statistics
    curPar <- get(load(file = curParPat))
    # Sets the Current Sufficient Statistics
    XX <- curPar$XX
    Xy <- curPar$Xy
    yy <- curPar$yy
    sN <- curPar$sN
  }

  # Variable Initialization
  ## Checks if it is the first Batch
  if(staBat == 1){
    # Initializes the Variables
    hb <- staPar$hb
    ht <- staPar$ht
    hl <- staPar$hl
    hs <- staPar$hs
    hd <- staPar$hd
    hg <- staPar$hg
  } else {
    # Gets the Current Point Estimates
    curPar <- get(load(file = curParPat))
    # Sets the Sufficient Statistics
    hb <- curPar$hb
    ht <- curPar$ht
    hl <- curPar$hl
    hs <- curPar$hs
    hd <- curPar$hd
    hg <- curPar$hg
  }

  # Performs Spike and Lasso with DFP
  for(i in staBat:endBat){
    # Obtains the Data Batch
    ## Regression Matrix
    X <- get(load(file=paste0(XDatPat, per[i], '.RData')))
    ## Dependent Variable
    y <- get(load(file=paste0(yDatPat, per[i], '.RData')))

    # Time Start
    beg = Sys.time()

    # Updates the Sufficient Statistics
    XX <- XX + t(X) %*% X
    Xy <- Xy + t(X) %*% y
    yy <- yy + sum(y^2)
    sN <- sN + length(y)

    # Performs Spike & Lasso with DFP
    samOut <- dfpSpiLasSte(XX    = XX,
                           Xy    = Xy,
                           yy    = yy,
                           sN    = sN,
                           hb    = hb,
                           ht    = ht,
                           hl    = hl,
                           hs    = hs,
                           hd    = hd,
                           hg    = hg,
                           hlsh  = hlsh,
                           hlsc  = hlsc,
                           hds1  = hds1,
                           hds2  = hds2,
                           nmcmc = nmcmc,
                           c     = c)
    # Obtains the Samples
    sb <- samOut$sb
    st <- samOut$st
    sl <- samOut$sl
    ss <- samOut$ss
    sd <- samOut$sd
    sg <- samOut$sg
    # Updates the Point Estimates
    hb <- colMeans(sb)
    ht <- colMeans(st)
    hl <- mean(sl)
    hs <- mean(ss)
    hd <- mean(sd)
    hg <- colMeans(sg)

    # Runtime
    tim <- Sys.time() - beg

    # Saves the Requested Coefficients and Associated Parameters
    ## Checks if saved is requested
    if(savBol){
      savPar[[i]] <- list(sb     = sb[, savCoe],
                          st     = st[, savCoe],
                          sl     = sl,
                          ss     = ss,
                          sl     = sl,
                          sd     = sd,
                          sg     = sg[, savCoe],
                          hg     = hg,
                          savCoe = savCoe)
    }
    ## Saves to File
    save(savPar, file = savPat)

    if(is.na(hs)){
      break
    }

    # Saves the Current Point Estimates and Partition
    curPar <- list(hb = hb,
                   ht = ht,
                   hl = hl,
                   hs = hs,
                   hd = hd,
                   hg = hg,
                   XX = XX,
                   Xy = Xy,
                   yy = yy,
                   sN = sN)

    ## Saves to File
    if(savCurPar){
      save(curPar, file = curParPat)
    }

    # Computes the Relevant Statistics for Evaluation
    ## Checks if the Statistics are required
    if(staFla){
      # Gets the Stats for period i
      perSta <- pre.sta(y, X, sb, ss)
      # Saves the stats
      sta[[i]] <- c(perSta, tim)
      # Saves to File
      save(sta, file = staPat)
    }

    ### Progress Bar Display
    if(progress){
      setTxtProgressBar(pb    = pb,
                        value = (i - staBat + 1) / (endBat - staBat + 1))
    }
  }
  ### Closes the Progress Bar
  if(progress){
    close(pb)
  }

  return(list(sb = sb, st = st, sl = sl, ss = ss, sd = sd, sg = sg, XX = XX, Xy = Xy, yy = yy, sN = sN))
}
