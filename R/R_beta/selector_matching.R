#' @include helper_functions.R
#' This is included to have the get get_CV_sizes.R function

#' gof_matching
#' @title gof_matching
#' @name gof_matching
#' @description Match each unit in one group(treatment/control) to units in
#' another group(treatment/control) based on the chosen estimand. For each iteration
#' of cross validation, assign treatment lablels (0 or 1) with probability given 
#' by averaged propensity score for each pair of match and estimate CATE. Return
#' the mean squared error and its standard deviation of CATE from the true difference. 
#' @param estimand ATE, ATT or ATC
#' @param replace logical. If TRUE then samples will be replaced while matching
#' @inheritParams gof_transformed
#' @return mean(error) and sd(error)
#' @import Matching
#' @import ranger
#' @export gof_matching
gof_matching <- function(feat,
                         yobs,
                         tr,
                         estimator,
                         estimand = 'ATT',
                         k = 2,
                         replace = TRUE,
                         emin = 1e-5,
                         verbose = FALSE) {

  # catch nonsensible specifications
  if (emin <= 0 | emin >= 0.5) {
    stop("0 < emin < 0.5")
  }
  catch_error(feat, yobs, tr, k)
  
  # -------------------------------------------------------------------------
  # Estimate propensity score
  pscore_pred <- estimate_pscore(feat, tr, emin)
  
  # --------------------------------------------------------------------------
  # Create the matched data set. We match on the features and we create a
  # data set which has features|Y(0)|Y(1)|pscore
  m <- Matching::Match(
    Y = yobs,
    Tr = tr,
    X = feat,
    estimand = estimand,
    replace = replace,
    ties = FALSE
  )
  
  Y1 <- m$mdata$Y[m$mdata$Tr == 1]
  Y0 <- m$mdata$Y[m$mdata$Tr == 0]
  pscore_m <-
    (pscore_pred[m$index.treated] + pscore_pred[m$index.control]) / 2
  if (estimand == 'ATT') {
    feat_m <- feat[m$index.treated, ]
  }
  if (estimand == 'ATC') {
    feat_m <- feat[m$index.control, ]
  }
  if (estimand == 'ATE') {
    feat_m <- feat[m$index.control, ]
  }
  
  # ----------------------------------------------------------------------------
  # CV compute CATEs
  # Situation:
  # feat_m : contains the covariates, Y1, Y0 :potential outcomes, pscore_m : ps
  n_matched <- length(Y1)
  
  # Create CV idxes
  cv_idx <- sample(rep(x = 1:k, times = get_CV_sizes(n_matched, k)))
  
  cate_est <- rep(NA, n_matched) # will contain the estimates
  for (b in 1:k) {
    if (verbose) {
      print(paste("Running", b, "out of", k, "CV fold."))
    }
    # get train and test set -- training set is everything but fold i
    train_idx <- which(cv_idx != b)
    test_idx <- which(cv_idx == b)
    
    # construct the treatment column
    feat_b <- feat_m[train_idx, ]
    tr_b <- rbinom(length(train_idx), 1, pscore_m[train_idx])
    yobs_b <- ifelse(tr_b == 1, Y1[train_idx], Y0[train_idx])
    
    if (sum(tr_b) == 0 | sum(tr_b) == length(tr_b)) {
      stop("In the CV all units were in the treated group or all units are in
           the control group")
      # TODO: We might want to do conditional sampling here such that the 
      # proportion of treated units is equal to the poroportion of treated units
      # in the original data set.
    }
    
    # Estimate CATE with the given learner function
    estimator_trained <- estimator(feat = feat_b,
                                   tr = tr_b,
                                   yobs = yobs_b)
    cate_est[test_idx] <- EstimateCate(estimator_trained, feat_m[test_idx, ])
  }
  
  # --------------------------------------------------------------------------
  # Compute the ITE = Y(1) - Y(0)
  ITE <- Y1 - Y0
  
  # Calcualte the Goodness-of-Fit
  mse <- mean((ITE - cate_est) ^ 2)
  sd_err <- sd((ITE - cate_est) ^ 2) / sqrt(n_matched)
  
  # --------------------------------------------------------------------------
  return(c(mse, sd_err))
}
