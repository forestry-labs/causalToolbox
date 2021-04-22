#' @include CATE_estimators.R
#' @include helper_functions.R
#' @include MRF.R
#' @import BART

## the standard M_BART object
setClass(
  "M_BART",
  contains = "MetaLearner",
  slots = list(
    feature_train = "data.frame",
    tr_train = "numeric",
    yobs_train = "numeric",
    ndpost = "numeric",
    ntree = "numeric",
    nthread = "numeric",
    mu0.BART = "list", 
    mu1.BART = "list", 
    e.BART = "list",
    tau.BART = "list",
    creator = "function"
  )
)

#' @rdname Mlearners
#' @description \code{M_BART} is an implementation of the Modified Outcome Estimator with
#'   Bayesian Additive Regression Trees (Chipman et al. 2010) as the base learner.
#' @param mu.BART,e.BART,tau.BART Hyperparameters of the BART functions for the
#'   control and treated group. (Use \code{?BART::mc.wbart} for a detailed
#'   explanation of their effects.)
#' @export M_BART
#' @inherit M_RF details
#' @inherit T_BART
#' @export
M_BART <-
  function(feat,
           tr,
           yobs,
           ndpost = 1200,
           ntree = 200,
           nthread = 1,
           mu.BART = list(
             sparse = FALSE,
             theta = 0,
             omega = 1,
             a = 0.5,
             b = 1,
             augment = FALSE,
             rho = NULL,
             usequants = FALSE,
             cont = FALSE,
             sigest = NA,
             sigdf = 3,
             sigquant = 0.90,
             k = 2.0,
             power = 2.0,
             base = .95,
             sigmaf = NA,
             lambda = NA,
             numcut = 100L,
             nskip = 100L
           ),
           e.BART = list(
             sparse = FALSE,
             theta = 0,
             omega = 1,
             a = 0.5,
             b = 1,
             augment = FALSE,
             rho = NULL,
             usequants = FALSE,
             cont = FALSE,
             sigest = NA,
             sigdf = 3,
             sigquant = 0.90,
             k = 2.0,
             power = 2.0,
             base = .95,
             sigmaf = NA,
             lambda = NA,
             numcut = 100L,
             nskip = 100L
           ),
           tau.BART = list(
             sparse = FALSE,
             theta = 0,
             omega = 1,
             a = 0.5,
             b = 1,
             augment = FALSE,
             rho = NULL,
             usequants = FALSE,
             cont = FALSE,
             sigest = NA,
             sigdf = 3,
             sigquant = 0.90,
             k = 2.0,
             power = 2.0,
             base = .95,
             sigmaf = NA,
             lambda = NA,
             numcut = 100L,
             nskip = 100L
           )) {
    feat <- as.data.frame(feat)
    
    new(
      "M_BART",
      feature_train = feat,
      tr_train = tr,
      yobs_train = yobs,
      ndpost = ndpost,
      ntree = ntree,
      nthread = nthread,
      mu0.BART = mu.BART, 
      mu1.BART = mu.BART,
      e.BART = e.BART,
      tau.BART = tau.BART,
      creator = function(feat, tr, yobs) {
        M_BART(
          feat = feat,
          tr = tr,
          yobs = yobs,
          ndpost = ndpost,
          ntree = ntree,
          nthread = nthread,
          mu.BART = mu.BART,
          e.BART = e.BART,
          tau.BART = tau.BART
        )
      }
    )
  }


#' EstimateCate-M_BART
#' @rdname  EstimateCate-M_BART
#' @inherit EstimateCate
#' @exportMethod EstimateCate
#' @import stats
setMethod(
  f = "EstimateCate",
  signature = "M_BART",
  definition = function(theObject,
                        feature_new,
                        verbose = FALSE,
                        return_CI = FALSE) {
    yobs <- theObject@yobs_train
    feat <- theObject@feature_train
    tr <- theObject@tr_train
    ## Constructing regression adjustment --------------------------------------

    pred_matrix_f_1 <- get_BART_pred(x.train = feat[as.logical(tr),], 
                                     y.train = yobs[as.logical(tr)], 
                                     x.test = feat, 
                                     ndpost = theObject@ndpost, 
                                     ntree = theObject@ntree, 
                                     nthread = theObject@nthread, 
                                     hyperparam = theObject@mu0.BART)
    
    mu_hat_1 <- colMeans(pred_matrix_f_1)
    
    pred_matrix_f_0 <- get_BART_pred(x.train = feat[!tr,], 
                                     y.train = yobs[!tr],
                                     x.test = feat,
                                     ndpost = theObject@ndpost,
                                     ntree = theObject@ntree,
                                     nthread = theObject@nthread,
                                     hyperparam = theObject@mu1.BART)
    mu_hat_0 <- colMeans(pred_matrix_f_0)
    
    prop_matrix <- get_BART_pred(x.train = feat, 
                                 y.train = as.numeric(factor(tr)),
                                 x.test = feat,
                                 ndpost = theObject@ndpost,
                                 ntree = theObject@ntree,
                                 nthread = theObject@nthread,
                                 hyperparam = theObject@e.BART)
    propensity_score_hat <- pnorm(apply(prop_matrix, 2, mean))
    
    modified_outcome_ra <- 
      (tr - propensity_score_hat) / 
      (propensity_score_hat * (1 - propensity_score_hat)) *
      (yobs - mu_hat_1 * (1 - propensity_score_hat) - 
         mu_hat_0 * propensity_score_hat)
    
    if (verbose)
      print("Done estimating regression adjustment.")
    
    ## Computing tau  ----------------------------------------------------------
    
    pred_matrix_tau <- get_BART_pred(x.train = feat, 
                                     y.train = modified_outcome_ra,
                                     x.test = feature_new,
                                     ndpost = theObject@ndpost,
                                     ntree = theObject@ntree,
                                     nthread = theObject@nthread,
                                     hyperparam = theObject@tau.BART)
    
    
    if (verbose) print("Done with the second stage.")
    
    ## Generating the output ---------------------------------------------------
    
    tau_hat <- apply(pred_matrix_tau, 2, mean)
    
    if (return_CI) {
      CI <- t(apply(pred_matrix_tau, 
                    2, 
                    function(x) {quantile(x, probs = c(.025, 0.975))}))
      
      to_return <- as.data.frame(cbind(tau_hat, CI))
      row.names(to_return) <- 1:nrow(to_return)
      colnames(to_return) <- c('pred', 'X5.', 'X95.')
      return(to_return)
    } else{
      return(tau_hat)
    }
  }
)


#' CateCI-M_BART
#' @rdname CateCI-M_BART
#' @export
setMethod(
  f = "CateCI",
  signature = "M_BART",
  definition = function(theObject,
                        feature_new,
                        verbose = FALSE) {
    return(
      EstimateCate(
        theObject,
        feature_new,
        verbose = verbose,
        return_CI = TRUE
      )
    )
  }
)


