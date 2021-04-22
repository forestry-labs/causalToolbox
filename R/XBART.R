#' @include CATE_estimators.R
#' @include helper_functions.R
#' @import BART

## the standard Xlearner object with random forest
setClass(
  "X_BART",
  contains = "MetaLearner",
  slots = list(
    feature_train = "data.frame",
    tr_train = "numeric",
    yobs_train = "numeric",
    predmode = "character",
    ndpost = "numeric",
    ntree = "numeric",
    nthread = "numeric",
    mu0.BART = "list", 
    mu1.BART = "list",
    tau0.BART = "list", 
    tau1.BART = "list",
    e.BART = "list",
    creator = "function"
  )
)


#' @rdname Xleaners
#' @description X_BART is an implementation of the X-learner with Bayesian
#'   Additive Regression Trees (Chipman et al. 2010) at the first and second stage
#' @param mu.BART,tau.BART,e.BART hyperparameters of the BART functions for the
#'   estimates of the first and second stage and the propensity score. Use
#'   \code{?BART::mc.wbart} for a detailed explanation of their effects.
#' @inherit X-Learner
#' @inheritParams T_BART
#' @export 
X_BART <-
  function(feat,
           tr,
           yobs,
           predmode = "pscore",
           nthread = 1,
           ndpost = 1200,
           ntree = 200,
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
           )) {
    feat <- as.data.frame(feat)
    new(
      "X_BART",
      feature_train = feat,
      tr_train = tr,
      yobs_train = yobs,
      predmode = predmode,
      ndpost = ndpost,
      ntree = ntree,
      nthread = nthread,
      mu0.BART = mu.BART, 
      mu1.BART = mu.BART,
      tau0.BART = tau.BART, 
      tau1.BART = tau.BART,
      e.BART = e.BART,
      creator = function(feat, tr, yobs) {
        X_BART(feat,
               tr,
               yobs,
               predmode = predmode,
               ndpost = ndpost,
               ntree = ntree, 
               mu.BART = mu.BART,
               tau.BART = tau.BART, 
               e.BART = e.BART)
      }
    )
  }


#' EstimateCate-X_BART
#' @rdname EstimateCate-X_BART
#' @inherit EstimateCate
#' @import stats
#' @export
setMethod(
  f = "EstimateCate",
  signature = "X_BART",
  definition = function(theObject,
                        feature_new,
                        verbose = FALSE,
                        return_CI = FALSE)
  {
    yobs <- theObject@yobs_train
    feat <- theObject@feature_train
    tr <- theObject@tr_train
    ndpost <- theObject@ndpost
    predmode <- theObject@predmode

    yobs_0 <- yobs[tr == 0]
    X_0 <- feat[tr == 0, ]
    yobs_1 <- yobs[tr == 1]
    X_1 <- feat[tr == 1, ]

    # First stage --------------------------------------------------------------

    n_1 <- sum(tr)
    n_0 <- sum(1 - tr)

    if (return_CI) {
      # if CI should be returned, we also need to estimate the uncertainty, we
      # had in the first stage.
      f_0_test_set <- rbind(X_1, feature_new)
      f_1_test_set <- rbind(X_0, feature_new)
    } else{
      f_0_test_set <- X_1
      f_1_test_set <- X_0
    }


    pred_matrix_f_0 <- get_BART_pred(x.train = X_0, 
                                     y.train = yobs_0, 
                                     x.test = f_0_test_set, 
                                     ndpost = theObject@ndpost, 
                                     ntree = theObject@ntree, 
                                     nthread = theObject@nthread, 
                                     hyperparam = theObject@mu0.BART)


    mu_hat_1 <- apply(pred_matrix_f_0[ ,1:n_1], 2, mean)

    pred_matrix_f_1  <- get_BART_pred(x.train = X_1, 
                                      y.train = yobs_1, 
                                      x.test = f_1_test_set, 
                                      ndpost = theObject@ndpost, 
                                      ntree = theObject@ntree, 
                                      nthread = theObject@nthread, 
                                      hyperparam = theObject@mu1.BART)

    mu_hat_0 <- apply(pred_matrix_f_1[ ,1:n_0], 2, mean)

    if (verbose)
      print("Done with the first stage.")

    
    # second stage -------------------------------------------------------------
    
    D_1 <- yobs_1 - mu_hat_1
    D_0 <- mu_hat_0 - yobs_0

    pred_matrix_s_1 <- get_BART_pred(x.train = X_1, 
                                     y.train = D_1, 
                                     x.test = feature_new, 
                                     ndpost = theObject@ndpost, 
                                     ntree = theObject@ntree, 
                                     nthread = theObject@nthread, 
                                     hyperparam = theObject@tau1.BART)

    tau_hat_1 <- apply(pred_matrix_s_1, 2, mean)

    pred_matrix_s_0  <- get_BART_pred(x.train = X_0, 
                                      y.train = D_0, 
                                      x.test = feature_new, 
                                      ndpost = theObject@ndpost, 
                                      ntree = theObject@ntree, 
                                      nthread = theObject@nthread, 
                                      hyperparam = theObject@tau0.BART)

    tau_hat_0 <- apply(pred_matrix_s_0, 2, mean)

    if (verbose)
      print("Done with the second stage.")
    
    # Combining the two --------------------------------------------------------
    if (predmode == "pscore") {
      prop_matrix   <- get_BART_pred(
        x.train = feat,
        y.train = as.numeric(factor(tr)),
        x.test = feature_new,
        ndpost = theObject@ndpost,
        ntree = theObject@ntree,
        nthread = theObject@nthread,
        hyperparam = theObject@e.BART
      )
      
      g_weights <- pnorm(apply(prop_matrix, 2, mean))
      if (verbose)
        print("Done with the propensity score estimation.")
    } else if (predmode == "1/2") {
      g_weights <- 1 / 2
    } else if (predmode == "constant p-score") {
      g_weights <- sum(tr) / length(tr)
    } else if (predmode == "only control") {
      g_weights <- 0
    } else if (predmode == "only treated") {
      g_weights <- 1
    } else if (predmode == "variance") {
      var_s_0 <- apply(pred_matrix_s_0, 2, var) / ndpost
      var_s_1 <- apply(pred_matrix_s_1, 2, var) / ndpost
      g_weights <- var_s_1 / (var_s_1 + var_s_0)
    }

    # Combining the two --------------------------------------------------------

    pred <- g_weights * tau_hat_0 + (1 - g_weights) * tau_hat_1


    if (return_CI) {
      # Variance from the first stage:

      #TODO : This is a very concervatice way of getting CI, one could directly
      # use the MCMC samples and combine them or look at the convoultion of the
      # empericals.

      n_new <- nrow(feature_new)

      get_CI_mu0 <- t(apply(pred_matrix_f_0[ ,(n_1 + 1):(n_1 + n_new)], 2,
                            function(x) quantile(x, probs = c(.05, 0.95))))
      get_CI_mu1 <- t(apply(pred_matrix_f_1[ ,(n_0 + 1):(n_0 + n_new)], 2,
                            function(x) quantile(x, probs = c(.05, 0.95))))

      mu0_hat_feature_new <- apply(pred_matrix_f_0[ ,(n_1 + 1):(n_1 + n_new)],
                                   2, mean)
      mu1_hat_feature_new <- apply(pred_matrix_f_1[ ,(n_0 + 1):(n_0 + n_new)],
                                   2, mean)


      # Variance from the second stage:
      get_CI_0 <- t(apply(pred_matrix_s_0, 2, function(x)
        quantile(x, probs = c(.05, 0.95))))
      get_CI_1 <- t(apply(pred_matrix_s_1, 2, function(x)
        quantile(x, probs = c(.05, 0.95))))

      CI_comb <-
        g_weights *       (get_CI_0 - get_CI_mu1[ ,2:1] + mu1_hat_feature_new) +
        (1 - g_weights) * (get_CI_1 - get_CI_mu0[ ,2:1] + mu0_hat_feature_new)

      to_return <- as.data.frame(cbind(pred, CI_comb))
      row.names(to_return) <- 1:nrow(to_return)
      colnames(to_return) <- c('pred','X5.','X95.')
      return(to_return)
    } else{
      return(pred)
    }
  }
)


#' CateCI-X_BART
#' @rdname CateCI-X_BART
#' @aliases CateCI,X_BART-method
#' @exportMethod CateCI
setMethod(
  f = "CateCI",
  signature = "X_BART",
  definition = function(theObject,
                        feature_new,
                        verbose = FALSE)
  {
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


