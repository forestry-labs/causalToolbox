#' @include XRF.R


#' @title Autotuning for X-Learner with honest RF for both stages
#' @name X_RF_autotune_hyperband
#' @rdname X_RF_autotune_hyperband
#' @description This function tunes
#' @param feat A data frame of all the features.
#' @param tr A numeric vector contain 0 for control and 1 for treated variables.
#' @param yobs A numeric vector containing the observed outcomes.
#' @param sample.fraction TODO: Add Description
#' @param num_iter Number of iterations.
#' @param eta Downsampling rate. Default value is 3.
#' @param verbose if tuning process in verbose mode
#' @param seed A random seed.
#' @param nthread Number of threads used to work in parallel. 
#' @seealso \code{\link{X_RF_autotune_simple}}, \code{\link{X_RF_autotune_gpp}},
#' @examples
#'   set.seed(14236142)
#'   feat <- iris[, -1]
#'   tr <- rbinom(nrow(iris), 1, .5)
#'   yobs <- iris[, 1]
#'   # train a
#'   xl_gpp <- X_RF_autotune_hyperband(feat, tr, yobs, verbose = FALSE)
#'   # computes the CATE and confidence intervals for CATE
#'   EstimateCate(xl_gpp, feat)
#'   CateCI(xl_gpp, feat, B = 5, verbose = FALSE)
#' @export X_RF_autotune_hyperband
X_RF_autotune_hyperband <-
  function(feat,
           tr,
           yobs,
           sample.fraction = 0.75,
           num_iter = 3 ^ 8,
           eta = 3,
           verbose = TRUE,
           seed = 24750371,
           nthread = 0) {
    feat <- as.data.frame(feat)

    hyperparameter_list <- list()
    base_learners <- list()

    yobs_0 <- yobs[tr == 0]
    yobs_1 <- yobs[tr == 1]

    X_0 <- feat[tr == 0,]
    X_1 <- feat[tr == 1,]

    m_0 <-
      forestry::autoforestry(
        x = X_0,
        y = yobs_0,
        sampsize = floor(nrow(X_0) * sample.fraction),
        num_iter = num_iter,
        eta = eta,
        verbose = verbose,
        seed = seed,
        nthread = nthread
      )

    m_1 <-
      forestry::autoforestry(
        x = X_1,
        y = yobs_1,
        sampsize = floor(nrow(X_1) * sample.fraction),
        num_iter = num_iter,
        eta = eta,
        verbose = verbose,
        seed = seed,
        nthread = nthread
      )

    if (verbose) {
      print("Done with the first stage.")
    }
    r_0 <- predict(m_1, X_0) - yobs_0
    r_1 <- yobs_1 - predict(m_0, X_1)

    m_tau_0 <-
      forestry::autoforestry(
        x = X_0,
        y = r_0,
        sampsize = floor(nrow(X_0) * sample.fraction),
        num_iter = num_iter,
        eta = eta,
        verbose = verbose,
        seed = seed,
        nthread = nthread
      )

    m_tau_1 <-
      forestry::autoforestry(
        x = X_1,
        y = r_1,
        sampsize = floor(nrow(X_1) * sample.fraction),
        num_iter = num_iter,
        eta = eta,
        verbose = verbose,
        seed = seed,
        nthread = nthread
      )

    if (verbose) {
      print("Done with the second stage.")
    }

    m_prop <-
      forestry::forestry(x = feat,
               y = tr,
               ntree = 500)
    if (verbose) {
      print("Done with the propensity score estimation.")
    }

    hyperparameter_list <- get_hyper_parameter_list(m_0, m_1, m_tau_0,
                                                   m_tau_1, m_prop, feat,
                                                   nthread)
    return(
      new(
        "X_RF",
        feature_train = feat,
        tr_train = tr,
        yobs_train = yobs,
        m_0 = m_0,
        m_1 = m_1,
        m_tau_0 = m_tau_0,
        m_tau_1 = m_tau_1,
        m_prop = m_prop,
        hyperparameter_list = hyperparameter_list,
        creator = function(feat, tr, yobs) {
          X_RF_fully_specified(feat = feat,
                               tr = tr,
                               yobs = yobs,
                               hyperparameter_list = hyperparameter_list,
                               verbose = verbose)
        }
      )
    )
  }



get_hyper_parameter_list <-
  function(m_0, m_1, m_tau_0, m_tau_1, m_prop, feat, nthread) {
    hyperparameter_list <- list(
      "general" = list("predmode" = "propmean", "nthread" = nthread),
      "l_first_0" = get_hyper_parameter_list_for_this_learner(m_0, feat),
      "l_first_1" = get_hyper_parameter_list_for_this_learner(m_1, feat),
      "l_second_0" = get_hyper_parameter_list_for_this_learner(m_tau_0, feat),
      "l_second_1" = get_hyper_parameter_list_for_this_learner(m_tau_1, feat),
      "l_prop" = get_hyper_parameter_list_for_this_learner(m_prop, feat)
    )
    return(hyperparameter_list)
  }

get_hyper_parameter_list_for_this_learner <-
  function(rfm, feat) {
    return(
      list(
        "relevant_Variable" = 1:ncol(feat),
        "ntree" = rfm@ntree,
        "replace" = rfm@replace,
        "sampsize" = rfm@sampsize,
        "mtry" = rfm@mtry,
        "nodesizeSpl" = rfm@nodesizeSpl,
        "nodesizeAvg" = rfm@nodesizeAvg,
        "splitratio" = rfm@splitratio,
        "middleSplit" = rfm@middleSplit
      )
    )
  }
