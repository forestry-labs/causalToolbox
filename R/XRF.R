# This file implements the X-Learner (https://arxiv.org/pdf/1706.03461.pdf)
# with the forestry implementation (https://github.com/soerenkuenzel/forestry)
# as base learner.
#' @include CATE_estimators.R
#' @include helper_functions.R
NULL

# X-RF class -------------------------------------------------------------------
setClass(
  "X_RF",
  contains = "MetaLearner",
  slots = list(
    m_0 = "forestry",
    m_1 = "forestry",
    m_tau_0 = "forestry",
    m_tau_1 = "forestry",
    m_prop = "forestry",
    hyperparameter_list = "list"
  )
)

# X_RF generator ---------------------------------------------------------------
#' @title X-Learners
#' @rdname Xleaners
#' @name X-Learner
#' @description X_RF is an implementation of the X-learner with Random Forests
#'   (Breiman 2001) in the first and second stage.
#' @details 
#' The X-Learner estimates the CATE in three steps:
#' \enumerate{
#'  \item
#'     Estimate the response functions 
#'     \deqn{\mu_0(x) = E[Y(0) | X = x]}
#'     \deqn{\mu_1(x) = E[Y(1) | X = x]} 
#'     using the base learner and denote the estimates as \eqn{\hat \mu_0} and
#'     \eqn{\hat \mu_1}.
#'  \item
#'     Impute the treatment effects for the individuals in the treated group,
#'     based on the control outcome estimator, and the treatment effects for the
#'     individuals in the control group, based on the treatment outcome
#'     estimator, that is,
#'     \deqn{D^1_i = Y_i(1) - \hat \mu_0(X_i)}
#'     \deqn{D^0_i = \hat \mu_1(X_i) - Y_i(0).}
#'     Now employ the base learner in two ways: using \eqn{D^1_i} as the
#'     dependent variable to obtain \eqn{\hat \tau_1(x)}, and using \eqn{D^0_i}
#'     as the dependent variable to obtain \eqn{\hat \tau_0(x)}.
#'  \item 
#'     Define the CATE estimate by a weighted average of the two estimates at
#'     Stage 2: 
#'     \deqn{\tau(x) = g(x) \hat \tau_0(x) + (1 - g(x)) \hat \tau_1(x).} 
#'     If \code{predmode = "propmean"}, then \eqn{g(x) = e(x)}, where
#'     \eqn{e(x)} is an estimate of the propensity score using the 
#'     \href{https://github.com/forestry-labs/Rforestry}{\code{Rforestry}} Random Forests
#'     version with the hyperparameters specified in \code{e.forestry}.
#'     If \code{predmode = "control"}, then \eqn{g(x) = 1}, and if 
#'     \code{predmode = "treated"}, then \eqn{g(x) = 0}.
#' }
#' @param feat A data frame containing the features.
#' @param tr A numeric vector with 0 for control and 1 for treated variables.
#' @param yobs A numeric vector containing the observed outcomes.
#' @param predmode Specifies how the two estimators of the second stage should
#'   be aggregated. Possible types are "propmean," "control," and "treated." The
#'   default is "propmean," which refers to propensity score weighting.
#' @param nthread Number of threads which should be used to work in parallel.
#' @param verbose TRUE for detailed output, FALSE for no output.
#' @param mu.forestry,tau.forestry,e.forestry A list containing the
#'   hyperparameters for the \code{Rforestry} package that are used for
#'   estimating the response functions, the CATE, and the propensity score.
#'   These hyperparameters are passed to the \code{Rforestry} package. (Please
#'   refer to the \href{https://github.com/forestry-labs/Rforestry}{Rforestry}
#'   package for a more detailed documentation of the hyperparamters.)
#'   \itemize{
#'      \item \code{relevant.Variable} Variables that are only used in the first 
#'            stage.
#'      \item \code{ntree} Numbers of trees used in the first stage.
#'      \item \code{replace} Sample with or without replacement in the first 
#'            stage.
#'      \item \code{sample.fraction} The size of total samples to draw for the 
#'            training data in the first stage.
#'      \item \code{mtry} The number of variables randomly selected in each 
#'            splitting point.
#'      \item \code{nodesizeSpl} Minimum nodesize in the first stage for 
#'            the observations in the splitting set. (See the details of the 
#'            \code{forestry} package)
#'      \item \code{nodesizeAvg} Minimum nodesize in the first stage for 
#'            the observations in the averaging set.
#'      \item \code{nodesizeStrictSpl} Minimum nodesize in the first stage for 
#'            the observations in the splitting set. (See the details of the 
#'            \code{forestry} package)
#'      \item \code{nodesizeStrictAvg} Minimum nodesize in the first stage for 
#'            the observations in the averaging set.
#'      \item \code{splitratio} Proportion of the training data used as the 
#'            splitting dataset in the first stage.
#'      \item \code{middleSplit} If true, the split value will be exactly in the 
#'            middle of two observations. Otherwise, it will take a point 
#'            based on a uniform distribution between the two observations. 
#'      \item \code{OOBhonest} If true, forestry object will use the Out of Bag
#'            honesty implemented in the \code{Rforestry} package.
#'   }
#' @return An object from a class that contains the \code{CATEestimator}
#'   class. It should be used with one of the following functions:
#'   \code{EstimateCATE}, \code{CateCI}, and \code{CateBIAS}. The object has at least the
#'   following slots:
#'   \item{\code{feature_train}}{A copy of feat.}
#'   \item{\code{tr_train}}{A copy of tr.}
#'   \item{\code{yobs_train}}{A copy of yobs.}
#'   \item{\code{creator}}{Function call that creates the CATE estimator. This
#'   is used for different bootstrap procedures.}
#' @author Soeren R. Kuenzel
#' @references
#' \itemize{
#'   \item Sören Künzel, Jasjeet Sekhon, Peter Bickel, and Bin Yu (2017). 
#'     MetaLearners for Estimating Heterogeneous Treatment Effects using
#'     Machine Learning. 
#'     \url{https://www.pnas.org/content/116/10/4156}
#'   \item 
#'     Sören Künzel, Simon Walter, and Jasjeet Sekhon (2018).
#'     Causaltoolbox---Estimator Stability for Heterogeneous Treatment Effects.
#'     \url{https://arxiv.org/pdf/1811.02833.pdf}
#'   \item Sören Künzel, Bradly Stadie, Nikita Vemuri, Varsha Ramakrishnan, 
#'     Jasjeet Sekhon, and Pieter Abbeel (2018). 
#'     Transfer Learning for Estimating Causal Effects using Neural Networks. 
#'     \url{https://arxiv.org/pdf/1808.07804.pdf}
#'   }
#' @family metalearners
#' @examples
#' require(causalToolbox)
#' 
#' # create example data set
#' simulated_experiment <- simulate_causal_experiment(
#'   ntrain = 1000,
#'   ntest = 1000,
#'   dim = 10
#' )
#' feat <- simulated_experiment$feat_tr
#' tr <- simulated_experiment$W_tr
#' yobs <- simulated_experiment$Yobs_tr
#' feature_test <- simulated_experiment$feat_te
#' 
#' # create the CATE estimator using Random Forests (RF)
#' xl_rf <- X_RF(feat = feat, tr = tr, yobs = yobs)
#' tl_rf <- T_RF(feat = feat, tr = tr, yobs = yobs)
#' sl_rf <- S_RF(feat = feat, tr = tr, yobs = yobs)
#' ml_rf <- M_RF(feat = feat, tr = tr, yobs = yobs)
#' xl_bt <- X_BART(feat = feat, tr = tr, yobs = yobs)
#' tl_bt <- T_BART(feat = feat, tr = tr, yobs = yobs)
#' sl_bt <- S_BART(feat = feat, tr = tr, yobs = yobs)
#' ml_bt <- M_BART(feat = feat, tr = tr, yobs = yobs)
#'   
#' cate_esti_xrf <- EstimateCate(xl_rf, feature_test)
#'
#' # evaluate the performance.
#' cate_true <- simulated_experiment$tau_te
#' mean((cate_esti_xrf - cate_true) ^ 2)
#' \dontrun{
#' # create confidence intervals via bootstrapping. 
#' xl_ci_rf <- CateCI(xl_rf, feature_test, B = 500)
#' }
#' @export 
X_RF <-
  function(feat,
           tr,
           yobs,
           predmode = "propmean",
           nthread = 0,
           verbose = FALSE,
           mu.forestry =
             list(
               relevant.Variable = 1:ncol(feat),
               ntree = 1000,
               replace = TRUE,
               sample.fraction = 0.8,
               mtry = round(ncol(feat) * 13 / 20),
               nodesizeSpl = 2,
               nodesizeAvg = 1,
               nodesizeStrictSpl = 2,
               nodesizeStrictAvg = 1,
               splitratio = 1,
               middleSplit = TRUE,
               OOBhonest = TRUE
             ),
           tau.forestry =
             list(
               relevant.Variable = 1:ncol(feat),
               ntree = 1000,
               replace = TRUE,
               sample.fraction = 0.7,
               mtry = round(ncol(feat) * 17 / 20),
               nodesizeSpl = 5,
               nodesizeAvg = 6,
               nodesizeStrictSpl = 3,
               nodesizeStrictAvg = 1,
               splitratio = 1,
               middleSplit = TRUE,
               OOBhonest = TRUE
             ),
           e.forestry =
             list(
               relevant.Variable = 1:ncol(feat),
               ntree = 500,
               replace = TRUE,
               sample.fraction =  0.5,
               mtry = ncol(feat),
               nodesizeSpl = 11,
               nodesizeAvg = 33,
               nodesizeStrictSpl = 2,
               nodesizeStrictAvg = 1,
               splitratio = 1,
               middleSplit = FALSE,
               OOBhonest = TRUE
             )) {
    
    # Cast input data to a standard format -------------------------------------
    feat <- as.data.frame(feat)
  
    # Catch misspecification erros ---------------------------------------------
    if (!(nthread - round(nthread) == 0) | nthread < 0) {
      stop("nthread must be a positive integer!")
    }
    
    if (!is.logical(verbose)) {
      stop("verbose must be either TRUE or FALSE.")
    }
    
    if (!predmode %in% c("propmean", "extreme", "control", "treated")) {
      stop("predmode should be one of propmean, extreme, control, or treated.")
    }
    
    catch_input_errors(feat, yobs, tr)
    
    # Set relevant relevant.Variable -------------------------------------------
    # User often sets the relevant variables by column names and not numerical
    # values. We translate it here to the index of the columns.
    
    if (is.null(mu.forestry$relevant.Variable)) {
      mu.forestry$relevant.Variable <- 1:ncol(feat)
    } else{
      if (is.character(mu.forestry$relevant.Variable))
        mu.forestry$relevant.Variable <-
          which(colnames(feat) %in% mu.forestry$relevant.Variable)
    }
    
    if (is.null(tau.forestry$relevant.Variable)) {
      tau.forestry$relevant.Variable <- 1:ncol(feat)
    } else{
      if (is.character(tau.forestry$relevant.Variable))
        tau.forestry$relevant.Variable <-
          which(colnames(feat) %in% tau.forestry$relevant.Variable)
    }
    
    if (is.null(e.forestry$relevant.Variable)) {
      e.forestry$relevant.Variable <- 1:ncol(feat)
    } else{
      if (is.character(e.forestry$relevant.Variable))
        e.forestry$relevant.Variable <-
          which(colnames(feat) %in% e.forestry$relevant.Variable)
    }
    
    # Translate the settings to a feature list ---------------------------------
    general_hyperpara <- list("predmode" = predmode,
                              "nthread" = nthread)
    
    hyperparameter_list <- list(
      "general" = general_hyperpara,
      "l_first_0" = mu.forestry,
      "l_first_1" = mu.forestry,
      "l_second_0" = tau.forestry,
      "l_second_1" = tau.forestry,
      "l_prop" = e.forestry
    )
    
    return(
      X_RF_fully_specified(
        feat = feat,
        tr = tr,
        yobs = yobs,
        hyperparameter_list = hyperparameter_list,
        verbose = verbose
      )
    )
  }

# X-RF basic constructor -------------------------------------------------------
X_RF_fully_specified <-
  function(feat,
           tr,
           yobs,
           hyperparameter_list,
           verbose) {
    
    # First stage --------------------------------------------------------------
    yobs_0 <- yobs[tr == 0]
    yobs_1 <- yobs[tr == 1]
    
    X_0 <- feat[tr == 0, ]
    X_1 <- feat[tr == 1, ]
    
    m_0 <-
      Rforestry::forestry(
        x = X_0[, hyperparameter_list[["l_first_0"]]$relevant.Variable],
        y = yobs_0,
        ntree = hyperparameter_list[["l_first_0"]]$ntree,
        replace = hyperparameter_list[["l_first_0"]]$replace,
        sample.fraction = hyperparameter_list[["l_first_0"]]$sample.fraction,
        mtry = hyperparameter_list[["l_first_0"]]$mtry,
        nodesizeSpl = hyperparameter_list[["l_first_0"]]$nodesizeSpl,
        nodesizeAvg = hyperparameter_list[["l_first_0"]]$nodesizeAvg,
        nodesizeStrictSpl = hyperparameter_list[["l_first_0"]]$nodesizeStrictSpl,
        nodesizeStrictAvg = hyperparameter_list[["l_first_0"]]$nodesizeStrictAvg,
        nthread = hyperparameter_list[["general"]]$nthread,
        splitrule = "variance",
        splitratio = hyperparameter_list[["l_first_0"]]$splitratio,
        OOBhonest = hyperparameter_list[["l_first_0"]]$OOBhonest
      )
    
    m_1 <-
      Rforestry::forestry(
        x = X_1[, hyperparameter_list[["l_first_1"]]$relevant.Variable],
        y = yobs_1,
        ntree = hyperparameter_list[["l_first_1"]]$ntree,
        replace = hyperparameter_list[["l_first_1"]]$replace,
        sample.fraction = hyperparameter_list[["l_first_1"]]$sample.fraction,
        mtry = hyperparameter_list[["l_first_1"]]$mtry,
        nodesizeSpl = hyperparameter_list[["l_first_1"]]$nodesizeSpl,
        nodesizeAvg = hyperparameter_list[["l_first_1"]]$nodesizeAvg,
        nodesizeStrictSpl = hyperparameter_list[["l_first_1"]]$nodesizeStrictSpl,
        nodesizeStrictAvg = hyperparameter_list[["l_first_1"]]$nodesizeStrictAvg,
        nthread = hyperparameter_list[["general"]]$nthread,
        splitrule = "variance",
        splitratio = hyperparameter_list[["l_first_1"]]$splitratio,
        OOBhonest = hyperparameter_list[["l_first_1"]]$OOBhonest
      )
    
    if (verbose) {
      print("Done with the first stage.")
    }
    
    # Second Stage -------------------------------------------------------------
    r_0 <- 
      predict(m_1, 
              X_0[, hyperparameter_list[["l_first_0"]]$relevant.Variable]) - 
      yobs_0
    r_1 <-
      yobs_1 - 
      predict(m_0, X_1[, hyperparameter_list[["l_first_1"]]$relevant.Variable])
    
    m_tau_0 <-
      Rforestry::forestry(
        x = X_0[, hyperparameter_list[["l_second_0"]]$relevant.Variable],
        y = r_0,
        ntree = hyperparameter_list[["l_second_0"]]$ntree,
        replace = hyperparameter_list[["l_second_0"]]$replace,
        sample.fraction = hyperparameter_list[["l_second_0"]]$sample.fraction,
        mtry = hyperparameter_list[["l_second_0"]]$mtry,
        nodesizeSpl = hyperparameter_list[["l_second_0"]]$nodesizeSpl,
        nodesizeAvg = hyperparameter_list[["l_second_0"]]$nodesizeAvg,
        nodesizeStrictSpl = hyperparameter_list[["l_second_0"]]$nodesizeStrictSpl,
        nodesizeStrictAvg = hyperparameter_list[["l_second_0"]]$nodesizeStrictAvg,
        nthread = hyperparameter_list[["general"]]$nthread,
        splitrule = "variance",
        splitratio = hyperparameter_list[["l_second_0"]]$splitratio,
        OOBhonest = hyperparameter_list[["l_second_0"]]$OOBhonest
      )
    
    m_tau_1 <-
      Rforestry::forestry(
        x = X_1[, hyperparameter_list[["l_second_1"]]$relevant.Variable],
        y = r_1,
        ntree = hyperparameter_list[["l_second_1"]]$ntree,
        replace = hyperparameter_list[["l_second_1"]]$replace,
        sample.fraction = hyperparameter_list[["l_second_1"]]$sample.fraction,
        mtry = hyperparameter_list[["l_second_1"]]$mtry,
        nodesizeSpl = hyperparameter_list[["l_second_1"]]$nodesizeSpl,
        nodesizeAvg = hyperparameter_list[["l_second_1"]]$nodesizeAvg,
        nodesizeStrictSpl = hyperparameter_list[["l_second_1"]]$nodesizeStrictSpl,
        nodesizeStrictAvg = hyperparameter_list[["l_second_1"]]$nodesizeStrictAvg,
        nthread = hyperparameter_list[["general"]]$nthread,
        splitrule = "variance",
        splitratio = hyperparameter_list[["l_second_1"]]$splitratio,
        OOBhonest = hyperparameter_list[["l_second_1"]]$OOBhonest
      )
    if (verbose) {
      print("Done with the second stage.")
    }
    
    # Prop score estimation ----------------------------------------------------
    m_prop <-
      Rforestry::forestry(
        x = feat[, hyperparameter_list[["l_prop"]]$relevant.Variable],
        y = tr,
        ntree = hyperparameter_list[["l_prop"]]$ntree,
        replace = hyperparameter_list[["l_prop"]]$replace,
        sample.fraction = hyperparameter_list[["l_prop"]]$sample.fraction,
        mtry = hyperparameter_list[["l_prop"]]$mtry,
        nodesizeSpl = hyperparameter_list[["l_prop"]]$nodesizeSpl,
        nodesizeAvg = hyperparameter_list[["l_prop"]]$nodesizeAvg,
        nodesizeStrictSpl = hyperparameter_list[["l_prop"]]$nodesizeStrictSpl,
        nodesizeStrictAvg = hyperparameter_list[["l_prop"]]$nodesizeStrictAvg,
        nthread = hyperparameter_list[["general"]]$nthread,
        splitrule = "variance",
        splitratio = hyperparameter_list[["l_prop"]]$splitratio,
        OOBhonest = hyperparameter_list[["l_prop"]]$OOBhonest
      )
    if (verbose) {
      print("Done with the propensity score estimation.")
    }
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
          X_RF_fully_specified(feat,
                               tr,
                               yobs,
                               hyperparameter_list,
                               verbose)
        }
      )
    )
  }

# Estimate CATE Method ---------------------------------------------------------
#' EstimateCate-X_hRF
#' @rdname EstimateCate
#' @inherit EstimateCate
#' @exportMethod EstimateCate
setMethod(
  f = "EstimateCate",
  signature = "X_RF",
  definition = function(theObject, feature_new)
  {
    feature_new <- as.data.frame(feature_new)

    catch_feat_input_errors(feature_new)
    
    predmode <- theObject@hyperparameter_list[["general"]]$predmode
    prop_scores <- predict(theObject@m_prop, feature_new)
    
    
    
    if (predmode == "propmean") {
      return(
        prop_scores        * predict(theObject@m_tau_0, feature_new) +
          (1 - prop_scores)  * predict(theObject@m_tau_1, feature_new)
      )
    }
    
    if (predmode == "extreme") {
      return(ifelse(
        prop_scores > .5,
        predict(theObject@m_tau_0, feature_new),
        predict(theObject@m_tau_1, feature_new)
      ))
    }
    
    if (predmode == "control") {
      return(predict(theObject@m_tau_0, feature_new))
    }
    
    if (predmode == "treated") {
      return(predict(theObject@m_tau_1, feature_new))
    }
    
    stop("predmode should be one of propmean, extreme, control, or treated.")
    
  }
)



