#' @include CATE_estimators.R

# M-RF class -------------------------------------------------------------------
setClass(
  "M_RF",
  contains = "MetaLearner",
  slots = list(
    m_0 = "forestry",
    m_1 = "forestry",
    m_prop = "forestry",
    m_tau = "forestry",
    hyperparameter_list = "list"
  ),
  validity = function(object)
  {
    if (!all(object@tr_train %in% c(0, 1))) {
      return("TR is the treatment and must be either 0 or 1")
    }
    return(TRUE)
  }
)

#' @title M-Learners
#' @rdname Mlearners
#' @name M-Learner
#' @description \code{M_RF} is an implementation of the Modified Outcome Estimator with
#'   Random Forest (Breiman 2001) as the base learner.
#' @details 
#' The M-Learner estimates the CATE in two steps:
#' \enumerate{
#'  \item
#'     Estimate the response functions and the propensity score,
#'     \deqn{\mu_0(x) = E[Y(0) | X = x]}
#'     \deqn{\mu_1(x) = E[Y(1) | X = x]} 
#'     \deqn{e(x) = E[W | X = x]} 
#'     using the base learner and denote the estimates as \eqn{\hat \mu_0},
#'     \eqn{\hat \mu_1}, and \eqn{\hat e}.
#'  \item
#'     Define the adjusted modified outcomes as 
#'     \deqn{ R _i = (Z_i - \hat e(x_i)) / (\hat e(x_i)[1 - \hat e(x_i)]) 
#'     (Y_i - \hat \mu_1(x_i) [1 - \hat e(x_i)] - \hat \mu_0(x_i)\hat e(x_i)).}
#'     Now employ the base learner to estimate
#'     \deqn{\tau(x) = E[R | X = x].} 
#'     The result is the CATE estimator.
#'     }
#' @param mu.forestry,tau.forestry,e.forestry A list containing the
#'   hyperparameters for the \code{forestry} package that are used for
#'   estimating the response functions, the CATE, and the propensity score.
#'   These hyperparameters are passed to the \code{forestry} package. (Please
#'   refer to the \href{https://github.com/soerenkuenzel/forestry}{forestry}
#'   package for a more detailed documentation of the hyperparamters.)
#'   \itemize{
#'      \item \code{relevant.Variable} Variables that are only used in the first 
#'            stage.
#'      \item \code{ntree} Numbers of trees used in the first stage.
#'      \item \code{replace} Sample with or without replacement in the first 
#'            stage.
#'      \item \code{sample.fraction} Size of total samples drawn for the 
#'            training data in the first stage.
#'      \item \code{mtry} Number of variables randomly selected in each 
#'            splitting point.
#'      \item \code{nodesizeSpl} Minimum nodesize in the first stage for 
#'            the observations in the splitting set. (See the details of the 
#'            \code{forestry} package)
#'      \item \code{nodesizeAvg} Minimum nodesize in the first stage for 
#'            the observations in the averaging set.
#'      \item \code{splitratio} Proportion of the training data used as the 
#'            splitting dataset in the first stage.
#'      \item \code{middleSplit} If true, the split value will be exactly in the 
#'            middle of two observations. Otherwise, it will take a point 
#'            based on a uniform distribution between the two observations. 
#'   }
#' @references
#' \itemize{
#'   \item Sören Künzel, Jasjeet Sekhon, Peter Bickel, and Bin Yu (2017). 
#'     MetaLearners for Estimating Heterogeneous Treatment Effects Using
#'     Machine Learning. 
#'     \url{https://www.pnas.org/content/116/10/4156}
#'   \item 
#'     Sören Künzel, Simon Walter, and Jasjeet Sekhon (2018).
#'     Causaltoolbox---Estimator Stability for Heterogeneous Treatment Effects.
#'     \url{https://arxiv.org/pdf/1811.02833.pdf}
#'   \item Daniel Rubin and Mark J van der Laan (2007). A Doubly Robust
#'   Censoring Unbiased Transformation.
#'     \url{https://www.ncbi.nlm.nih.gov/pubmed/22550646}
#'   }
#' @inherit X-Learner
#' @family metalearners
#' @export
M_RF <-
  function(feat,
           tr,
           yobs,
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
               splitratio = 1,
               middleSplit = TRUE
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
               splitratio = .5,
               middleSplit = FALSE
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
               splitratio = 0.8,
               middleSplit = TRUE
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
    general_hyperpara <- list("nthread" = nthread)
    
    hyperparameter_list <- list(
      "general" = general_hyperpara,
      "l_first_0" = mu.forestry,
      "l_first_1" = mu.forestry,
      "l_second" = tau.forestry,
      "l_prop" = e.forestry
    )
    
    return(
      M_RF_fully_specified(
        feat = feat,
        tr = tr,
        yobs = yobs,
        hyperparameter_list = hyperparameter_list,
        verbose = verbose
      )
    )
  }
    
# M-RF basic constructor -------------------------------------------------------
M_RF_fully_specified <-
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
        nthread = hyperparameter_list[["general"]]$nthread,
        splitrule = "variance",
        splitratio = hyperparameter_list[["l_first_0"]]$splitratio
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
        nthread = hyperparameter_list[["general"]]$nthread,
        splitrule = "variance",
        splitratio = hyperparameter_list[["l_first_1"]]$splitratio
      )
    
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
        nthread = hyperparameter_list[["general"]]$nthread,
        splitrule = "variance",
        splitratio = hyperparameter_list[["l_prop"]]$splitratio
      )
    
    if (verbose) {
      print("Done with the first stage.")
    }
    
    # Second Stage -------------------------------------------------------------

      m0_hat <- predict(m_0, feat)
      m1_hat <- predict(m_1, feat)

      propensity_score_hat <- predict(m_prop, feat)
      propensity_score_hat <- pmin(pmax(propensity_score_hat, 0.01), 0.99)

      modified_outcome_ra <-
        (tr - propensity_score_hat) / 
        (propensity_score_hat * (1 - propensity_score_hat)) *
        (yobs - m1_hat * (1 - propensity_score_hat) - 
           m0_hat * propensity_score_hat)

      m_tau <- Rforestry::forestry(
        x = feat,
        y = modified_outcome_ra,
        ntree = hyperparameter_list[["l_second"]]$ntree,
        replace = hyperparameter_list[["l_second"]]$replace,
        sample.fraction = hyperparameter_list[["l_second"]]$sample.fraction,
        mtry = hyperparameter_list[["l_second"]]$mtry,
        nodesizeSpl = hyperparameter_list[["l_second"]]$nodesizeSpl,
        nodesizeAvg = hyperparameter_list[["l_second"]]$nodesizeAvg,
        nthread = hyperparameter_list[["general"]]$nthread,
        splitrule = "variance",
        splitratio = hyperparameter_list[["l_second"]]$splitratio
      )

    new(
      "M_RF",
      feature_train = feat,
      tr_train = tr,
      yobs_train = yobs,
      m_0 = m_0,
      m_1 = m_1,
      m_prop = m_prop,
      m_tau = m_tau,
      hyperparameter_list = hyperparameter_list,
      creator = function(feat, tr, yobs) {
        M_RF_fully_specified(feat,
                             tr,
                             yobs,
                             hyperparameter_list = hyperparameter_list,
                             verbose = verbose)
      }
    )
  }


#' EstimateCate-S_hRF
#' @rdname EstimateCate-S_hRF
#' @export
setMethod(
  f = "EstimateCate",
  signature = "M_RF",
  definition = function(theObject, feature_new)
  {
    feature_new <- as.data.frame(feature_new)
    return(predict(theObject@m_tau, feature_new))
  }
)
