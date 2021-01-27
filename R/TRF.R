# This file implements the T-Learner (https://arxiv.org/pdf/1706.03461.pdf)
# with the forestry implementation (https://github.com/soerenkuenzel/forestry)
# as base learner.
#' @include CATE_estimators.R
#' @include helper_functions.R
#' @include XRF.R
NULL

# T-RF class -------------------------------------------------------------------
setClass(
  "T_RF",
  contains = "MetaLearner",
  slots = list(
    feature_train = "data.frame",
    tr_train = "numeric",
    yobs_train = "numeric",
    m_0 = "forestry",
    m_1 = "forestry",
    hyperparameter_list = "list", 
    creator = "function"
  ),
  validity = function(object)
  {
    if (!all(object@tr_train %in% c(0, 1))) {
      return("TR is the treatment and must be either 0 or 1")
    }
    return(TRUE)
  }
)

# T_RF generator ---------------------------------------------------------------
#' @title T-Learners
#' @rdname Tlearners
#' @name T-Learner
#' @description \code{T_RF} is an implementation of the T-learner combined with Random
#'   Forest (Breiman 2001) for both response functions.
#' @details 
#' The CATE is estimated using two estimators:
#' \enumerate{
#'  \item
#'     Estimate the response functions 
#'     \deqn{\mu_0(x) = E[Y(0) | X = x]}
#'     \deqn{\mu_1(x) = E[Y(1) | X = x]} 
#'     using the base leaner and denote the estimates as \eqn{\hat \mu_0} and
#'     \eqn{\hat \mu_1}.
#'  \item 
#'     Define the CATE estimate as
#'     \deqn{\tau(x) = \hat \mu_1 - \hat \mu_0.}
#' }
#' @param mu0.forestry,mu1.forestry Lists containing the hyperparameters for the
#'   \code{forestry} package that are used in \eqn{\hat \mu_0} and \eqn{\hat
#'   \mu_1}, respectively. These hyperparameters are passed to the
#'   \code{forestry} package. (Please refer to the
#'   \href{https://github.com/soerenkuenzel/forestry}{\code{forestry}} package
#'   for a more detailed documentation of the hyperparamters.)
#' @return Object of class \code{T_RF}. It should be used with one of the
#'   following functions \code{EstimateCATE}, \code{CateCI}, and
#'   \code{CateBIAS}. The object has the following slots:
#'   \item{\code{feature_train}}{A copy of feat.}
#'   \item{\code{tr_train}}{A copy of tr.}
#'   \item{\code{yobs_train}}{A copy of yobs.}
#'   \item{\code{m_0}}{An object of class forestry that is fitted with the 
#'      observed outcomes of the control group as the dependent variable.}
#'   \item{\code{m_1}}{An object of class forestry that is fitted with the 
#'      observed outcomes of the treated group as the dependent variable.}
#'   \item{\code{hyperparameter_list}}{List containting the hyperparameters of 
#'      the three random forest algorithms used.}
#'   \item{\code{creator}}{Function call of T_RF. This is used for different 
#'      bootstrap procedures.}
#' @inherit X-Learner
#' @family metalearners
#' @export
T_RF <-
  function(feat,
           tr,
           yobs, 
           nthread = 0,
           verbose = TRUE,
           mu0.forestry =
             list(
               relevant.Variable = 1:ncol(feat),
               ntree = 1000,
               replace = TRUE,
               sample.fraction = 0.9,
               mtry = ncol(feat),
               nodesizeSpl = 1,
               nodesizeAvg = 3,
               splitratio = .5,
               middleSplit = FALSE
             ),
           mu1.forestry =
             list(
               relevant.Variable = 1:ncol(feat),
               ntree = 1000,
               replace = TRUE,
               sample.fraction = 0.9,
               mtry = ncol(feat),
               nodesizeSpl = 1,
               nodesizeAvg = 3,
               splitratio = .5,
               middleSplit = FALSE
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
    
    if (is.null(mu0.forestry$relevant.Variable)) {
      mu0.forestry$relevant.Variable <- 1:ncol(feat)
    } else{
      if (is.character(mu0.forestry$relevant.Variable))
        mu0.forestry$relevant.Variable <-
          which(colnames(feat) %in% mu0.forestry$relevant.Variable)
    }
    
    if (is.null(mu1.forestry$relevant.Variable)) {
      mu1.forestry$relevant.Variable <- 1:ncol(feat)
    } else{
      if (is.character(mu1.forestry$relevant.Variable))
        mu1.forestry$relevant.Variable <-
          which(colnames(feat) %in% mu1.forestry$relevant.Variable)
    }
    
    # Translate the settings to a feature list ---------------------------------
    general_hyperpara <- list("nthread" = nthread)
    
    hyperparameter_list <- list(
      "general" = general_hyperpara,
      "l_first_0" = mu0.forestry,
      "l_first_1" = mu1.forestry
    )
    
    return(
     T_RF_fully_specified(
        feat = feat,
        tr = tr,
        yobs = yobs,
        hyperparameter_list = hyperparameter_list,
        verbose = verbose
      )
    )
  }

# T-RF basic constructor -------------------------------------------------------
T_RF_fully_specified <-
  function(feat,
           tr,
           yobs,
           hyperparameter_list,
           verbose) {
    yobs_0 <- yobs[tr == 0]
    yobs_1 <- yobs[tr == 1]
    
    X_0 <- feat[tr == 0,]
    X_1 <- feat[tr == 1,]
    
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
    
    
    new(
      "T_RF",
      feature_train = feat,
      tr_train = tr,
      yobs_train = yobs,
      m_0 = m_0,
      m_1 = m_1,
      hyperparameter_list = hyperparameter_list,
      creator = function(feat, tr, yobs) {
        T_RF_fully_specified(feat,
                             tr,
                             yobs,
                             hyperparameter_list,
                             verbose)
      }
    )
  }

# Estimate CATE Method ---------------------------------------------------------
#' EstimateCate-T_hRF
#' @rdname EstimateCate
#' @inherit EstimateCate
#' @exportMethod EstimateCate
setMethod(
  f = "EstimateCate",
  signature = "T_RF",
  definition = function(theObject, feature_new)
  {
    feature_new <- as.data.frame(feature_new)
    catch_feat_input_errors(feature_new)

    return(
      predict(theObject@m_1, feature_new) -
        predict(theObject@m_0, feature_new)
    )
  }
)

