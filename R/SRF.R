#' @include CATE_estimators.R

# S-RF class -------------------------------------------------------------------
setClass(
  "S_RF",
  contains = "MetaLearner",
  slots = list(
    feature_train = "data.frame",
    tr_train = "numeric",
    yobs_train = "numeric",
    forest = "forestry",
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

# S_RF generator ---------------------------------------------------------------
#' @title S-Learners
#' @description \code{S_RF} is an implementation of the S-Learner combined with
#'   Random Forests (Breiman 2001).
#' @rdname Slearners
#' @name S-Learner
#' @details 
#' In the S-Learner, the outcome is estimated using all of the features and the 
#' treatment indicator without giving the treatment indicator a special role. 
#' The predicted CATE for an individual unit is then the difference between the 
#' predicted values when the treatment assignment indicator is changed from 
#' control to treatment:
#' \enumerate{
#'  \item
#'     Estimate the joint response function \deqn{\mu(x, w) = E[Y | X = x, W =
#'     w]} using the
#'     base learner.
#'     We denote the estimate as \eqn{\hat \mu}.
#'  \item 
#'     Define the CATE estimate as
#'     \deqn{\tau(x) = \hat \mu_1(x, 1) - \hat \mu_0(x, 0).}
#' }
#' @param mu.forestry A list containing the hyperparameters for the
#'   \code{forestry} package that are used in \eqn{\hat \mu_0}.
#'   These hyperparameters are passed to the \code{forestry} package. 
#' @return Object of class \code{S_RF}. It should be used with one of the
#'   following functions \code{EstimateCATE}, \code{CateCI}, and
#'   \code{CateBIAS}. The object has the following slots:
#'   \item{\code{feature_train}}{A copy of feat.}
#'   \item{\code{tr_train}}{A copy of tr.}
#'   \item{\code{yobs_train}}{A copy of yobs.}
#'   \item{\code{m_0}}{An object of class forestry that is fitted with the 
#'      observed outcomes of the control group as the dependent variable.}
#'   \item{\code{m_1}}{An object of class forestry that is fitted with the 
#'      observed outcomes of the treated group as the dependent variable.}
#'   \item{\code{hyperparameter_list}}{A list containting the hyperparameters of 
#'      the three random forest algorithms used.}
#'   \item{\code{creator}}{Function call of \code{S_RF}. This is used for different 
#'      bootstrap procedures.}
#' @inherit X-Learner
#' @family metalearners
#' @export
S_RF <-
  function(feat,
           tr,
           yobs, 
           nthread = 0,
           verbose = TRUE,
           mu.forestry =
             list(
               relevant.Variable = 1:ncol(feat),
               ntree = 1000,
               replace = TRUE,
               sample.fraction = 0.9,
               mtry = ncol(feat),
               nodesizeSpl = 1,
               nodesizeAvg = 3,
               nodesizeStrictSpl = 3,
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
    
    # Translate the settings to a feature list ---------------------------------
    general_hyperpara <- list("nthread" = nthread)
    
    hyperparameter_list <- list(
      "general" = general_hyperpara,
      "mu.forestry" = mu.forestry
    )
    
    return(
      S_RF_fully_specified(
        feat = feat,
        tr = tr,
        yobs = yobs,
        hyperparameter_list = hyperparameter_list,
        verbose = verbose
      )
    )
  }
    
# S-RF basic constructor -------------------------------------------------------
S_RF_fully_specified <-
  function(feat,
           tr,
           yobs,
           hyperparameter_list,
           verbose) {
    
    m <- Rforestry::forestry(
      x = cbind(feat[, hyperparameter_list[["mu.forestry"]]$relevant.Variable], 
                tr),
      y = yobs,
      ntree = hyperparameter_list[["mu.forestry"]]$ntree,
      replace = hyperparameter_list[["mu.forestry"]]$replace,
      sample.fraction = hyperparameter_list[["mu.forestry"]]$sample.fraction,
      mtry = hyperparameter_list[["mu.forestry"]]$mtry,
      nodesizeSpl = hyperparameter_list[["mu.forestry"]]$nodesizeSpl,
      nodesizeAvg = hyperparameter_list[["mu.forestry"]]$nodesizeAvg,
      nodesizeStrictSpl = hyperparameter_list[["mu.forestry"]]$nodesizeStrictSpl,
      nodesizeStrictAvg = hyperparameter_list[["mu.forestry"]]$nodesizeStrictAvg,
      nthread = hyperparameter_list[["general"]]$nthread,
      splitrule = "variance",
      splitratio = hyperparameter_list[["mu.forestry"]]$splitratio,
      OOBhonest = hyperparameter_list[["mu.forestry"]]$OOBhonest
    )
    
    new(
      "S_RF",
      feature_train = feat,
      tr_train = tr,
      yobs_train = yobs,
      forest = m,
      hyperparameter_list = hyperparameter_list,
      creator = function(feat, tr, yobs) {
        S_RF_fully_specified(feat,
                             tr,
                             yobs,
                             hyperparameter_list,
                             verbose)
      }
    )
  }

############################
### Estimate CATE Method ###
############################
#' EstimateCate-S_hRF
#' @rdname EstimateCate-S_hRF
#' @inherit EstimateCate
#' @exportMethod EstimateCate
setMethod(
  f = "EstimateCate",
  signature = "S_RF",
  definition = function(theObject, 
                        feature_new,
                        ...)
  {
    feature_new <- as.data.frame(feature_new)
    catch_feat_input_errors(feature_new)

    # Check if we want to do bias correction predictions
    if ("aggregation" %in% ls() && substr(aggregation,1,2) == "bc") {
      if (aggregation == "bc1") {
        return(
          correctedPredict(theObject@forest, cbind(feature_new, tr = 1)) -
            correctedPredict(theObject@forest, cbind(feature_new, tr = 0))
        )
      } else if (aggregation == "bc2") {
        return(
          correctedPredict(theObject@forest, cbind(feature_new, tr = 1), nrounds = 1) -
            correctedPredict(theObject@forest, cbind(feature_new, tr = 0), nrounds = 1)
        )
      } else if (aggregation == "bc3") {
        return(
          correctedPredict(theObject@forest, cbind(feature_new, tr = 1), nrounds = 1, monotone = TRUE) -
            correctedPredict(theObject@forest, cbind(feature_new, tr = 0), nrounds = 1, monotone = TRUE)
        )
      } else if (aggregation == "bc4") {
        return(
          correctedPredict(theObject@forest, cbind(feature_new, tr = 1), simple=FALSE) -
            correctedPredict(theObject@forest, cbind(feature_new, tr = 0), simple=FALSE)
        )
      } else {
        stop(paste0("Aggregation not found: ",aggregation))
      }
    }
    
    return(
      predict(theObject@forest, cbind(feature_new, tr = 1), ...) -
        predict(theObject@forest, cbind(feature_new, tr = 0), ...)
    )
  }
)
