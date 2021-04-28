#' @include CATE_estimators.R
#' @include helper_functions.R
NULL

# X-RF class -------------------------------------------------------------------
setClass(
  "DR_RF",
  contains = "MetaLearner",
  slots = list(
    pseudo_out = "forestry",
    m_tau_0 = "forestry",
    m_tau_1 = "forestry",
    m_prop = "forestry",
    hyperparameter_list = "list"
  )
)


# DR_RF generator ---------------------------------------------------------------
#' @title DR-Learners
#' @rdname DRleaners
#' @name DR-Learner
#' @description DR_RF is an implementation of the DR-learner with Random Forests
#'   (Breiman 2001) as the base learners.
#' @param feat A data frame containing the features.
#' @param tr A numeric vector with 0 for control and 1 for treated variables.
#' @param yobs A numeric vector containing the observed outcomes.
#' @param predmode Specifies how the two estimators of the second stage should
#'   be aggregated. Possible types are "propmean," "control," and "treated." The
#'   default is "propmean," which refers to propensity score weighting.
#' @param trunc_level Level at which to truncate the estimated propensity scores
#'   this ensures that the predicted propensity scores are bounded between
#'   trunc_level < p_score < 1-trunc_level. Default is .02.
#' @param nthread Number of threads which should be used to work in parallel.
#' @param verbose TRUE for detailed output, FALSE for no output.
#' @param prop.forestry,tau.forestry,mu.forestry,pseu.forestry A list containing the
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
#'   \item Edward Kennedy (2020).
#'     Optimal doubly robust estimation of heterogeneous causal effects.
#'     \url{https://arxiv.org/abs/2004.14497}
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
DR_RF <-
  function(feat,
           tr,
           yobs,
           predmode = "propmean",
           nthread = 0,
           verbose = FALSE,
           trunc_level = .02,
           prop.forestry =
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
           mu.forestry =
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
           pseu.forestry =
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
             ))
{

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

    if (is.null(prop.forestry$relevant.Variable)) {
      prop.forestry$relevant.Variable <- 1:ncol(feat)
    } else{
      if (is.character(prop.forestry$relevant.Variable))
        prop.forestry$relevant.Variable <-
          which(colnames(feat) %in% prop.forestry$relevant.Variable)
    }

    if (is.null(pseu.forestry$relevant.Variable)) {
      pseu.forestry$relevant.Variable <- 1:ncol(feat)
    } else{
      if (is.character(pseu.forestry$relevant.Variable))
        pseu.forestry$relevant.Variable <-
          which(colnames(feat) %in% pseu.forestry$relevant.Variable)
    }

    # Translate the settings to a feature list ---------------------------------
    general_hyperpara <- list("predmode" = predmode,
                              "nthread" = nthread,
                              "trunc_level" = trunc_level)

    hyperparameter_list <- list(
      "general" = general_hyperpara,
      "l_first_0" = mu.forestry,
      "l_first_1" = tau.forestry,
      "l_pseudo" = pseu.forestry,
      "l_prop" = prop.forestry
    )

    return(
      DR_RF_fully_specified(
        feat = feat,
        tr = tr,
        yobs = yobs,
        hyperparameter_list = hyperparameter_list,
        verbose = verbose
      )
    )
  }


# X-RF basic constructor -------------------------------------------------------
DR_RF_fully_specified <-
  function(feat,
           tr,
           yobs,
           hyperparameter_list,
           verbose) {

    if (verbose) {
      print("Running DR Learner- RF base learners")
    }

    trunc_level <- hyperparameter_list[["general"]][["trunc_level"]]
    # First do data splitting
    initial_split <-  caret::createFolds(yobs, k = 3)

    # Nuisance training

    # Propensity Score regression using fold 1
    if (verbose) {
      print("Fitting Propensity Score Regression")
    }


    propensity_reg <- forestry(x = feature_train[initial_split$Fold1,
                                                 hyperparameter_list[["l_prop"]]$relevant.Variable],
                               y = w_train[initial_split$Fold1],
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
                               OOBhonest = hyperparameter_list[["l_prop"]]$OOBhonest)

    # Outcome Regression using fold 2
    outcome_reg_data <- feature_train[initial_split$Fold2,]
    outcome_reg_outcome <- yobs_train[initial_split$Fold2]

    treat_idx <- which(w_train[initial_split$Fold2] == 1)
    contr_idx <- which(w_train[initial_split$Fold2] == 0)

    if (verbose) {
      print("Fitting Treatment Outcome Regression")
    }

    treat_reg <- forestry(x = outcome_reg_data[treat_idx,
                                               hyperparameter_list[["l_first_1"]]$relevant.Variable],
                          y = outcome_reg_outcome[treat_idx],
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
                          OOBhonest = hyperparameter_list[["l_first_1"]]$OOBhonest)

    if (verbose) {
      print("Fitting Control Outcome Regression")
    }

    contr_reg <- forestry(x = outcome_reg_data[contr_idx,
                                               hyperparameter_list[["l_first_0"]]$relevant.Variable],
                          y = outcome_reg_outcome[contr_idx],
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
                          OOBhonest = hyperparameter_list[["l_first_0"]]$OOBhonest)

    # Phase 3 Pseudo outcome regression
    pseudo_reg_data <- feature_train[initial_split$Fold3,]
    pseudo_out <- yobs_train[initial_split$Fold3]
    pseudo_ind <- w_train[initial_split$Fold3]


    reg_p_score <- predict(propensity_reg, feature.new = pseudo_reg_data)

    # Make sure we truncate the estimated propensity score for stability
    data.frame(X1 = reg_p_score) %>%
      dplyr::mutate(X1 = case_when(X1 < trunc_level ~ trunc_level,
                                   X1 > 1-trunc_level ~ 1-trunc_level,
                                   TRUE ~ X1)) %>%
      select(X1) -> truncated_reg_pscore

    trunc_p <- truncated_reg_pscore[,1]
    c_pred <- predict(contr_reg, feature.new = pseudo_reg_data)
    t_pred <- predict(treat_reg, feature.new = pseudo_reg_data)
    tailored_pred <- ifelse(pseudo_ind,
                            t_pred,
                            c_pred)

    # Now create the pseudo outcome
    if (verbose) {
      print("Fitting Pseudo Outcome Regression")
    }

    pseudo_outcome <- ((pseudo_ind - trunc_p) / (trunc_p*(1-trunc_p)))*(pseudo_out - tailored_pred) + t_pred - c_pred
    pseudo_outcome_reg <- forestry(x = pseudo_reg_data,
                                   y = pseudo_outcome,
                                   ntree = hyperparameter_list[["l_pseudo"]]$ntree,
                                   replace = hyperparameter_list[["l_pseudo"]]$replace,
                                   sample.fraction = hyperparameter_list[["l_pseudo"]]$sample.fraction,
                                   mtry = hyperparameter_list[["l_pseudo"]]$mtry,
                                   nodesizeSpl = hyperparameter_list[["l_pseudo"]]$nodesizeSpl,
                                   nodesizeAvg = hyperparameter_list[["l_pseudo"]]$nodesizeAvg,
                                   nodesizeStrictSpl = hyperparameter_list[["l_pseudo"]]$nodesizeStrictSpl,
                                   nodesizeStrictAvg = hyperparameter_list[["l_pseudo"]]$nodesizeStrictAvg,
                                   nthread = hyperparameter_list[["general"]]$nthread,
                                   splitrule = "variance",
                                   splitratio = hyperparameter_list[["l_pseudo"]]$splitratio,
                                   OOBhonest = hyperparameter_list[["l_pseudo"]]$OOBhonest)

    return(
      new(
        "DR_RF",
        feature_train = feat,
        tr_train = tr,
        yobs_train = yobs,
        pseudo_out = pseudo_outcome_reg,
        m_tau_0 = contr_reg,
        m_tau_1 = treat_reg,
        m_prop = propensity_reg,
        hyperparameter_list = hyperparameter_list,
        creator = function(feat, tr, yobs) {
          DR_RF_fully_specified(feat,
                                tr,
                                yobs,
                                hyperparameter_list,
                                verbose)
        }
      )
    )
  }

# Estimate CATE Method ---------------------------------------------------------
#' EstimateCate-DR_hRF
#' @rdname EstimateCate-DR_hRF
#' @inherit EstimateCate
#' @exportMethod EstimateCate
setMethod(
  f = "EstimateCate",
  signature = "DR_RF",
  definition = function(theObject, feature_new)
  {
    feature_new <- as.data.frame(feature_new)

    catch_feat_input_errors(feature_new)

    predmode <- theObject@hyperparameter_list[["general"]]$predmode
    prop_scores <- predict(theObject@m_prop, feature_new)



    if (predmode == "propmean") {
      return(
        prop_scores * predict(theObject@m_tau_0, feature_new) +
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
