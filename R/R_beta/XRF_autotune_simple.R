#' @include XRF.R


#' @title Simple Autotuning for X-Learner with RF
#' @name X_RF_autotune_simple
#' @rdname X_RF_autotune_simple
#' @description This function tunes the X-Learner withrandom forest by
#' testing which of 11 prespecified settings seems to be the best
#' @param feat A data frame of all the features.
#' @param tr A numeric vector containing 0 for control and 1 for treated 
#' variables.
#' @param yobs A numeric vector containing the observed outcomes.
#' @param ntree Number of trees used
#' @param nthread Number of threads which can run in parallel. If set 0, then
#' the maximum amount of possible threads is determined automatically. If set to
#' 1 then the algorithm is absolutely deterministic (after specifying a seed).
#' @param ntree_testing TODO: Add Description
#' @param verbose if tuning process in verbose mode
#' @seealso \code{\link{X_RF_autotune_gpp}},
#' \code{\link{X_RF_autotune_hyperband}}
#' @examples
#'   set.seed(14236142)
#'   feat <- iris[, -1]
#'   tr <- rbinom(nrow(iris), 1, .5)
#'   yobs <- iris[, 1]
#'   # train a
#'   xl_gpp <- X_RF_autotune_simple(feat, tr, yobs, ntree = 100, nthread = 0,
#'   verbose = FALSE)
#'   # computes the CATE and confidence intervals for CATE
#'   EstimateCate(xl_gpp, feat)
#'   CateCI(xl_gpp, feat, B = 5, verbose = FALSE)
#'
#' @export X_RF_autotune_simple
X_RF_autotune_simple <-
  function(feat,
           tr,
           yobs,
           ntree = 2000,
           ntree_testing = 600,
           nthread = 0,
           verbose = TRUE) {
    starting_settings <- get_starting_settings(feat = feat, tr = tr,
                    ntree = ntree_testing, nthread = nthread)
    #   list(
    #   "start_setting_1" = get_setting_strong(feat, ntree, nthread),
    #   "start_setting_2" = get_setting_weak(feat, ntree, nthread)
    # )

    setup_eval <-
      check_setups(starting_settings = starting_settings,
                   feat = feat,
                   tr = tr,
                   yobs = yobs,
                   ntree = ntree,
                   nthread = nthread,
                   verbose = verbose)
    if (verbose)
      print(paste(names(starting_settings)[which.min(setup_eval$comb)],
                  "is the best."))

    best_setting <- starting_settings[[which.min(setup_eval$comb)]]

    best_setting$l_first_0$ntree <- ntree
    best_setting$l_first_1$ntree <- ntree
    best_setting$l_second_0$ntree <- ntree
    best_setting$l_second_1$ntree <- ntree
    best_setting$l_prop$ntree <- ntree

    return(
      X_RF_fully_specified(
        feat = feat,
        tr = tr,
        yobs = yobs,
        hyperparameter_list = best_setting,
        verbose = verbose
      )
    )
  }

# This function checks the starting_settings using the OOB errors from those
# It returns a table containing the expected performance of each of the setups
check_setups <-
  function(starting_settings,
           feat,
           tr,
           yobs,
           ntree = 2000,
           nthread = 0,
           verbose = TRUE) {
    # define starting points 1 and 2

    OOB_errors <- data.frame()
    for (i in 1:length(starting_settings)) {
      if (verbose)
        print(paste0("Checking setup ", i, " out of ", length(starting_settings)))
      starting_setting <- starting_settings[[i]]
      starting_setting_name <- names(starting_settings)[i]
      OOB_errors <-
        tryCatch({
          rbind(OOB_errors,
                evaluate_setting(starting_setting, feat, tr, yobs))
        },
        error = function(e) {
          warning(
            paste(
              "Setting",
              starting_setting_name,
              "cannot be checked,",
              "since an error occured:"
            )
          )
          print(e)
          OOB_errors
        })
    }
    OOB_errors$tau_0 <- OOB_errors$l_first_1 + OOB_errors$l_second_0
    OOB_errors$tau_1 <- OOB_errors$l_first_0 + OOB_errors$l_second_1
    OOB_errors$comb <- best_MSE_constant(OOB_errors$tau_0, OOB_errors$tau_1)
      # (OOB_errors$tau_0 + OOB_errors$tau_1) / 2
    OOB_errors <-
      cbind(data.frame("name" = names(starting_settings)),
            OOB_errors)

    return(OOB_errors)
  }

# This function evaluates a setting by and returns the OOB error of the four
# base learners.
evaluate_setting <- function(setting, feat, tr, yobs) {
  # setting <- starting_setting
  x_eval <- X_RF_fully_specified(feat,
                                 tr,
                                 yobs,
                                 hyperparameter_list = setting,
                                 verbose = FALSE)

  return(
    data.frame(
      "l_first_0" = forestry::getOOB(x_eval@m_0, noWarning = FALSE),
      "l_first_1" = forestry::getOOB(x_eval@m_1, noWarning = FALSE),
      "l_second_0" = forestry::getOOB(x_eval@m_tau_0, noWarning = FALSE),
      "l_second_1" = forestry::getOOB(x_eval@m_tau_1, noWarning = FALSE)
    )
  )
}

################################################################################
################################################################################
# Define serveral standard settings:

get_setting_strong <- function(feat,
                               ntree,
                               nthread,
                               relevant_Variable_first = 1:ncol(feat),
                               relevant_Variable_second = 1:ncol(feat),
                               relevant_Variable_prop = 1:ncol(feat)) {
  hyperparameter_list <- list(
    "general" = list("predmode" = "propmean",
                     "nthread" = nthread),
    "l_first_0" = list(
      "relevant_Variable" = relevant_Variable_first,
      "ntree" = ntree,
      "replace" = TRUE,
      "sample.fraction" = .8,
      "mtry" = round(ncol(feat) / 2),
      "nodesizeSpl" = 1,
      "nodesizeAvg" = 5,
      "splitratio" = .5,
      "middleSplit" = FALSE
    ),
    "l_first_1" = list(
      "relevant_Variable" = relevant_Variable_first,
      "ntree" = ntree,
      "replace" = TRUE,
      "sample.fraction" = .8,
      "mtry" = round(ncol(feat) / 2),
      "nodesizeSpl" = 1,
      "nodesizeAvg" = 5,
      "splitratio" = .5,
      "middleSplit" = FALSE
    ),
    "l_second_0" = list(
      "relevant_Variable" = relevant_Variable_second,
      "ntree" = ntree,
      "replace" = TRUE,
      "sample.fraction" = .9,
      "mtry" = ncol(feat),
      "nodesizeSpl" = 5,
      "nodesizeAvg" = 3,
      "splitratio" = .5,
      "middleSplit" = FALSE
    ),
    "l_second_1" = list(
      "relevant_Variable" = relevant_Variable_second,
      "ntree" = ntree,
      "replace" = TRUE,
      "sample.fraction" = .9,
      "mtry" = ncol(feat),
      "nodesizeSpl" = 5,
      "nodesizeAvg" = 3,
      "splitratio" = .5,
      "middleSplit" = FALSE
    ),
    "l_prop" = list(
      "relevant_Variable" = relevant_Variable_prop,
      "ntree" = ntree,
      "replace" = TRUE,
      "sample.fraction" = 0.9,
      "mtry" = ncol(feat),
      "nodesizeSpl" = 5,
      "nodesizeAvg" = 3,
      "splitratio" = .5,
      "middleSplit" = FALSE
    )
  )
  return(hyperparameter_list)
}


get_setting_weak <- function(feat,
                             ntree,
                             nthread,
                             relevant_Variable_first = 1:ncol(feat),
                             relevant_Variable_second = 1:ncol(feat),
                             relevant_Variable_prop = 1:ncol(feat)) {
  hyperparameter_list <- list(
    "general" = list("predmode" = "propmean",
                     "nthread" = nthread),
    "l_first_0" = list(
      "relevant_Variable" = relevant_Variable_first,
      "ntree" = ntree,
      "replace" = FALSE,
      "sample.fraction" = .7,
      "mtry" = round(ncol(feat) / 2),
      "nodesizeSpl" = 1,
      "nodesizeAvg" = 1,
      "splitratio" = .4,
      "middleSplit" = FALSE
    ),
    "l_first_1" = list(
      "relevant_Variable" = relevant_Variable_first,
      "ntree" = ntree,
      "replace" = FALSE,
      "sample.fraction" = .7,
      "mtry" = round(ncol(feat) / 2),
      "nodesizeSpl" = 1,
      "nodesizeAvg" = 1,
      "splitratio" = .4,
      "middleSplit" = FALSE
    ),
    "l_second_0" = list(
      "relevant_Variable" = relevant_Variable_second,
      "ntree" = ntree,
      "replace" = FALSE,
      "sample.fraction" = 0.8,
      "mtry" = max(1, round(ncol(feat) / 5)),
      "nodesizeSpl" = 10,
      "nodesizeAvg" = 3,
      "splitratio" = .75,
      "middleSplit" = FALSE
    ),
    "l_second_1" = list(
      "relevant_Variable" = relevant_Variable_second,
      "ntree" = ntree,
      "replace" = FALSE,
      "sample.fraction" = 0.8,
      "mtry" = max(1, round(ncol(feat) / 5)),
      "nodesizeSpl" = 10,
      "nodesizeAvg" = 3,
      "splitratio" = .75,
      "middleSplit" = FALSE
    ),
    "l_prop" = list(
      "relevant_Variable" = relevant_Variable_prop,
      "ntree" = ntree,
      "replace" = FALSE,
      "sample.fraction" = 0.8,
      "mtry" = max(1, round(ncol(feat) / 5)),
      "nodesizeSpl" = 10,
      "nodesizeAvg" = 3,
      "splitratio" = .75,
      "middleSplit" = FALSE
    )
  )
  return(hyperparameter_list)
}



get_starting_settings <- function(feat,
                                  tr,
                                  ntree,
                                  nthread,
                                  relevant_Variable_first = 1:ncol(feat),
                                  relevant_Variable_second = 1:ncol(feat),
                                  relevant_Variable_prop = 1:ncol(feat)) {

  # this will contain the list of good parameter settings:
  hyperparameter_list_list <- list()
  for(i in 1:nrow(starting_values)){
    # read in setting i
    this_setting <- starting_values[i,]

    hyperparameter_list <- list(
      "general" = list("predmode" = "propmean",
                       "nthread" = nthread),
      "l_first_0" = list(
        "relevant_Variable" = relevant_Variable_first,
        "ntree" = ntree,
        "replace" = this_setting$replace_first,
        "sample.fraction" = .8,
        "mtry" = max(1, round(ncol(feat) * this_setting$mtry_first / 20)),
        "nodesizeSpl" = this_setting$nodesizeSpl_first,
        "nodesizeAvg" = this_setting$nodesizeAvg_first,
        "splitratio" = this_setting$splitratio_first,
        "middleSplit" = this_setting$middleSplit_first
      ),
      "l_first_1" = list(
        "relevant_Variable" = relevant_Variable_first,
        "ntree" = ntree,
        "replace" = this_setting$replace_first,
        "sample.fraction" = .8,
        "mtry" = max(1, round(ncol(feat) * this_setting$mtry_first / 20)),
        "nodesizeSpl" = this_setting$nodesizeSpl_first,
        "nodesizeAvg" = this_setting$nodesizeAvg_first,
        "splitratio" = this_setting$splitratio_first,
        "middleSplit" = this_setting$middleSplit_first
      ),
      "l_second_0" = list(
        "relevant_Variable" = relevant_Variable_second,
        "ntree" = ntree,
        "replace" = this_setting$replace_second,
        "sample.fraction" = .8,
        "mtry" = max(1, round(ncol(feat) * this_setting$mtry_second / 20)),
        "nodesizeSpl" = this_setting$nodesizeSpl_second,
        "nodesizeAvg" = this_setting$nodesizeAvg_second,
        "splitratio" = this_setting$splitratio_second,
        "middleSplit" = this_setting$middleSplit_second
      ),
      "l_second_1" = list(
        "relevant_Variable" = relevant_Variable_second,
        "ntree" = ntree,
        "replace" = this_setting$replace_second,
        "sample.fraction" = .8,
        "mtry" = max(1, round(ncol(feat) * this_setting$mtry_second / 20)),
        "nodesizeSpl" = this_setting$nodesizeSpl_second,
        "nodesizeAvg" = this_setting$nodesizeAvg_second,
        "splitratio" = this_setting$splitratio_second,
        "middleSplit" = this_setting$middleSplit_second
      ),
      "l_prop" = list(
        "relevant_Variable" = relevant_Variable_prop,
        "ntree" = ntree,
        "replace" = this_setting$replace_prop,
        "sample.fraction" = .8,
        "mtry" = max(1, round(ncol(feat) * this_setting$mtry_prop / 20)),
        "nodesizeSpl" = this_setting$nodesizeSpl_prop,
        "nodesizeAvg" = this_setting$nodesizeAvg_prop,
        "splitratio" = this_setting$splitratio_prop,
        "middleSplit" = this_setting$middleSplit_prop
      )
    )

    # make sure that the nodesizes are not too large:
    nodesizebounds <- get_upper_bounds_for_nodesize(
      starting_point = hyperparameter_list, tr)

    hyperparameter_list$l_first_0$nodesizeAvg <-
      hyperparameter_list$l_first_1$nodesizeAvg <- min(
        hyperparameter_list$l_first_0$nodesizeAvg,
        nodesizebounds["nodesizeAvg_first_upper"])

    hyperparameter_list$l_first_0$nodesizeSpl <-
      hyperparameter_list$l_first_1$nodesizeSpl <- min(
        hyperparameter_list$l_first_0$nodesizeSpl,
        nodesizebounds["nodesizeSpl_first_upper"])

    hyperparameter_list$l_second_0$nodesizeAvg <-
      hyperparameter_list$l_second_1$nodesizeAvg <- min(
        hyperparameter_list$l_second_0$nodesizeAvg,
        nodesizebounds["nodesizeAvg_second_upper"])

    hyperparameter_list$l_second_0$nodesizeSpl <-
      hyperparameter_list$l_second_1$nodesizeSpl <- min(
        hyperparameter_list$l_second_0$nodesizeSpl,
        nodesizebounds["nodesizeSpl_second_upper"])

    hyperparameter_list$l_prop$nodesizeAvg <-
      min(hyperparameter_list$l_prop$nodesizeAvg,
      floor(length(tr) * hyperparameter_list$l_prop$sample.fraction *
      ifelse(hyperparameter_list$l_prop$splitratio == 1, 1,
             (1 - hyperparameter_list$l_prop$splitratio))))

    hyperparameter_list$l_prop$nodesizeSpl <-
      min(hyperparameter_list$l_prop$nodesizeSpl,
          floor(length(tr) * hyperparameter_list$l_prop$sample.fraction *
                  ifelse(hyperparameter_list$l_prop$splitratio == 0, 1,
                  hyperparameter_list$l_prop$splitratio)))

    hyperparameter_list_list[[as.data.frame(starting_values)[i,1]]] <- hyperparameter_list
  }
  return(hyperparameter_list_list)
}


get_upper_bounds_for_nodesize <- function(starting_point, tr) {
  n_obs <- data.frame()
  for (this_learner in c("l_first_0",  "l_first_1", "l_second_0",
                         "l_second_1")) {
    # total amount of data possible for the tree:
    if (this_learner %in% c("l_first_0", "l_second_0")) {
      n <- sum(1 - tr)
    } else if (this_learner %in% c("l_first_1", "l_second_1")) {
      n <- sum(tr)
    }
    # fraction of observations per tree:
    sff <- starting_point[[this_learner]]$sample.fraction
    # ratio of those in the splitting set:
    sr <- starting_point[[this_learner]]$splitratio
    # amount of data for tree
    n_pertree <- floor(n * sff)
    # amount of data in the spl and avg set for m
    n_obs <- rbind(n_obs,
                   data.frame(
                     learner = this_learner,
                     n_avg = floor(n_pertree * ifelse(sr == 1, 1, (1 - sr))),
                     n_spl = floor(n_pertree * ifelse(sr == 0, 1, sr))
                   ))
  }
  return(
    c(
      nodesizeAvg_first_upper =  max(1, min(n_obs$n_avg[1:2])),
      nodesizeAvg_second_upper = max(1, min(n_obs$n_avg[3:4])),
      nodesizeSpl_first_upper =  max(1, min(n_obs$n_spl[1:2])),
      nodesizeSpl_second_upper = max(1, min(n_obs$n_spl[3:4]))
    )
  )
}

best_MSE_constant <- function(MSE_1, MSE_2){
  # Given two independent and unbiased estimators for theta with MSE_1 and
  # MSE_2, this function computes the MSE of the optimally combined estimator
  # which was combined in the following way:
  # alpha <- best_average_constant(MSE_1, MSE_2)
  # combined_estimator <- alpha theta_1 + (1 - alpha) theta_2
  alpha <- best_average_constant(MSE_1, MSE_2)
  return(alpha ^ 2 * MSE_1 + (1 - alpha) ^ 2 * MSE_2)
}

best_average_constant <- function(MSE_1, MSE_2){
  # Given two independent and unbiased estimators for theta with MSE_1 and
  # MSE_2, this function computes the optimal coefficient for combining those two
  # estimators in the following way:
  # alpha <- best_average_constant(MSE_1, MSE_2)
  # combined_estimator <- alpha theta_1 + (1 - alpha) theta_2
  alpha <- MSE_2 / (MSE_1 + MSE_2)
  return(alpha)
}

