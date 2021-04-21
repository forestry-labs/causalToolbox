#' @include XRF.R

#' @title tuneStageOne-forestry
#' @name tuneStageOne-forestry
#' @rdname tuneStageOne-forestry
#' @description Tune the estimator in the first stage given a fixed second stage
#' estimator
#' @param x A data frame of all training predictors.
#' @param y A vector of all training responses.
#' @param m_tau_init The initialize m_tau estimator.
#' @param sampsize The size of total samples to draw for the training data.
#' @param num_iter Maximum iterations/epochs per configuration. Default is 1024.
#' @param eta Downsampling rate. Default value is 2.
#' @param verbose if tuning process in verbose mode
#' @param seed A random seed
#' @param nthread Number of threads to train and predict thre forest. The
#' default number is 0 which represents using all cores.
#' @return A `forestry` object
#' @export tuneStageOne
tuneStageOne <- function(x,
                         y,
                         m_tau_init,
                         sampsize = as.integer(nrow(x) * 0.75),
                         num_iter = 1024,
                         eta = 2,
                         verbose = FALSE,
                         seed = 24750371,
                         nthread = 0) {

  # Creat a dummy tree just to reuse its data.
  dummy_tree <- Rforestry::forestry(x, y, ntree=1, nodesizeSpl=nrow(x), nodesizeAvg=nrow(x))

  # Number of unique executions of Successive Halving (minus one)
  s_max <- as.integer(log(num_iter) / log(eta))

  # Total number of iterations (without reuse) per execution of
  # successive halving (n,r)
  B <- (s_max + 1) * num_iter

  if (verbose) {
    print(
      paste(
        "Hyperband will run successive halving in",
        s_max,
        "times, with",
        B,
        "iterations per execution."
      )
    )
  }

  # Begin finite horizon hyperband outlerloop
  models <- vector("list", s_max + 1)
  models_OOB <- vector("list", s_max + 1)

  set.seed(seed)

  for (s in s_max:0) {
    if (verbose) {
      print(paste("Hyperband successive halving round", s_max + 1 - s))
    }

    # Initial number of configurations
    n <- as.integer(ceiling(B / num_iter / (s + 1) * eta ^ s))

    # Initial number of iterations to run configurations for
    r <- num_iter * eta ^ (-s)

    if (verbose) {
      print(paste(">>> Total number of configurations:", n))
      print(paste(
        ">>> Number of iterations per configuration:",
        as.integer(r)
      ))
    }

    # Begin finite horizon successive halving with (n,r)
    # Generate parameters:
    allConfigs <- data.frame(
      mtry = sample(1:ncol(x), n, replace = TRUE),
      min_node_size_spl = NA, #sample(1:min(30, nrow(x)), n, replace = TRUE),
      min_node_size_ave = NA, #sample(1:min(30, nrow(x)), n, replace = TRUE),
      splitratio = runif(n, min = 0.1, max = 1),
      replace = sample(c(TRUE, FALSE), n, replace = TRUE),
      middleSplit = sample(c(TRUE, FALSE), n, replace = TRUE)
    )

    min_node_size_spl_raw <- floor(allConfigs$splitratio * sampsize *
                                     rbeta(n, 1, 3))
    allConfigs$min_node_size_spl <- ifelse(min_node_size_spl_raw == 0, 1,
                                           min_node_size_spl_raw)
    min_node_size_ave <- floor((1 - allConfigs$splitratio) * sampsize *
                                 rbeta(n, 1, 3))
    allConfigs$min_node_size_ave <- ifelse(min_node_size_ave == 0, 1,
                                           min_node_size_ave)

    if (verbose) {
      print(paste(">>>", n, " configurations have been generated."))
    }

    val_models <- vector("list", nrow(allConfigs))
    r_old <- 1
    for (j in 1:nrow(allConfigs)) {
      tryCatch({
        val_models[[j]] <- Rforestry::forestry(
          x = x,
          y = y,
          ntree = r_old,
          mtry = allConfigs$mtry[j],
          nodesizeSpl = allConfigs$min_node_size_spl[j],
          nodesizeAvg = allConfigs$min_node_size_ave[j],
          splitratio = allConfigs$splitratio[j],
          replace = allConfigs$replace[j],
          sampsize = sampsize,
          nthread = nthread,
          middleSplit = allConfigs$middleSplit[j],
          reuseforestry=dummy_tree
        )

      }, error = function(err) {
        val_models[[j]] <- NULL
      })
    }

    if (s != 0) {
      for (i in 0:(s - 1)) {
        # Run each of the n_i configs for r_i iterations and keep best
        # n_i/eta
        n_i <- as.integer(n * eta ^ (-i))
        r_i <- as.integer(r * eta ^ i)
        r_new <- r_i - r_old

        # if (verbose) {
        #   print(paste("Iterations", i))
        #   print(paste("Total number of configurations:", n_i))
        #   print(paste("Number of iterations per configuration:", r_i))
        # }

        val_losses <- vector("list", nrow(allConfigs))

        # Iterate to evaluate each parameter combination and cut the
        # parameter pools in half every iteration based on its score
        for (j in 1:nrow(allConfigs)) {
          if (r_new > 0 && !is.null(val_models[[j]])) {
            val_models[[j]] <- Rforestry::addTrees(val_models[[j]], r_new)
          }
          if (!is.null(val_models[[j]])) {
            # If the model is available, get its OOB error
            val_losses[[j]] <- Rforestry::getOOB(val_models[[j]], noWarning = TRUE)
            # Calculate residuals
            res <- predict(val_models[[j]], x) - y
            # Train an forestry for tau

            m_tau <-
              Rforestry::forestry(
                x = x,
                y = res,
                ntree = max(r_i, 1),
                replace = m_tau_init@replace,
                sampsize = m_tau_init@sampsize,
                mtry = m_tau_init@mtry,
                nodesizeSpl = m_tau_init@nodesizeSpl,
                nodesizeAvg = m_tau_init@nodesizeAvg,
                splitratio = m_tau_init@splitratio,
                seed = seed,
                verbose = FALSE,
                nthread = nthread,
                splitrule = "variance",
                middleSplit = m_tau_init@middleSplit
              )
            # If the tau model is valid, adding its OOB to the existing OOB
            if (!is.null(m_tau)) {
              tau_oob <- Rforestry::getOOB(m_tau, noWarning = TRUE)
              val_losses[[j]] <- val_losses[[j]] + tau_oob
            } else {
              val_losses[[j]] <- NA
            }

            if (is.na(val_losses[[j]])) {
              val_losses[[j]] <- Inf
            }
          } else {
            val_losses[[j]] <- Inf
          }
        }

        r_old <- r_i

        val_losses_idx <-
          sort(unlist(val_losses), index.return = TRUE)
        val_top_idx <- val_losses_idx$ix[0:as.integer(n_i / eta)]
        allConfigs <- allConfigs[val_top_idx,]
        val_models <- val_models[val_top_idx]
        gc()
        rownames(allConfigs) <- 1:nrow(allConfigs)

        # if (verbose) {
        #   print(paste(length(val_losses_idx$ix) - nrow(allConfigs),
        #               "configurations have been eliminated."))
        # }

      }

    }
    # End finite horizon successive halving with (n,r)
    if (!is.null(val_models[[1]])) {
      best_OOB <- Rforestry::getOOB(val_models[[1]], noWarning = TRUE)
      res <- predict(val_models[[1]], x) - y
      m_tau <-
        Rforestry::forestry(
          x = x,
          y = res,
          ntree = m_tau_init@ntree,
          replace = m_tau_init@replace,
          sampsize = m_tau_init@sampsize,
          mtry = m_tau_init@mtry,
          nodesizeSpl = m_tau_init@nodesizeSpl,
          nodesizeAvg = m_tau_init@nodesizeAvg,
          splitratio = m_tau_init@splitratio,
          seed = seed,
          verbose = FALSE,
          nthread = nthread,
          splitrule = "variance",
          middleSplit = m_tau_init@middleSplit
        )

      # If the tau model is valid, adding its OOB to the existing OOB
      if (!is.null(m_tau)) {
        best_OOB <- best_OOB + Rforestry::getOOB(m_tau, noWarning = TRUE)
      } else {
        best_OOB <- NA
      }
      if (is.na(best_OOB)) {
        stop()
        best_OOB <- Inf
      }
    } else {
      stop()
      best_OOB <- Inf
    }
    if (verbose) {
      print(paste(">>> Successive halving ends and the best model is saved."))
      print(paste(">>> OOB:", best_OOB))
    }

    if (!is.null(val_models[[1]]))
      models[[s + 1]] <- val_models[[1]]
    models_OOB[[s + 1]] <- best_OOB

  }

  # End finite horizon hyperband outlerloop and sort by performance
  model_losses_idx <- sort(unlist(models_OOB), index.return = TRUE)

  if (verbose) {
    print(
      paste(
        "Best model is selected from best-performed model in",
        s_max,
        "successive halving, with OOB",
        models_OOB[model_losses_idx$ix[1]]
      )
    )
  }

  return(models[[model_losses_idx$ix[1]]])

}

#' @title Autotuning for X-Learner with honest RF for both stages
#' @name autoJointforestry
#' @rdname autoJointforestry
#' @description [TBD]
#' @param feat A data frame of all the features.
#' @param tr A numeric vector contain 0 for control and 1 for treated variables.
#' @param yobs A numeric vector containing the observed outcomes.
#' @param sample.fraction ..
#' @param num_iter number of iterations.
#' @param eta ..
#' @param verbose ..
#' @param seed ..
#' @param nthread ..
#' @param tune_iter ..
#' @seealso \code{\link{X_RF_autotune_simple}}, \code{\link{X_RF_autotune_gpp}},
#' @export autoJointforestry
autoJointforestry <-
  function(feat,
           tr,
           yobs,
           sample.fraction = 0.75,
           num_iter = 3 ^ 8,
           eta = 3,
           verbose = TRUE,
           seed = 24750371,
           nthread = 0,
           tune_iter = 10) {
    feat <- as.data.frame(feat)

    hyperparameter_list <- list()
    base_learners <- list()

    yobs_0 <- yobs[tr == 0]
    yobs_1 <- yobs[tr == 1]

    X_0 <- feat[tr == 0,]
    X_1 <- feat[tr == 1,]

    # First, find the best configurations for both estimators
    m_0 <-
      Rforestry::autoforestry(
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
      Rforestry::autoforestry(
        x = X_1,
        y = yobs_1,
        sampsize = floor(nrow(X_1) * sample.fraction),
        num_iter = num_iter,
        eta = eta,
        verbose = verbose,
        seed = seed,
        nthread = nthread
      )

    r_0 <- predict(m_1, X_0) - yobs_0
    r_1 <- yobs_1 - predict(m_0, X_1)

    m_tau_0 <-
      Rforestry::autoforestry(
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
      Rforestry::autoforestry(
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
      print("Initialize configurations using hyperband.")
    }


    for (i in 1:tune_iter) {
      print(paste("Start joint tuning, iteration", i))

      m_0 <- tuneStageOne(
        x = X_0,
        y = yobs_0,
        m_tau_init = m_tau_0,
        sampsize = floor(nrow(X_0) * sample.fraction),
        num_iter = num_iter,
        eta = eta,
        verbose = verbose,
        seed = seed + i,
        nthread = nthread
      )

      m_1 <- tuneStageOne(
        x = X_1,
        y = yobs_1,
        m_tau_init = m_tau_1,
        sampsize = floor(nrow(X_1) * sample.fraction),
        num_iter = num_iter,
        eta = eta,
        verbose = verbose,
        seed = seed + i,
        nthread = nthread
      )

      r_0 <- predict(m_1, X_0) - yobs_0
      r_1 <- yobs_1 - predict(m_0, X_1)

      m_tau_0 <-
        Rforestry::autoforestry(
          x = X_0,
          y = r_0,
          sampsize = floor(nrow(X_0) * sample.fraction),
          num_iter = num_iter,
          eta = eta,
          verbose = verbose,
          seed = seed + i,
          nthread = nthread
        )

      m_tau_1 <-
        Rforestry::autoforestry(
          x = X_1,
          y = r_1,
          sampsize = floor(nrow(X_1) * sample.fraction),
          num_iter = num_iter,
          eta = eta,
          verbose = verbose,
          seed = seed + i,
          nthread = nthread
        )

    }


    m_prop <-
      Rforestry::forestry(x = feat,
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
