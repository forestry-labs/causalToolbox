#' @import Rforestry
NULL

# Cate estimators --------------------------------------------------------------
setClass("CATEestimator",
         slots = list(
           feature_train = "data.frame",
           tr_train = "numeric",
           yobs_train = "numeric",
           creator = "function"
         ))
setOldClass("forestry::honestRF")

# MetaLearner -----------------------------------------------------------------
setClass(
  "MetaLearner",
  contains = "CATEestimator"
)
#'Method EstimateCate
#'@name EstimateCate
#'@rdname EstimateCate
#'@description Returns the estimated CATE. 
#'@param theObject A `MetaLearner` object.
#'@param feature_new A feature data frame.
#'@param ... Additional parameters that are specific for some MetaLearners
#'@export
setGeneric(
  name = "EstimateCate",
  def = function(theObject, feature_new, ...) {
    standardGeneric("EstimateCate")
  }
)
#'Method CateCI
#'@name CateCI
#'@rdname CateCI
#'@param theObject A `MetaLearner` object.
#'@param feature_new A feature data frame.
#'@param method Different versions of the bootstrap.
#'@param B Number of bootstrap samples.
#'@param B_Second Number of bootstrap samples to take in the second layer of the
#' double bootstrap (the calibration samples). By default this is equal to B, 
#' however in practice we suggest using a slightly smaller value as the runtime
#' is constrained by O(B * B_Second).
#'@param nthread Number of threads to be used in parallel.
#'@param verbose TRUE for detailed output, FALSE for no output.
#'@param bootstrapVersion Default is normalApprox, which will use the 
#' bootstrap normal approximation to get CI. Smoothed will use CI around the
#' smoothed bootstrap as introduced by Efron 2014. The third option is to use 
#' the doubleBootstrap option, which uses a double level bootstrap to calibrate 
#' the quantiles used in the bootstrap estimation of the intervals. For reference
#' see https://arxiv.org/pdf/1511.00273.pdf, although this is an older algorithm
#' which was introduced much earlier.
#'@export
#' @examples
#' \dontrun{
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
#' CateCI(xl_rf, feature_test, B = 500)
#' }
setGeneric(
  name = "CateCI",
  def = function(theObject,
                 feature_new,
                 method = "maintain_group_ratios",
                 bootstrapVersion = "normalApprox",
                 B = 2000,
                 B_Second = B,
                 nthread = 0,
                 verbose = TRUE) {
    standardGeneric("CateCI")
  }
)

# Confidence Interval Estimation -----------------------------------------------
# Estimating Confidence intervals

#' CateCI-MetaLearner
#' @rdname CateCI-MetaLearner
#' @description Returns the estimated confidence intervals for the CATE.
#' @return A data frame of estimated CATE confidence intervals.
#' @inherit CateCI
#' @export
setMethod(
  f = "CateCI",
  signature = "CATEestimator",
  definition = function(theObject,
                        feature_new,
                        method,
                        bootstrapVersion,
                        B,
                        B_Second,
                        nthread,
                        verbose) {
    ## shortcuts:
    feat <- theObject@feature_train
    tr <- theObject@tr_train
    yobs <- theObject@yobs_train
    creator <- theObject@creator
    ntrain <- length(tr)
    if ((bootstrapVersion == "smoothed") & 
       (as.double(nrow(feat)) * as.double(nrow(feature_new)) > 5e8)) {
      stop(paste("We would have to create a", nrow(feat), 
                 "by", nrow(feature_new), "matrix. This is too big to run in",
                 "a reasonable amount of time. Either decrease feature_new",
                 "or use `bootstrapVersion <- normalApprox` option.",
                 "The matrix should have less than 5e8 values."))
    }
    if ((bootstrapVersion == "smoothed") & B < 2000) {
      warning(
        paste(
          "We have found that when using smoothed intervals,",
          "B should be chosen to be bigger than 2000."
        )
      )
    }
    
    if (bootstrapVersion == "doubleBootstrap" & verbose & B >= 1000) {
      warning(
        paste("Warning, the doubleBootstrap uses B^2 = ", B^2, " bootstraps",
              "On larger data sets this may take some time")
      )
    }
    
    if (bootstrapVersion == "conformal") {
      # standard conformal intervals ------------------------------------------
      
      # The local conformal intervals make use of the weighting function of RF,
      # so for now, this can only be used with S_RF, trained with OOBhonest = TRUE
      if ((class(theObject)[1] != "S_RF") | (!theObject@hyperparameter_list$mu.forestry$OOBhonest)) {
        stop("For local conformal intervals, we must use S_RF")
      }
      if (theObject@forest@OOBhonest) {
        pred_insample <- EstimateCate(theObject = theObject,
                                      feature_new = theObject@feature_train,
                                      aggregation = "oob")
      } else {
        pred_insample <- EstimateCate(theObject = theObject,
                                      feature_new = theObject@feature_train,
                                      aggregation = "average")
      }
      
      
      quantiles <- quantile(pred_insample, probs = c(.025, .975))
      
      pred <- EstimateCate(theObject, feature_new = feature_new)
      
      return(data.frame(
        pred = pred,
        X5. =  pred + unname(quantiles)[1],
        X95. = pred + unname(quantiles)[2]
        # X5. =  pred - (CI_b$X95. - CI_b$X5.) / 2,
        # X95. = pred + (CI_b$X95. - CI_b$X5.) / 2
        # X5. =  2 * pred - CI_b$X95.,
        # X95. = 2 * pred - CI_b$X5.
      ))
      
    }
    
    if (method == "maintain_group_ratios") {
      createbootstrappedData <- function() {

        smpl_0 <- sample((1:ntrain)[tr == 0],
                       replace = TRUE,
                       size = sum(1 - tr))
        smpl_1 <- sample((1:ntrain)[tr == 1],
                         replace = TRUE,
                         size = sum(tr))
        smpl <- sample(c(smpl_0, smpl_1))

        return(list(
          feat_b = feat[smpl, ],
          tr_b = tr[smpl],
          yobs_b = yobs[smpl], 
          smpl = smpl
        ))
      }
      
      createSecondBootstrappedData <- function(
        original_indices
      ) {
        
        smpl_0 <- sample(original_indices[tr == 0],
                         replace = TRUE,
                         size = sum(1 - tr))
        smpl_1 <- sample(original_indices[tr == 1],
                         replace = TRUE,
                         size = sum(tr))
        smpl <- sample(c(smpl_0, smpl_1))
        
        return(list(
          feat_b = feat[smpl, ],
          tr_b = tr[smpl],
          yobs_b = yobs[smpl], 
          smpl = smpl
        ))
      }
    }

    # Run the bootstrap CI estimation #####################################

    # pred_B will contain for each simulation the prediction of each of the B
    # simulaions:
    pred_B <-
      as.data.frame(matrix(NA, nrow = nrow(feature_new), ncol = B))
  
    # S is needed for Efron's smooth bootstrap each column corresponse to one 
    # bootstrap sample and each row corresponse to one of the smple indexes
    S <- as.data.frame(matrix(0, nrow = length(yobs), ncol = B))
    row.names(S) <- 1:length(yobs)
    colnames(S) <- 1:B
    
    
    known_warnings <- c()
    # this is needed such that bootstrapped warnings are only printed once
    for (b in 1:B) { # b= 1
      if (verbose)
        print(b)
      went_wrong <- 0
      # if that is 100 we really cannot fit it and bootstrap
      # seems to be infeasible.

      while (is.na(pred_B[1, b])) {
        if (went_wrong == 100)
          stop("one of the groups might be too small to
               do valid inference.")
        S[, b] <- rep(0, nrow(S))
        
        pred_B[, b] <-
          tryCatch({
            bs <- createbootstrappedData()
            
            counts <- table(bs$smpl)
            if (bootstrapVersion == "doubleBootstrap") {
              S[, b] <- bs[["smpl"]]
            } else {
              S[names(counts), b] <- counts
            }

            withCallingHandlers(
              # this is needed such that bootstrapped warnings are only
              # printed once
              EstimateCate(
                creator(
                  feat = bs$feat_b,
                  tr = bs$tr_b,
                  yobs = bs$yobs_b
                ),
                feature_new = feature_new
              ),
              warning = function(w) {
                if (w$message %in% known_warnings) {
                  # message was already printed and can be ignored
                  invokeRestart("muffleWarning")
                } else{
                  # message is added to the known_warning list:
                  known_warnings <<- c(known_warnings, w$message)
                }
              }
            )
          },
          error = function(e) {
            return(NA)
          })
        went_wrong <- went_wrong + 1
      }
    }

    if (bootstrapVersion == "normalApprox") {
      
      # Normal Approximated Bootstrap -----------------------------------------
      
      pred <- EstimateCate(theObject, feature_new = feature_new)
      # the the 5% and 95% CI from the bootstrapped procedure
      CI_b <- data.frame(
        X5. =  apply(pred_B, 1, function(x)
          quantile(x, c(.025))),
        X95. = apply(pred_B, 1, function(x)
          quantile(x, c(.975))),
        sd = apply(pred_B, 1, function(x) sd(x))
      )
      
      return(data.frame(
        pred = pred,
        X5. =  pred - 1.96 * CI_b$sd,
        X95. = pred + 1.96 * CI_b$sd
        # X5. =  pred - (CI_b$X95. - CI_b$X5.) / 2,
        # X95. = pred + (CI_b$X95. - CI_b$X5.) / 2
        # X5. =  2 * pred - CI_b$X95.,
        # X95. = 2 * pred - CI_b$X5.
      ))
    } else if (bootstrapVersion == "smoothed") {
      # Smoothed Bootstrap -----------------------------------------------------
      smoothed_mean <- apply(pred_B, 1, mean)
      

      pred_term <- as.matrix(
        pred_B - 
          matrix( smoothed_mean,
                  nrow = length(smoothed_mean),
                  ncol = B,
                  byrow = FALSE
          ))
      
      S_term <- as.matrix(
        S -
          matrix(apply(S, 1, mean),
                 nrow = nrow(S),
                 ncol = B,
                 byrow = FALSE))
      var_sol <- apply(((pred_term %*% t(S_term)) / (B - 1) )^2, 1, sum)
      
      return(data.frame(
        pred = smoothed_mean,
        X5. =  smoothed_mean - 1.96 * sqrt(var_sol),
        X95. = smoothed_mean + 1.96 * sqrt(var_sol)))
      
    } else if (bootstrapVersion == "doubleBootstrap") {
      # Double Bootstrap -------------------------------------------------------
      library(dplyr)
      B_1 <- B
      B_2 <- B_Second
      # For doubleBootstrap, we now do B bootstrap resamples of each original
      # bootstrap sample. We then estimate the predictions using each resampled
      # double bootstraps. This gives us B^2 total prediction vectors.
      # Then we find the largest quantile Lambda such that in 95% of the B 
      # doubleBootstrap sets, the standard predictions are contained in the 
      # 1-lambda and lambda percentiles of the second layer bootstrap predicions.
      # One idea is maybe we use B_2 ~ O(B_1^q) with .5 < q < 1, this is around 
      # the asymptotic
      # number for Wager's Infinitesimal Jacknife, and follows Hall's heuristic
      # to use a smaller number of double bootstraps than first level bootstraps
      
      
      jobs <- as.data.frame(expand.grid(1:B_2, 1:B_1))
      colnames(jobs) <- c("B_2", "B_1")
      jobs <- cbind(jobs, data.frame(matrix(0,nrow = B_1*B_2, ncol = nrow(feature_new))))
      
      if (verbose) {
        print("Running second layer of bootstraps")
      }
      
      # For now do second layer of bootstraps in serial
      for (i in 1:B_1) {
        # Get current outside bootstrap indices
        current_bootstrap_idx <- S[,i]
        
        for (j in 1:B_2) {
          if (verbose) {
            print(paste0("Running ",((i-1)*B_2+j)," of ", B_1*B_2))
          }

          jobs[which(jobs$B_2 == j & jobs$B_1 == i), c(-1,-2)] <-
            tryCatch({
              second_bootstrap <- createSecondBootstrappedData(current_bootstrap_idx)
      
              withCallingHandlers(
                # this is needed such that bootstrapped warnings are only
                # printed once
                EstimateCate(
                  creator(
                    feat = second_bootstrap$feat_b,
                    tr = second_bootstrap$tr_b,
                    yobs = second_bootstrap$yobs_b
                  ),
                  feature_new = feature_new
                ),
                warning = function(w) {
                  if (w$message %in% known_warnings) {
                    # message was already printed and can be ignored
                    invokeRestart("muffleWarning")
                  } else{
                    # message is added to the known_warning list:
                    known_warnings <<- c(known_warnings, w$message)
                  }
                }
              )
            },
            error = function(e) {
              return(NA)
            })
        }
      }
      
      # Now we have the double bootstrap estimates, we need to calibrate the
      # CI based on 
      # We do this by starting with a quantile, .95, and pulling this for each
      # inner bootstrap. If the coverage is too low, we decrease the quantile, 
      # if it is too high, we increase the quantile until the calibration gets 
      # close to .95 across the outer layer of bootstraps
      if (verbose) {
        print("Calibrating quantiles")
      }
      
      lambda <- .95
      converged <- FALSE
      max_reps <- 30
      max_close_reps <- 3
      reps <- 0
      close_reps <- 0
      
      # Get CATE estimates
      pred <- EstimateCate(theObject, feature_new = feature_new)
      
      while (!converged) {
        # Get upper and lower bounds of each second layer bootstrap at levels lambda/2, (1-lambda)/2
        jobs %>%
          group_by(as.factor(B_1)) %>%
          summarise(across(starts_with("X"), 
                    list(upper = 
                           function(x){return(quantile(x, probs = c(1-(1-lambda)/2)))},
                         lower = 
                           function(x){return(quantile(x, probs = c((1-lambda)/2)))}
                         ))) -> quants
        
        quants %>% 
          dplyr::select(contains("lower")) %>% 
          t() -> lower_bounds
        
        quants %>% 
          dplyr::select(contains("upper")) %>% 
          t() -> upper_bounds
        
        # Get coverage of the estimated predictions across first level of bootstrap
        # using current quantiles quantile(probs = c((1-lambda)/2, 1-(1-lambda)/2)))
        coverage_lower <- apply(lower_bounds,
                          MARGIN = 2,
                          FUN = function(x){return(x < pred)})
        
        coverage_upper <- apply(upper_bounds,
                          MARGIN = 2,
                          FUN = function(x){return(x > pred)})
        
        (coverage_lower & coverage_upper) %>% 
          apply(MARGIN = 1, FUN = function(x){return(length(which(x))/length(x))}) %>% 
          unname() %>% 
          median() -> med_point_coverage
            
        # We have two options here, we can either adjust the quantiles separately
        # for each point, or we can adjust the quantiles for the set- using the 
        # median coverage. We do the second for now
        
        # If coverage > .95 increase lambda, if coverage < .95 reduce lambda, 
        # if coverage ~= .95, return the right quantiles from the first layer of bootstraps
        reps <- reps + 1
        if (reps > max_reps ){
          break
        }
        if (med_point_coverage < .94) {
          lambda <- min(lambda + .01, .995)
        } else if (med_point_coverage > .96) {
          lambda <- min(lambda - .01, .995)
        } else {
          close_reps <- close_reps + 1
          if (close_reps > max_close_reps){
            break
          }
          if (med_point_coverage < .95) {
            lambda <- min(lambda + .05, .995)
          } else {
            lambda <- min(lambda - .05, .995)
          }
        }
      }
      
      # Now we have calibrated the lambda, we want to return the quantiles
      # of the first layer of bootstrap estimates according to 1-(1-lambda)/2 and (1-lambda)/2
      
      return(list("CI" = data.frame(
        pred = pred,
        # Get the lower quantile based on lambda
        X5. = apply(pred_B, 
                    MARGIN = 1, 
                    FUN = function(x){return(quantile(x, probs = c((1-lambda)/2)))}),
        # Get the upper quantile based on lambda
        X95. = apply(pred_B, 
                     MARGIN = 1, 
                     FUN = function(x){return(quantile(x, probs = c(1-(1-lambda)/2)))})
      ), "lambda" = lambda))
      
    } else {
      stop("bootstrapVersion must be specified.")
    }
  }
)





