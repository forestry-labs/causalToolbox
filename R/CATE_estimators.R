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
#'@param nthread Number of threads to be used in parallel.
#'@param verbose TRUE for detailed output, FALSE for no output.
#'@param bootstrapVersion Default is normalApprox, which will use the 
#'bootstrap normal approximation to get CI. Smoothed will use CI around the
#'smoothed bootstrap as introduced by Efron 2014.
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
    }

    #### Run the bootstrap CI estimation #####################################

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
            S[names(counts), b] <- counts
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
      # normal Approximated Bootstrap ------------------------------------------
      # get the predictions from the original method
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
    } else {
      stop("bootstrapVersion must be specified.")
    }
  }
)





