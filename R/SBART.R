#'@import BART
#'@include CATE_estimators.R


## the standard Xlearner object with random forest
setClass(
  "S_BART",
  contains = "MetaLearner",
  slots = list(
    feature_train = "data.frame",
    tr_train = "numeric",
    yobs_train = "numeric",
    ndpost = "numeric",
    ntree = "numeric",
    nthread = "numeric",
    mu.BART = "list", 
    creator = "function"
  )
)

#' @rdname Slearners
#' @description \code{S_BART} is an implementation of the S-Learner with Bayesian
#'   Additive Regression Trees (Chipman et al. 2010).
#' @param mu.BART hyperparameters of the BART function. Use
#'   \code{?BART::mc.wbart} for a detailed explanation of their effects.
#' @import methods
#' @inherit S_RF details
#' @inherit T_BART
#' @family metalearners
#' @export
S_BART <-
  function(feat,
           tr,
           yobs,
           ndpost = 1200,
           ntree = 200,
           nthread = 1,
           verbose = FALSE, 
           mu.BART = list(
             sparse = FALSE,
             theta = 0,
             omega = 1,
             a = 0.5,
             b = 1,
             augment = FALSE,
             rho = NULL,
             usequants = FALSE,
             cont = FALSE,
             sigest = NA,
             sigdf = 3,
             sigquant = 0.90,
             k = 2.0,
             power = 2.0,
             base = .95,
             sigmaf = NA,
             lambda = NA,
             numcut = 100L,
             nskip = 100L
           )) {
    
    feat <- as.data.frame(feat)

    new(
      "S_BART",
      feature_train = feat,
      tr_train = tr,
      yobs_train = yobs,
      ndpost = ndpost,
      ntree = ntree,
      nthread = nthread,
      mu.BART = mu.BART, 
      creator = function(feat, tr, yobs) {
        S_BART(feat,
               tr,
               yobs,
               ndpost = ndpost,
               ntree = ntree,
               nthread = nthread,
               verbose = verbose, 
               mu.BART = mu.BART)
      }
    )
  }

#' EstimateCate-S_BART
#' @rdname EstimateCate-S_BART
#' @inherit EstimateCate
#' @exportMethod EstimateCate
setMethod(
  f = "EstimateCate",
  signature = "S_BART",
  definition = function(theObject, feature_new, verbose = FALSE)
  {
    ndpost <- theObject@ndpost
    n_feature_new <- nrow(feature_new)

    pred_matrix <- get_pred_mat(
      theObject = theObject,
      feature_new = feature_new,
      ndpost = ndpost, 
      nthread = theObject@nthread)

    pred_0 <- apply(pred_matrix, 2, mean)[1:n_feature_new]
    pred_1 <- 
      apply(pred_matrix, 2, mean)[(n_feature_new + 1):(2 * n_feature_new)]

    return(pred_1 - pred_0)
  }
)


#' CateCI-S_BART
#' @rdname CateCI-S_BART
#' @aliases CateCI,S_BART-method
#' @exportMethod CateCI
setMethod(
  f = "CateCI",
  signature = "S_BART",
  definition = function(theObject, feature_new, verbose = FALSE)
  {
    ndpost <- theObject@ndpost
    n_feature_new <- nrow(feature_new)
    
    pred_matrix <- get_pred_mat(
      theObject = theObject,
      feature_new = feature_new,
      ndpost = ndpost,
      nthread = theObject@nthread
    )

    ite_matrix <-
      pred_matrix[, (n_feature_new + 1):(2 * n_feature_new)] -
      pred_matrix[, 1:n_feature_new]


    pred <- apply(ite_matrix, 2, mean)

    CI <-
      t(apply(ite_matrix, 2, function(x)
        quantile(x, probs = c(.05, 0.95))))

    to_return <- as.data.frame(cbind(pred,  CI))
    row.names(to_return) <- 1:nrow(to_return)
    colnames(to_return) <- c('pred','X5.','X95.')
    return(to_return)
  }
)


# Helper Functions -------------------------------------------------------------

get_pred_mat <- function(theObject, feature_new, ndpost, nthread) {
  feature_new <- as.data.frame(feature_new)
  n_feature_new <- nrow(feature_new)
  pred_matrix <-
    BART::mc.wbart(
      x.train = cbind(theObject@feature_train, tr = theObject@tr_train),
      y.train = theObject@yobs_train,
      x.test = cbind(rbind(feature_new, feature_new),
                     tr = c(
                       rep(0, n_feature_new), rep(1, n_feature_new)
                     )),
      ndpost = ndpost,
      ntree = theObject@ntree,
      mc.cores = nthread,
      sparse =    theObject@mu.BART$sparse,
      theta =     theObject@mu.BART$theta,
      omega =     theObject@mu.BART$omega,
      a =         theObject@mu.BART$a,
      b =         theObject@mu.BART$b,
      augment =   theObject@mu.BART$augment,
      rho =       theObject@mu.BART$rho,
      usequants = theObject@mu.BART$usequants,
      cont =      theObject@mu.BART$cont,
      sigest =    theObject@mu.BART$sigest,
      sigdf =     theObject@mu.BART$sigdf,
      sigquant =  theObject@mu.BART$sigquant,
      k =         theObject@mu.BART$k,
      power =     theObject@mu.BART$power,
      base =      theObject@mu.BART$base,
      sigmaf =    theObject@mu.BART$sigmaf,
      lambda =    theObject@mu.BART$lambda,
      numcut =    theObject@mu.BART$numcut,
      nskip =     theObject@mu.BART$nskip
    )$yhat.test

  return(pred_matrix)
}


typeof(as.integer(round(.2)))
typeof(2L)
