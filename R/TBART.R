#' @include CATE_estimators.R
#' @include SBART.R
#' @import BART



## the standard Xlearner object with random forest
setClass(
  "T_BART",
  contains = "MetaLearner",
  slots = list(
    feature_train = "data.frame",
    tr_train = "numeric",
    yobs_train = "numeric",
    ndpost = "numeric",
    ntree = "numeric",
    nthread = "numeric",
    mu0.BART = "list", 
    mu1.BART = "list",
    creator = "function"
  )
)


# T_BART generator -------------------------------------------------------------
#' @rdname Tlearners
#' @description \code{T_BART} is an implementation of the T-learner with Bayesian
#'   Additive Regression Trees (Chipman et al. 2010) for both response
#'   functions.
#' @param ndpost Number of posterior draws.
#' @param ntree Number of trees.
#' @param mu0.BART,mu1.BART Hyperparameters of the BART functions for the
#'   control and treated group. (Use \code{?BART::mc.wbart} for a detailed 
#'   explanation of their effects.)
#' @inherit T_RF details
#' @inherit X-Learner
#' @family NULL
#' @export
T_BART <-
  function(feat,
           tr,
           yobs,
           ndpost = 1200,
           ntree = 200,
           nthread = 1,
           verbose = FALSE, 
           mu0.BART = list(
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
           ), 
           mu1.BART = list(
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
      "T_BART",
      feature_train = feat,
      tr_train = tr,
      yobs_train = yobs,
      ndpost = ndpost,
      ntree = ntree,
      nthread = nthread,
      mu0.BART = mu0.BART, 
      mu1.BART = mu1.BART,
      creator = function(feat, tr, yobs) {
        T_BART(feat,
               tr,
               yobs,
               ndpost = ndpost,
               ntree = ntree,
               nthread = nthread,
               mu0.BART = mu0.BART, 
               mu1.BART = mu1.BART)
      }
    )
  }

#' EstimateCate-T_BART-method
#' @rdname EstimateCate-T_BART-method
#' @inherit EstimateCate
#' @exportMethod EstimateCate
#' @param verbose TRUE for detailed output FALSE for no output
#' @param return_CI If TRUE, return predictions and their confidence intervals;
#' @import parallel
setMethod(
  f = "EstimateCate",
  signature = "T_BART",
  definition = function(theObject,
                        feature_new,
                        verbose = FALSE,
                        return_CI = FALSE)
  {
    # theObject = xb;  verbose = TRUE; ndpost = 100; return_CI = TRUE;
    # feature_new = feat[1:5,]
    ndpost <- theObject@ndpost
    yobs <- theObject@yobs_train
    if ((is.null(yobs)) | (is.data.frame(yobs) && nrow(yobs) == 0)) {
      stop('Observations cannot be empty or null.')
    }
    feat <- theObject@feature_train
    if ((is.null(feat)) | (is.data.frame(feat) && nrow(feat) == 0)) {
      stop('Features cannot be empty or null.')
    }
    tr <- theObject@tr_train
    if ((is.null(tr)) | (is.data.frame(tr) && nrow(tr) == 0)) {
      stop('Treatment assignments cannot be empty or null.')
    }
    nthread <- theObject@nthread
    if (nthread > parallel::detectCores()) {
      warning('nthread is chosen bigger than the number of cores. It is changed
              to be equal to the number of cores.')
      nthread <- parallel::detectCores()
    }

    yobs_0 <- yobs[tr == 0]
    if (is.data.frame(yobs_0) && nrow(yobs_0) == 0) {
      stop('There is no observation in the control group (labelled 0).')
    }
    X_0 <- feat[tr == 0, ]
    if (is.data.frame(X_0) && nrow(X_0) == 0) {
      stop('There is no feature in the control group (labelled 0).')
    }
    yobs_1 <- yobs[tr == 1]
    if (is.data.frame(yobs_1) && nrow(yobs_1) == 0) {
      stop('There is no observation in the treatment group (labelled 1).')
    }
    X_1 <- feat[tr == 1, ]
    if (is.data.frame(X_1) && nrow(X_1) == 0) {
      stop('There is no feature in the treatment group (labelled 1).')
    }

    pred_matrix_f_0 <- BART::mc.wbart(
      x.train = X_0,
      y.train = yobs_0,
      x.test =  feature_new,
      ndpost = ndpost,
      ntree = theObject@ntree,
      mc.cores = nthread, 
      sparse =    theObject@mu0.BART$sparse,
      theta =     theObject@mu0.BART$theta,
      omega =     theObject@mu0.BART$omega,
      a =         theObject@mu0.BART$a,
      b =         theObject@mu0.BART$b,
      augment =   theObject@mu0.BART$augment,
      rho =       theObject@mu0.BART$rho,
      usequants = theObject@mu0.BART$usequants,
      cont =      theObject@mu0.BART$cont,
      sigest =    theObject@mu0.BART$sigest,
      sigdf =     theObject@mu0.BART$sigdf,
      sigquant =  theObject@mu0.BART$sigquant,
      k =         theObject@mu0.BART$k,
      power =     theObject@mu0.BART$power,
      base =      theObject@mu0.BART$base,
      sigmaf =    theObject@mu0.BART$sigmaf,
      lambda =    theObject@mu0.BART$lambda,
      numcut =    theObject@mu0.BART$numcut,
      nskip =     theObject@mu0.BART$nskip
    )$yhat.test
    
    mu_hat_0 <- apply(pred_matrix_f_0, 2, mean)
    
    pred_matrix_f_1 <- BART::mc.wbart(
      x.train = X_1,
      y.train = yobs_1,
      x.test =  feature_new,
      ndpost = ndpost,
      ntree = theObject@ntree,
      mc.cores = nthread, 
      sparse =    theObject@mu1.BART$sparse,
      theta =     theObject@mu1.BART$theta,
      omega =     theObject@mu1.BART$omega,
      a =         theObject@mu1.BART$a,
      b =         theObject@mu1.BART$b,
      augment =   theObject@mu1.BART$augment,
      rho =       theObject@mu1.BART$rho,
      usequants = theObject@mu1.BART$usequants,
      cont =      theObject@mu1.BART$cont,
      sigest =    theObject@mu1.BART$sigest,
      sigdf =     theObject@mu1.BART$sigdf,
      sigquant =  theObject@mu1.BART$sigquant,
      k =         theObject@mu1.BART$k,
      power =     theObject@mu1.BART$power,
      base =      theObject@mu1.BART$base,
      sigmaf =    theObject@mu1.BART$sigmaf,
      lambda =    theObject@mu1.BART$lambda,
      numcut =    theObject@mu1.BART$numcut,
      nskip =     theObject@mu1.BART$nskip
    )$yhat.test
    
    mu_hat_1 <- apply(pred_matrix_f_1, 2, mean)

    predictions <- mu_hat_1 - mu_hat_0
    
    if (return_CI) {
      get_CI_0 <- t(apply(pred_matrix_f_0, 2, function(x)
        quantile(x, probs = c(.05, 0.95))))
      get_CI_1 <- t(apply(pred_matrix_f_1, 2, function(x)
        quantile(x, probs = c(.05, 0.95))))
      
      CI_comb <- get_CI_1 - get_CI_0[, 2:1]
      
      return(cbind(predictions, CI_comb))
    } else{
      return(predictions)
    }
  }
)


#' CateCI-T_BART
#' @rdname CateCI-T_BART
#' @aliases CateCI,T_BART-method
#' @exportMethod CateCI
setMethod(
  f = "CateCI",
  signature = "T_BART",
  definition = function(theObject,
                        feature_new,
                        verbose = FALSE)
  {
    ndpost <- theObject@ndpost
    yobs <- theObject@yobs_train
    if ((is.null(yobs)) | (is.data.frame(yobs) && nrow(yobs) == 0)) {
      stop('Observations cannot be empty or null.')
    }
    feat <- theObject@feature_train
    if ((is.null(feat)) | (is.data.frame(feat) && nrow(feat) == 0)) {
      stop('Features cannot be empty or null.')
    }
    tr <- theObject@tr_train
    if ((is.null(tr)) | (is.data.frame(tr) && nrow(tr) == 0)) {
      stop('Treatment assignments cannot be empty or null.')
    }
    nthread <- theObject@nthread
    if (nthread > parallel::detectCores()) {
      warning('nthread is chosen bigger than the number of cores. It is changed
              to be equal to the number of cores.')
      nthread <- parallel::detectCores()
    }
    
    yobs_0 <- yobs[tr == 0]
    if (is.data.frame(yobs_0) && nrow(yobs_0) == 0) {
      stop('There is no observation in the control group (labelled 0).')
    }
    X_0 <- feat[tr == 0,]
    if (is.data.frame(X_0) && nrow(X_0) == 0) {
      stop('There is no feature in the control group (labelled 0).')
    }
    yobs_1 <- yobs[tr == 1]
    if (is.data.frame(yobs_1) && nrow(yobs_1) == 0) {
      stop('There is no observation in the treatment group (labelled 1).')
    }
    X_1 <- feat[tr == 1,]
    if (is.data.frame(X_1) && nrow(X_1) == 0) {
      stop('There is no feature in the treatment group (labelled 1).')
    }                
    # will contain for each estimator the prediciton for x_to_predict:
    output <- list()
    # will contain the appropraite CI:
    CI <- list()
    for (this_learner in c("l_first_0",
                           "l_first_1")) {
      if (this_learner == "l_first_0") {
        x = X_0
        y = yobs_0
        x_to_predict = feature_new
        parameters = theObject@mu0.BART
      } else if (this_learner == "l_first_1") {
        x = X_1
        y = yobs_1
        x_to_predict = feature_new
        parameters = theObject@mu1.BART
      }
      
      # x_to_predict_subset <- split(x_to_predict, rep(x = 1:nthread,
      #                                                times = get_CV_sizes(
      #                                                  nrow(x_to_predict),
      #                                                 nthread)))
      
      pred_matrix <- BART::mc.wbart(
        x.train = x,
        y.train = y,
        x.test =  x_to_predict,
        ndpost = ndpost,
        ntree = theObject@ntree,
        mc.cores = nthread,
        sparse =    parameters$sparse,
        theta =     parameters$theta,
        omega =     parameters$omega,
        a =         parameters$a,
        b =         parameters$b,
        augment =   parameters$augment,
        rho =       parameters$rho,
        usequants = parameters$usequants,
        cont =      parameters$cont,
        sigest =    parameters$sigest,
        sigdf =     parameters$sigdf,
        sigquant =  parameters$sigquant,
        k =         parameters$k,
        power =     parameters$power,
        base =      parameters$base,
        sigmaf =    parameters$sigmaf,
        lambda =    parameters$lambda,
        numcut =    parameters$numcut,
        nskip =     parameters$nskip
      )$yhat.test
      

      output[[this_learner]] <- apply(pred_matrix, 2, mean)
      CI[[this_learner]] <-
        t(apply(pred_matrix, 2, function(x)
          quantile(x, probs = c(.05, 0.95))))
      
    }
    if (verbose) {
      print("Done with the propensity score estimation.")
    }
    pred <- output[["l_first_1"]] - output[["l_first_0"]]
    
    CI_comb <- CI[["l_first_1"]] - CI[["l_first_0"]][, 2:1]
    
    to_return <- as.data.frame(cbind(pred, CI_comb))
    row.names(to_return) <- 1:nrow(to_return)
    colnames(to_return) <- c('pred','X5.','X95.')
    return(to_return)
  }
)
