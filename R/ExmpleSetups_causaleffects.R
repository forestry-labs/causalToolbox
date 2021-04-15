# This file implements different heterogeneous treatment effect estimation
# problems

# Correlation Matrix Simulator -------------------------------------------------
#' @rdname causalExp
#' @name simulate Experiments
#' @description \code{simulate_correlation_matrix} uses the C-vine method for
#'   simulating correlation matrices. (Refer to the referenced paper for details.)
#' @return A correlation matrix.
simulate_correlation_matrix <- function(dim, alpha) {

  betaparam <- 1 / alpha

  P <- matrix(nrow = dim, ncol = dim)           # storing partial correlations
  S <- diag(dim)

  for (k in 1:(dim - 1)) {
    for (i in (k + 1):dim) {
      P[k, i] <- rbeta(1, betaparam, betaparam) # sampling from beta
      P[k, i] <- (P[k, i] - 0.5) * 2     # linearly shifting to [-1, 1]
      p <- P[k, i]
      if (k > 1) {
        for (l in (k - 1):1) {
          # converting partial correlation to raw correlation
          p <- p * sqrt((1 - P[l, i] ^ 2) * (1 - P[l, k] ^ 2)) + P[l, i] *
            P[l, k]
          # p  = p * sqrt((1-P(l,i)^2)*(1-P(l,k)^2)) + P(l,i)*P(l,k);
        }
      }
      S[k, i] <- p
      S[i, k] <- p
    }
  }

  # permuting the variables to make the distribution permutation-invariant
  permutation <- sample(1:dim)
  S <- S[permutation, permutation]
  return(S)
}


# Causal Experiment Simulator --------------------------------------------------
#' @title Simulate a Causal Experiment
#' @rdname causalExp
#' @description \code{simulate_causal_experiment} simulates an RCT or observational data for causal
#'   effect estimation. It is mainly used to test different heterogenuous
#'   treatment effect estimation strategies.
#' @param ntrain Number of training examples.
#' @param ntest Number of test examples.
#' @param dim Dimension of the data set.
#' @param alpha Only used if \code{given_features} is not set and \code{feat_distribution} is
#'   chosen to be normal. It specifies how correlated the features can be. If
#'   alpha = 0, then the features are independent. If alpha is very large, then
#'   the features can be very correlated. Use the
#'   \code{simulate_correlation_matrix} function to get a better understanding
#'   of the impact of alpha.
#' @param feat_distribution Only used if \code{given_features} is not specified. Either
#'   "normal" or "unif." It specifies the distribution of the features.
#' @param given_features This is used if we already have features and want to
#'   test the performance of different estimators for a particular set of
#'   features.
#' @param pscore,mu0,tau Parameters that determine the propensity score, the
#'   response function for the control units, and tau, respectively. The
#'   different options can be seen using
#'   \code{names(pscores.simulate_causal_experiment)},
#'   \code{names(mu0.simulate_causal_experiment)}, and
#'   \code{names(tau.simulate_causal_experiment)}. This is implemented in this manner,
#'   because it enables the user to easily loop through the different
#'   estimators.
#' @param testseed The seed used to generate the test data. If NULL, then
#'   the seed of the main session is used.
#' @param trainseed The seed used to generate the training data. If NULL,
#'   then the seed of the main session is used.
#' @return A list with the following elements:
#'   \item{\code{setup_name}}{Name of the setup.}
#'   \item{\code{m_t_truth}}{Function containing the response function of the
#'   treated units.}
#'   \item{\code{m_c_truth}}{Function containing the response function of the
#'   control units.}
#'   \item{\code{propscore}}{Propensity score function.}
#'   \item{\code{alpha}}{Chosen alpha.}
#'   \item{\code{feat_te}}{Data.frame containing the features of the test
#'   samples.}
#'   \item{\code{W_te}}{Numeric vector containing the treatment assignment of
#'   the test samples.}
#'   \item{\code{tau_te}}{Numeric vector containing the true conditional average
#'   treatment effects of the test samples.}
#'   \item{\code{Yobs_te}}{Numeric vector containing the observed Y values of
#'   the test samples.}
#'   \item{\code{feat_tr}}{Data.frame containing the features of the training
#'   samples.}
#'   \item{\code{W_tr}}{Numeric vector containing the treatment assignment of
#'   the training samples.}
#'   \item{\code{tau_tr}}{Numeric vector containing the true conditional average
#'   treatment effects of the training samples.}
#'   \item{\code{Yobs_tr}}{Numeric vector containing the observed Y values of
#'   the training samples.}
#' @details The function simulates causal experiments by generating the
#'   features, treatment assignment, observed Y values, and CATE for
#'   a test set and a training set. pscore, mu0, and tau define the response
#'   functions and the propensity score. For example, \code{pscore =
#'   "osSparse1Linear"} specifies that \deqn{e(x) = max(0.05, min(.95, x1 / 2 +
#'   1 / 4))} and \code{mu0 ="sparseLinearWeak"} specifies that the response
#'   function for the control units is given by the simple linear function,
#'   \deqn{mu0(x) = 3 x1 + 5 x2.}
#' @seealso \code{\link{X-Learner}}
#' @references \itemize{
#'   \item Daniel Lewandowskia, Dorota Kurowickaa, Harry Joe (2009). Generating
#'   Random Correlation Matrices Based on Vines and Extended Onion Method.
#'   \item Sören Künzel, Jasjeet Sekhon, Peter Bickel, and Bin Yu (2017).
#'     Meta-learners for Estimating Heterogeneous Treatment Effects Using
#'     Machine Learning.
#'     }
#' @examples
#' require(causalToolbox)
#'
#' ce_sim <- simulate_causal_experiment(
#'   ntrain = 20,
#'   ntest = 20,
#'   dim = 7
#' )
#'
#' ce_sim
#'
#' \dontrun{
#' estimators <- list(
#'   S_RF = S_RF,
#'   T_RF = T_RF,
#'   X_RF = X_RF,
#'   S_BART = S_BART,
#'   T_BART = T_BART,
#'   X_BARTT = X_BART)
#'
#' performance <- data.frame()
#' for(tau_n in names(tau.simulate_causal_experiment)){
#'   for(mu0_n in names(mu0.simulate_causal_experiment)) {
#'     ce <- simulate_causal_experiment(
#'       given_features = iris,
#'       pscore = "rct5",
#'       mu0 = mu0_n,
#'       tau = tau_n)
#'
#'     for(estimator_n in names(estimators)) {
#'       print(paste(tau_n, mu0_n, estimator_n))
#'
#'       trained_e <- estimators[[estimator_n]](ce$feat_tr, ce$W_tr, ce$Yobs_tr)
#'       performance <-
#'         rbind(performance,
#'               data.frame(
#'                 mu0 = mu0_n,
#'                 tau = tau_n,
#'                 estimator = estimator_n,
#'                 MSE = mean((EstimateCate(trained_e, ce$feat_te) -
#'                             ce$tau_te)^2)))
#'     }
#'   }
#' }
#'
#' reshape2::dcast(data = performance, mu0 + tau ~ estimator)
#' }
#' @import MASS
#' @export
simulate_causal_experiment <- function(ntrain = nrow(given_features),
                                       ntest = nrow(given_features),
                                       dim = ncol(given_features),
                                       alpha = .1,
                                       feat_distribution = "normal",
                                       given_features = NULL,
                                       pscore = "rct5",
                                       mu0 = "sparseLinearStrong",
                                       tau = "sparseLinearWeak",
                                       testseed = NULL,
                                       trainseed = NULL,
                                       noEffect = FALSE) {

  ## First we define a base function which will later be called with different
  # setups:
  # the following function is the base function when creating the features.
  # the different setups below specify the specific form of the potential
  # outcomes and the treatment assignment and then call this function
  createTrainAndTest_base <-
    function(ntrain,
             ntest,
             dim,
             m_t_truth,
             m_c_truth,
             propscore,
             alpha,
             feat_distribution,
             given_features) {
      # this is the base creator, which is called by all other setups:
      tau <- function(feat) {
        as.numeric(m_t_truth(feat) - m_c_truth(feat))
      }
      if (is.null(given_features)) {
        given_features_fkt <- function(n, dim) {
          if (feat_distribution == "normal") {
            correlation_mat <- simulate_correlation_matrix(dim, alpha = alpha)
            mu <- rep(0, dim)
            feat <-
              data.frame(MASS::mvrnorm(
                n = n,
                mu = mu,
                Sigma = correlation_mat
              ))
          }
          if (feat_distribution == "unif") {
            feat <- data.frame(matrix(runif(n * dim), nrow = n))
          }
          colnames(feat) <- paste0("x", 1:(dim))
          return(feat)
        }
      } else{
        given_features_fkt <- function(n, dim) {
          colnames(given_features) <- paste0("x", 1:(ncol(given_features)))
          return(given_features[1:n, ])
        }
      }

      getW <- function(feat) {
        nn <- nrow(feat)
        propS <- propscore(feat)
        return(rbinom(nn, 1, propS))
      }
      getYobs <- function(feat, W) {
        nn <- nrow(feat)
        return(ifelse(
          W == 1,
          m_t_truth(feat) + rnorm(nn, 0, 1),
          # treated
          m_c_truth(feat) + rnorm(nn, 0, 1)   # control
        ))
      }

      ## fixing the training and testing seed (if given) without changing the
      # seed of the seesion:
      # first save the current random seed, then set the seed to what we like,
      # and then set the seed back to what it was.

      if (!is.null(trainseed)) {
        current_seed <- .Random.seed  # saves the current random stage
        set.seed(trainseed)           # introduces a new seed to stay consistent
      }
      feat_tr <- given_features_fkt(ntrain, dim)
      feat_tr_numeric <- as.data.frame(sapply(feat_tr, function(x) {
        if (is.numeric(x)) {
          return(x)
        } else{
          as.numeric(as.factor(x))
        }
      }))

      W_tr <- getW(feat_tr_numeric)
      Yobs_tr <- getYobs(feat_tr_numeric, W_tr)
      if (!is.null(trainseed)) {
        .Random.seed <- current_seed  # sets back the current random stage
      }
      if (!is.null(testseed)) {
        current_seed <- .Random.seed  # saves the current random stage
        set.seed(testseed)          # introduces a new seed to stay consistent
      }
      feat_te <- given_features_fkt(ntest, dim)
      feat_te_numeric <- as.data.frame(sapply(feat_te, function(x) {
        if (is.numeric(x)) {
          return(x)
        } else{
          as.numeric(as.factor(x))
        }
      }))
      W_te <- getW(feat_te_numeric)
      Yobs_te <- getYobs(feat_te_numeric, W_te)
      if (!is.null(testseed)) {
        .Random.seed <- current_seed  # sets back the current random stage
      }
      return(
        list(
          alpha = alpha,
          feat_te = feat_te,
          W_te = W_te,
          tau_te = tau(feat_te_numeric),
          Yobs_te = Yobs_te,
          feat_tr = feat_tr,
          W_tr = W_tr,
          tau_tr = tau(feat_tr_numeric),
          Yobs_tr = Yobs_tr
        )
      )
    }

  # Definition of different setups ---------------------------------------------
  if(!mu0 %in% names(mu0.simulate_causal_experiment)){
    stop(paste("mu0 must be one of",
               paste(names(mu0.simulate_causal_experiment), collapse = ", ")))
  }
  if(!tau %in% names(tau.simulate_causal_experiment)){
    stop(paste("tau must be one of",
               paste(names(tau.simulate_causal_experiment), collapse = ", ")))
  }
  if(!pscore %in% names(pscores.simulate_causal_experiment)){
    stop(paste("pscore must be one of",
               paste(names(pscores.simulate_causal_experiment), collapse = ", ")))
  }

  m_c_truth <- mu0.simulate_causal_experiment[[mu0]]
  if (noEffect) {
    m_t_truth <- function(feat) m_c_truth(feat)
  } else {
    m_t_truth <- function(feat)
      m_c_truth(feat) + tau.simulate_causal_experiment[[tau]](feat)
  }
  propscore <- pscores.simulate_causal_experiment[[pscore]]

  return(
    c(list(
        setup_name = paste0("mu0=", mu0,", tau=", tau, ", pscore=", pscore),
        m_t_truth = m_t_truth,
        m_c_truth = m_c_truth,
        propscore = propscore),
      createTrainAndTest_base(
        ntrain,
        ntest,
        dim,
        m_t_truth,
        m_c_truth,
        propscore,
        alpha,
        feat_distribution,
        given_features)))
}

# Propensity score functions ---------------------------------------------------
#' @rdname causalExp
#' @format NULL
#' @export
pscores.simulate_causal_experiment <- list(
  rct5 = function(feat) {.5},
  rct1 = function(feat) {.1},
  rct01 = function(feat) {.01},
  osSparse1Linear = function(feat) {
    apply(feat, 1, function(x)
      max(0.05, min(.95, x[1] / 2 + 1 / 4)))
  },
  osSparse1Beta = function(feat) {0.25 + dbeta(feat$x1, 2, 4) / 4})

# mu0 functions ----------------------------------------------------------------
#' @rdname causalExp
#' @format NULL
#' @export
mu0.simulate_causal_experiment <- list(
  simple = function(feat) {2 * feat$x1 - 1},
  sparseLinearWeak = function(feat) {3 * feat$x1 + 5 * feat$x2},
  sparseLinearStrong = function(feat) {30 * feat$x1 + 50 * feat$x2},
  fullLinearWeak = function(feat) {
    oldSeed <- .Random.seed; on.exit( {.Random.seed <<- oldSeed} )
    set.seed(53979361)
    d <- ncol(feat)

    beta <- runif(d, -5, 5)
    as.matrix(feat) %*% beta
  },
  semiLinear = function(feat) {
    oldSeed <- .Random.seed; on.exit( {.Random.seed <<- oldSeed} )
    set.seed(53979361)
    d <- ncol(feat)

    beta <- runif(d, -5, 5)
    ifelse(feat[,1] > .5,
           as.matrix(feat) %*% beta + 5,
           as.matrix(feat) %*% beta)
  },
  fullLinearStrong = function(feat) {
    oldSeed <- .Random.seed; on.exit({.Random.seed <<- oldSeed})
    set.seed(53979361)
    d <- ncol(feat)

    beta <- runif(d, 1, 30)
    as.matrix(feat) %*% beta
  },
  fullLocallyLinear = function(feat) {
    oldSeed <- .Random.seed; on.exit({.Random.seed <<- oldSeed})
    set.seed(7020829)
    d <- ncol(feat)

    beta1 <- runif(d, -15, 15)
    beta2 <- runif(d, -15, 15)
    beta3 <- runif(d, -15, 15)

    ifelse(
      feat[, ncol(feat)] < -0.4,
      as.matrix(feat) %*% beta1,
      ifelse(
        feat[, ncol(feat)] < 0.4,
        as.matrix(feat) %*% beta2,
        as.matrix(feat) %*% beta3))
  },
  fullLinearWeakStep = function(feat) {
    oldSeed <- .Random.seed; on.exit({.Random.seed <<- oldSeed})
    set.seed(1496661)
    d <- ncol(feat)

    beta <- runif(d, -5, 5)
    as.matrix(feat) %*% beta + ifelse(feat$x1 > 0, 5, 0)
  },
  sparseNonLinear1 = function(feat) {
    sin(feat$x1) *
      sin(feat$x2) *
      sin(feat$x3) *
      sin(feat$x4)
  },
  sparseNonLinear2 = function(feat) {
    4 / ((1 + exp(-12 * (feat$x1 - 0.5))) *
           (1 + exp(-12 * (feat$x2 - 0.5))) *
           (1 + exp(-12 * (feat$x3 - 0.5))) *
           (1 + exp(-12 * (feat$x4 - 0.5))) *
           (1 + exp(-12 * (feat$x5 - 0.5))))
  },
  sparseNonLinear3 = function(feat) {
    (1 + 1 / (1 + exp(-20 * (feat$x1 - 1 / 3)))) *
      (1 + 1 / (1 + exp(-20 * (feat$x2 - 1 / 3))))
  },
  complexNonLinear = function(feat) {
    -2/((1 + exp(-12 * (feat$x1 - 0.5)))*
         (1 + exp(-12 * (feat$x2 - 0.5))))
  }
)

# tau functions ----------------------------------------------------------------
#' @rdname causalExp
#' @format NULL
#' @export
tau.simulate_causal_experiment <- list(
  no = function(feat) {0},
  const = function(feat) {10},
  sparseLinearWeak = function(feat) {3 * feat$x1 + 5 * feat$x2},
  fullLinearWeak = function(feat) {
    oldSeed <- .Random.seed; on.exit({.Random.seed <<- oldSeed})
    set.seed(53979361)
    d <- ncol(feat)

    beta <- runif(d, -5, 5)
    as.matrix(feat) %*% beta
  },
  sparseLinearWeak = function(feat) {3 * feat$x1 + 5 * feat$x2},
  fullLinearStrong = function(feat) {
    oldSeed <- .Random.seed; on.exit({.Random.seed <<- oldSeed})
    set.seed(53979361)
    d <- ncol(feat)

    beta <- runif(d, 1, 30)
    as.matrix(feat) %*% beta
  },
  semiLinear = function(feat) {
    oldSeed <- .Random.seed; on.exit( {.Random.seed <<- oldSeed} )
    set.seed(53979361)

    ifelse(feat[,2] > .5,
           8,
           0)
  },
  fullLocallyLinear = function(feat) {
    oldSeed <- .Random.seed; on.exit({.Random.seed <<- oldSeed})
    set.seed(6482481)
    d <- ncol(feat)

    beta1 <- runif(d, -15, 15)
    beta2 <- runif(d, -15, 15)
    beta3 <- runif(d, -15, 15)

    ifelse(
      feat[, ncol(feat)] < -0.4,
      as.matrix(feat) %*% beta1,
      ifelse(
        feat[, ncol(feat)] < 0.4,
        as.matrix(feat) %*% beta2,
        as.matrix(feat) %*% beta3))
  },
  sparseNonLinear3 = function(feat) {
    (1 + 1 / (1 + exp(-20 * (feat$x1 - 1 / 3)))) *
      (1 + 1 / (1 + exp(-20 * (feat$x2 - 1 / 3))))
  },
  complexNonLinear = function(feat) {
    2/((1 + exp(-12 * (feat$x1 - 0.5)))*
         (1 + exp(-12 * (feat$x2 - 0.5))))
  })

# Real data --------------------------------------------------------------------
#' @title Get Out To Vote
#' @description
#' This is an example data set, and it has been created by looking at a
#' certain subset of the "Social Pressure and Voter Turnout: Evidence from a
#' Large-Scale Field Experiment" study that tested the impact of social pressure
#' on voter turnout. A precise description and the full data set can be found at
#' \url{https://isps.yale.edu/research/data/d001}.
#'
#' The study consists of seven key, individual-level covariates, most of which
#' are discrete: gender, age, and whether the registered individual voted in the
#' primary elections in 2000, 2002, and 2004 or the general elections in 2000 and
#' 2002. The sample was restricted to voters who had voted in the 2004 general
#' election. The outcome of interest was the turnout in the 2006 primary election,
#' which was an indicator variable, and the treatment was whether or not the
#' subjects were elected to receive a mailer.
#' @references
#' \itemize{
#'   \item Gerber AS, Green DP, Larimer CW (2008) Social Pressure and Voter
#'   Turnout: Evidence from a Large-Scale Field Experiment. Am Polit Sci Rev
#'   102:33–48. \url{https://isps.yale.edu/research/publications/isps08-001}
#'  }
"gotv"

