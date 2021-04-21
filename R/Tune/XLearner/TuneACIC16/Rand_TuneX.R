# This file is used to tune X, T, and S learner as well as possible.
library(hte)

# We start by loading the setup and
args = (commandArgs(TRUE))
print(args)
setup_i <- -as.numeric(args)
print(setup_i) # setup_i <- 4

setup_grid <- c(
    "RespSparseTau1strong",
    "RsparseT2weak",
    "complexTau",
    "Conf1",
    "rare1",
    "STMpp",
    "Ufail",
    "Usual1",
    "Wager1",
    "Wager2",
    "Wager3"
  )
setup <- setup_grid[setup_i]
print(setup)

set.seed(59906)
nthread <- 0
nsamples <- 1e5
dim <- 20
ntrain <- 5000
ntest <- 10000
alpha <- .01


# Sample each parameter
Rand_tune <- data.frame(
  setup = setup,
  dim = dim,
  ntrain = ntrain,
  replace_first  = sample(c(TRUE, FALSE), nsamples, replace = TRUE),
  replace_second = sample(c(TRUE, FALSE), nsamples, replace = TRUE),
  replace_prop   = sample(c(TRUE, FALSE), nsamples, replace = TRUE),
  sample.fraction_first  = sample(seq(from = 0.1, to = 1, by = 0.1), nsamples,replace = TRUE),
  sample.fraction_second = sample(seq(from = 0.1, to = 1, by = 0.1), nsamples,replace = TRUE),
  sample.fraction_prop   = sample(seq(from = 0.1, to = 1, by = 0.1), nsamples,replace = TRUE),
  mtry_first  = sample(1:dim, nsamples, replace = TRUE),
  mtry_second = sample(1:dim, nsamples, replace = TRUE),
  mtry_prop   = sample(1:dim, nsamples, replace = TRUE),
  nodesizeSpl_first  = round(rbeta(nsamples, 1, 4) * 100) + 1,
  nodesizeSpl_second = round(rbeta(nsamples, 1, 4) * 100) + 1,
  nodesizeSpl_prop   = round(rbeta(nsamples, 1, 4) * 100) + 1,
  nodesizeAvg_first  = round(rbeta(nsamples, 1, 4) * 100) + 1,
  nodesizeAvg_second = round(rbeta(nsamples, 1, 4) * 100) + 1,
  nodesizeAvg_prop   = round(rbeta(nsamples, 1, 4) * 100) + 1,
  splitratio_first  = sample(seq(from = 0.1, to = 1, by = 0.1), nsamples, replace = TRUE),
  splitratio_second = sample(seq(from = 0.1, to = 1, by = 0.1), nsamples, replace = TRUE),
  splitratio_prop   = sample(seq(from = 0.1, to = 1, by = 0.1), nsamples, replace = TRUE),
  middleSplit_first  = sample(c(TRUE, FALSE), nsamples, replace = TRUE),
  middleSplit_second = sample(c(TRUE, FALSE), nsamples, replace = TRUE),
  middleSplit_prop   = sample(c(TRUE, FALSE), nsamples, replace = TRUE)
)


## Setting up where the data should be saved:
data_folder_name <- "sim_data/"
if (!dir.exists(data_folder_name))
  dir.create(data_folder_name)
filename <-
  paste0(data_folder_name, "Rand_tune_", setup, "_XRF_", packageVersion("hte"),
         "_", Sys.Date(), ".csv")

# loop through rows and fill in the NA columns:
for (i in 1:nsamples) {
  print(paste("Starting with", i, "of", nsamples))
  stat.time <- Sys.time()
  # run the experiment:
  experiment <- simulate_causal_experiment(
    ntrain = ntrain,
    ntest = ntest,
    dim = dim,
    alpha = alpha,
    feat_distribution = "normal",
    setup = setup,
    testseed = 293901,
    trainseed = 93007
  )

  set.seed(78136)
  # train the X-learner:
  general_hyperpara <- list("predmode" = "propmean",
                            "nthread" = 0)

  first_stage_hyperpara <- list(
    "relevant_Variable" = 1:dim,
    "ntree" = 1000,
    "replace" = Rand_tune$replace_first[i],
    "sample.fraction" = Rand_tune$sample.fraction_first[i],
    "mtry" = Rand_tune$mtry_first[i],
    "nodesizeSpl" = Rand_tune$nodesizeSpl_first[i],
    "nodesizeAvg" = Rand_tune$nodesizeAvg_first[i],
    "splitratio" = Rand_tune$splitratio_first[i],
    "middleSplit" = Rand_tune$middleSplit_first[i]
  )
  second_stage_hyperpara <- list(
    "relevant_Variable" = 1:dim,
    "ntree" = 1000,
    "replace" = Rand_tune$replace_second[i],
    "sample.fraction" = Rand_tune$sample.fraction_second[i],
    "mtry" = Rand_tune$mtry_second[i],
    "nodesizeSpl" = Rand_tune$nodesizeSpl_second[i],
    "nodesizeAvg" = Rand_tune$nodesizeAvg_second[i],
    "splitratio" = Rand_tune$splitratio_second[i],
    "middleSplit" = Rand_tune$middleSplit_second[i]
  )
  prop_hyperpara <- list(
    "relevant_Variable" = 1:dim,
    "ntree" = 1000,
    "replace" = Rand_tune$replace_prop[i],
    "sample.fraction" = Rand_tune$sample.fraction_prop[i],
    "mtry" = Rand_tune$mtry_prop[i],
    "nodesizeSpl" = Rand_tune$nodesizeSpl_prop[i],
    "nodesizeAvg" = Rand_tune$nodesizeAvg_prop[i],
    "splitratio" = Rand_tune$splitratio_prop[i],
    "middleSplit" = Rand_tune$middleSplit_prop[i]
  )

  hyperparameter_list <- list(
    "general" = general_hyperpara,
    "l_first_0" = first_stage_hyperpara,
    "l_first_1" = first_stage_hyperpara,
    "l_second_0" = second_stage_hyperpara,
    "l_second_1" = second_stage_hyperpara,
    "l_prop" = prop_hyperpara
  )

  EMSE <-
    tryCatch({
      L <- X_RF_fully_specified(
        feat = experiment$feat_tr,
        tr = experiment$W_tr,
        yobs = experiment$Yobs_tr,
        hyperparameter_list = hyperparameter_list,
        verbose = FALSE
        )
      mean((experiment$tau_te - EstimateCate(L, experiment$feat_te)) ^ 2)
    },
    error = function(e) {
      print(e)
      NA
    })

  Rand_tune[i, setup] <- EMSE
  print(paste(
    "   Done with",
    setup,
    "Error was:",
    EMSE,
    "   and it took",
    as.numeric(Sys.time() - stat.time, units = "mins"),
    " minutes."
  ))

  col.names <- !file.exists(filename)
  write.table(
    Rand_tune[i, ],
    file = filename,
    append = TRUE,
    col.names = col.names,
    row.names = FALSE,
    sep = ","
  )
}
