# This file is used to tune X, T, and S learner as well as possible.
library(hte)

# We start by loading the setup and
args = (commandArgs(TRUE))
print(args)
setup_i <- -as.numeric(args)
print(setup_i)

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

#setup <- "RespSparseTau1strong"

set.seed(59906)
nthread <- 24
nsamples <- 10000
dim <- 20
ntrain <- 10000
ntest <- 10000
alpha <- .01

# Sample each parameter
Rand_tune <- data.frame(
  setup = setup,
  dim = dim,
  ntrain = ntrain,
  mtry = sample(1:dim, nsamples, replace = TRUE),
  nodesizeSpl = sample(c(1, 3, 5, 10, 30, 100), nsamples, replace = TRUE),
  nodesizeAvg = sample(c(1, 3, 5, 10, 30, 100), nsamples, replace = TRUE),
  splitratio = sample(seq(from = 0.1, to = 1, by = 0.1), nsamples, replace = TRUE),
  replace = sample(c(TRUE, FALSE), nsamples, replace = TRUE),
  sample_fraction = sample(seq(from = 0.1, to = 1, by = 0.1), nsamples,
                                 replace = TRUE)
)


## Setting up where the data should be saved:
data_folder_name <- "sim_data/"
if (!dir.exists(data_folder_name))
  dir.create(data_folder_name)
filename <-
  paste0(data_folder_name, "Rand_tune_", setup, "_TRF_", Sys.Date(), ".csv")
# if (file.exists(filename))
#   file.remove(filename)

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
  EMSE <-
    tryCatch({
      L <- T_RF(
        feat = experiment$feat_tr,
        tr = experiment$W_tr,
        yobs = experiment$Yobs_tr,
        mtry = Rand_tune$mtry[i],
        nodesizeSpl = Rand_tune$nodesizeSpl[i],
        nodesizeAvg = Rand_tune$nodesizeAvg[i],
        replace = Rand_tune$replace[i],
        ntree = 500,
        sample_fraction = Rand_tune$sample_fraction[i],
        nthread = nthread,
        splitratio = Rand_tune$splitratio[i]
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
