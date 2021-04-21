# this file gets the best tunring parameters from our simulations

################################################################################
# READ IN MSE AND CI DATA --> file_list_MSE, file_list_CI
################################################################################
library(dplyr)
## read in the data as list for MSE
file_folder_MSE <- "Tune/XLearner/TuneOurDataSets/sim_data/"

file_list_MSE <- list()
i <- 1
for(file in dir(file_folder_MSE)){
  print(file)
  if(substr(file, 1,3) == "old") next
  file_list_MSE[[i]] <- tbl_df(read.csv(paste0(file_folder_MSE, file), as.is = TRUE))
  i <- i + 1
}


################################################################################
# SAVE BEST SETTING IN --> starting_values
################################################################################

## output for each file the two best settings for CI and MSE
best_setups <- data.frame()

for (file in file_list_MSE) {
  # file = file_list_MSE[[1]]
  file <- file[order(as.numeric(as.data.frame(file)[, 25])),]
  file$setup <- paste0(file$setup, "_MSE")
  best_setups <- rbind(best_setups, file[1, -25])
}

starting_values <- as.data.frame(best_setups)
devtools::use_data(starting_values, internal = TRUE, overwrite = TRUE)
# This will save the data in R/sysdata.rda and will only be available for our
# function



# Find CI default setting ------------------------------------------------------
tuning_settings <- file_list_CI[[1]][,-c(1,25,26)] # tuning_setting = all settings
mergingvars <- names(tuning_settings)
for(file in file_list_CI){
  print(nrow(file))
  file <- file[ ,-1]
  tuning_settings <- merge(tuning_settings, file, by = mergingvars)
}
tuning_settings <- tbl_df(tuning_settings)


names(tuning_settings)
results <- tuning_settings[ , 24:45]
perc_na <- apply(results, 2, function(x) mean(is.na(x)))
results_selected <- results[ , perc_na <.25]
results_selected_CI <- results_selected[ , (1:7)*2]

table(apply(results_selected_CI >.6,1, mean))

# this has the best coverage:
good_coverage <- which(apply(results_selected_CI >.6,1, all))
as.numeric(results_selected_CI[good_coverage, ])


# Find best in terms of MSE:

results_selected_MSE <- results_selected[ , (1:7)*2 - 1]
apply(results_selected_MSE, 2, function(x) quantile(x, probs = c(.05, .1, .2), na.rm = TRUE))


results_selected_MSE[good_coverage, ]

good_MSE <- which.min(apply(apply(results_selected_MSE, 2, function(x) rank(x, na.last = TRUE)), 1, sum))

results_selected_MSE[good_MSE, ]
as.numeric(results_selected_CI[good_MSE, ])

as.data.frame(tuning_settings[good_MSE, 1:23]) #taken
################################################################################
# FIND GLOBALLY BEST STARTING SETTINGS:
################################################################################
## see which setting is similar to which by looking at correlations:
tuning_settings <- file_list[[1]][,-c(1,25)] # tuning_setting = all settings
mergingvars <- names(tuning_settings)
for(file in file_list){
  file <- file[ ,-1]
  tuning_settings <- merge(tuning_settings, file, by = mergingvars)
}

## look at pairwise correlations:
library(gpairs)
tuning_settings$STMpp[tuning_settings$STMpp > 400] <- NA
if(FALSE) gpairs(tuning_settings[ ,15:25])

# important groups are : "complexTau", "Conf1", "rare1", "STMpp", "Wager1", "Wager3"

## Select good features:
# devtools::install_github("rstudio/crosstalk")
# devtools::install_github("jcheng5/d3scatter")
# devtools::install_github("rstudio/leaflet")

# Example
library(crosstalk)
library(d3scatter)

shared_iris <- SharedData$new(iris)
bscols(
  d3scatter(shared_iris, ~Petal.Length, ~Petal.Width, width="100%", height=300),
  d3scatter(shared_iris, ~Sepal.Length, ~Sepal.Width, width="100%", height=300)
)

# find for each all setups which is under 10% best
quantiles10 <- apply(tuning_settings[ ,15:25],
                     2,
                     function(x) quantile(x, probs = c(.1), na.rm = TRUE))

apply(t(t(tuning_settings[ ,15:25]) < quantiles10), 2, function(x) mean(x, na.rm = TRUE))

howoftenbest <- apply(t(t(tuning_settings[ ,15:25]) < quantiles10), 1, function(x) sum(x, na.rm = TRUE))
ss <- cbind(t(t(tuning_settings[ ,15:25]) < quantiles10), howoftenbest)
relevant_high <- ss[,12] > 7
table(relevant_high)
apply(ss[relevant_high,], 2, function(x) table(x))



shared_tuning_settings <- SharedData$new(tuning_settings[relevant_high, ])
bscols(
  d3scatter(shared_tuning_settings, ~Wager3, ~Wager2, width="100%", height=300 ),
  d3scatter(shared_tuning_settings, ~Usual1, ~Ufail, width="100%", height=300),
  d3scatter(shared_tuning_settings, ~RsparseT2weak, ~RespSparseTau1strong, width="100%", height=300),
  d3scatter(shared_tuning_settings, ~Conf1, ~complexTau, width="100%", height=300)
)

shared_tuning_settings_all <- SharedData$new(cbind(tuning_settings, aa = factor(relevant_high, levels = c("TRUE", "FALSE"))))
bscols(
  d3scatter(shared_tuning_settings_all, ~rare1, ~STMpp, ~aa, width="100%", height=300 ),
  d3scatter(shared_tuning_settings_all, ~Wager1, ~STMpp,~aa, width="100%", height=300 )
)

#starting point 1
apply(tuning_settings[ , 1:14] , 2, function(x) table(x, relevant_high))


sp_1_firstStage <- data.frame(
  mtry_first = round(ncol(feat) / 2),
  min_node_size_spl_first = 1,
  min_node_size_ave_first = 5,
  splitratio_first = .5,
  replace_first = TRUE,
  sample_fraction_first = 0.8
)
sp_1_secondStage <- data.frame(
  predmode = "propmean",
  mtry_second = ncol(feat),
  min_node_size_spl_second = 5,
  min_node_size_ave_second = 3,
  splitratio_second = .5,
  replace_second = TRUE,
  sample_fraction_second = 0.9
)

#starting point 2

ss <- as.data.frame(ss)
ss$rare1
relevant_high2 <- !is.na(ss$rare1) & (ss$rare1 == 1 | ss$STMpp == 1 | ss$Wager1 == 1)

apply(tuning_settings[ , 1:14] , 2, function(x) table(x, relevant_high2))


sp_2_firstStage <- data.frame(
    mtry_first = round(ncol(feat) / 2),
    min_node_size_spl_first = 1,
    min_node_size_ave_first = 1,
    splitratio_first = .4,
    replace_first = FALSE,
    sample_fraction_first = 0.7
  )
sp_2_secondStage <- data.frame(
  predmode = "propmean",
  mtry_second = max(1, round(ncol(feat) / 5)),
  min_node_size_spl_second = 10,
  min_node_size_ave_second = 3,
  splitratio_second = .75,
  replace_second = FALSE,
  sample_fraction_second = 0.8
)




