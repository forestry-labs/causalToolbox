# this file gets the best tunring parameters from our simulations

library(dplyr)
## read in the data as list:
file_folder <- "Tune/XLearner/sim_data/"
file_list <- list()
i <- 1
for(file in dir(file_folder)){
  print(file)
  file_list[[i]] <- tbl_df(read.csv(paste0(file_folder, file), as.is = TRUE))
  i <- i + 1
}

## see which setting is similar to which by looking at correlations:
tuning_settings <- file_list[[1]][,-c(1,16)] # tuning_setting = all settings
mergingvars <- names(tuning_settings)
for(file in file_list){
  file <- file[ ,-1]
  tuning_settings <- merge(tuning_settings, file, by = mergingvars)
}

## look at pairwise correlations:
library(gpairs)
tuning_settings$STMpp[tuning_settings$STMpp >400] <- NA
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


apply(tuning_settings[ , 1:14] , 2, function(x) table(x, relevant_high))









