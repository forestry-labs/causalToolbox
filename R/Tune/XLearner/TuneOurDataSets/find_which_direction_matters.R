library(dplyr)
## read in the data as list:
file_folder <- "Tune/XLearner/sim_data/"
file_list <- list()
i <- 1
for(file in dir(file_folder)){
  if(substr(file, nchar(file)-5, nchar(file) - 4) != 22) next
  print(file)
  file_list[[i]] <- tbl_df(read.csv(paste0(file_folder, file), as.is = TRUE))
  i <- i + 1
}

## see which setting is similar to which by looking at correlations:
tuning_settings <- file_list[[1]][,-c(1,25)] # tuning_setting = all settings
mergingvars <- names(tuning_settings)
for(file in file_list){
  file <- file[ ,-1]
  tuning_settings <- merge(tuning_settings, file, by = mergingvars)
}

tuning_settings %>% tbl_df()

is_rel <- c()
for(file in file_list){
  file <- file[ ,-1]
  file$outcome <- as.numeric(as.data.frame(file)[,24])
  file <- file[,-24]
  file <- as.data.frame(file)
  pvals <- summary(lm(outcome ~ ., data = file))$coef[,4]
  is_rel <- c(is_rel, names(pvals)[pvals < .001])
}
table(is_rel)
# mtry_first
# mtry_second
# nodesizeAvg_first
# nodesizeAvg_second
# nodesizeSpl_first
# nodesizeSpl_second
# sampsize_first
# sampsize_second
# splitratio_second



