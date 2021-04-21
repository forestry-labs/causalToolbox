#' @include helper_functions.R
# ------------------------------------------------------------------------------

#' gof_subset
#' @title gof_subset
#' @name gof_subset
#' @description This function estimates ATE for each subset, average the CATEs 
#' over each of the subsets and see how close it is to the truth. 
#' @param important.features names of features which should be used in 
#' quickmatch to find the relevant subgroups. Only features specified here will 
#' be used to create the subgroups
#' @param min.treat.size.per.group the minimum size of each treatment group 
#' per matched subset. For example, if it is chosen to be 25 (the default), then 
#' at least 25 units in each subgroup are treated and 25 are in the control
#' group
#' @param normalize Specifies how the distance should be normalized. 
#' "none": no normalization. 
#' "mahalanobize": normalization by var(data)
#' "studentize" (default): normalization is done with the diagonal of var(data)
#' @inheritParams gof_transformed
#' @return mean(error) and sd(error)
#' @import quickmatch
#' @import distances
#' @export gof_subset
gof_subset <- function(feat, yobs, tr, estimator, 
                       important.features = colnames(feat), 
                       min.treat.size.per.group = 25,
                       normalize = "studentize",
                       k = 5,
                       verbose = FALSE) {
  n_obs <- length(tr)
  
  # ----------------------------------------------------------------------------
  # Catch Errors
  catch_error(feat, yobs, tr, k)
  if (sum(tr == 1) < min.treat.size.per.group) {
    stop("The treatment size must be greater than the min.treat.size.per.group.")
  }
  
  # ----------------------------------------------------------------------------
  # get subgroups
  feat_distances <- distances::distances(feat, 
                                         dist_variables = important.features, 
                                         normalize = normalize)
  unit_match <- quickmatch::quickmatch(distances = feat_distances, 
       treatments = tr, 
       treatment_constraints = c('0' = min.treat.size.per.group, 
                                 '1' = min.treat.size.per.group),
       target = NULL,
       caliper = NULL)
  # table(unit_match)
  
  # ----------------------------------------------------------------------------
  # Run a k-fold CV to estimate the CATE for each unit
  cate_est <- compute_CATE_estimates(feat, yobs, tr, estimator, k, verbose)
  
  # ----------------------------------------------------------------------------
  #  Evaluate how close the average CATE is to the matching estimated ATE
  eval_df <- data.frame(subset = sort(unique(unit_match)),
                        ATE_matching = NA,
                        ATE_estimator = NA)
  for (subset_n in unique(unit_match)){
    idx_subset <- unit_match == subset_n
    col_subset <- eval_df$subset == subset_n
    eval_df$ATE_matching[col_subset] <- mean(yobs[tr == 1 & idx_subset]) - 
      mean(yobs[tr == 0 & idx_subset])
    eval_df$ATE_estimator[col_subset] <- mean(cate_est[idx_subset])
  }
  
  gof_norm <- mean((eval_df$ATE_estimator - eval_df$ATE_matching)^2)
  gof_sd <- sd((eval_df$ATE_estimator - eval_df$ATE_matching)^2)  / sqrt(nrow(eval_df))
  # ----------------------------------------------------------------------------
  return(c(gof_norm, gof_sd))
}
