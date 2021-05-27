library(testthat)
test_that("Tests that ShRF is working correctly", {
  context('S-RF')
  set.seed(1423614230)

  feat <- iris[, -1]
  tr <- rbinom(nrow(iris), 1, .5)
  yobs <- iris[, 1]

  sl <- S_RF(
    feat = feat,
    tr = tr,
    yobs = yobs)
  
  context("Test conformal intervals")
  CI <- CateCI(sl, 
               feature_new = feat[1:75,],
               bootstrapVersion = "conformal")
  
  expect_equal(nrow(CI),75)
  
  
  context("Test double bootstrap intervals")
  # CI <- CateCI(sl, 
  #              feature_new = feat[1:75,],
  #              bootstrapVersion = "doubleBootstrap",
  #              B = 50,
  #              B_Second = 10)
  
  context("Test disallowed settings")

  # expect_equal(EstimateCate(sl, feat)[1], 0.0369185, tolerance = 1e-4)

  set.seed(432)
  cate_problem <-
    simulate_causal_experiment(
      ntrain = 400,
      ntest = 100,
      dim = 20,
      alpha = .1,
      feat_distribution = "normal",
      testseed = 543,
      trainseed = 234
    )

  sl <- S_RF(
    feat = cate_problem$feat_tr,
    yobs = cate_problem$Yobs_tr,
    tr = cate_problem$W_tr,
    mu.forestry =
      list(
        relevant.Variable = 1:ncol(cate_problem$feat_tr),
        ntree = 20,
        replace = TRUE,
        sample.fraction = 1,
        mtry = ncol(feat),
        nodesizeSpl = 1,
        nodesizeAvg = 3,
        nodesizeStrictSpl = 1,
        nodesizeStrictAvg = 3,
        splitratio = 1,
        middleSplit = FALSE,
        OOBhonest = TRUE))

  context("Allow passing aggregation options to predict")
  preds_oob <- EstimateCate(sl, feature_new = sl@feature_train,aggregation = "oob")
  
  expect_equal(length(preds_oob), 400)
  
  preds <- EstimateCate(sl, feature_new = sl@feature_train,aggregation = "average")
  
  expect_equal(length(preds), 400)
  
  expect_equal(mean((
    EstimateCate(sl, cate_problem$feat_te) - cate_problem$tau_te
  ) ^ 2),
  22.5,
  tolerance = 1)

})
