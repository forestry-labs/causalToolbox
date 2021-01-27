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
    yobs = yobs, 
    mu.forestry =
      list(
        relevant.Variable = 1:ncol(feat),
        ntree = 20,
        replace = TRUE,
        sample.fraction = 0.9,
        mtry = ncol(feat),
        nodesizeSpl = 1,
        nodesizeAvg = 3,
        splitratio = .5,
        middleSplit = FALSE
      ))

  # expect_equal(EstimateCate(sl, feat)[1], 0.0491662, tolerance = 1e-4)

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
        relevant.Variable = 1:ncol(feat),
        ntree = 20,
        replace = TRUE,
        sample.fraction = 0.9,
        mtry = ncol(feat),
        nodesizeSpl = 1,
        nodesizeAvg = 3,
        splitratio = .5,
        middleSplit = FALSE))

  # expect_equal(mean((
  #   EstimateCate(sl, cate_problem$feat_te) - cate_problem$tau_te
  # ) ^ 2),
  # 31.19558,
  # tolerance = 1)

})
