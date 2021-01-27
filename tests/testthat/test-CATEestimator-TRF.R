library(testthat)
test_that("Tests that TRF is working correctly", {
  context('T-RF')
  set.seed(1423614230)

  feat <- iris[, -1]
  tr <- rbinom(nrow(iris), 1, .5)
  yobs <- iris[, 1]

  tl <- T_RF(feat = feat, tr = tr, yobs = yobs)
  EstimateCate(tl, feat)[1]
  # expect_equal(EstimateCate(tl, feat)[1], 0.1143112, tolerance = 1e-2)

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

  tl <- T_RF(
    feat = cate_problem$feat_tr,
    yobs = cate_problem$Yobs_tr,
    tr = cate_problem$W_tr
  )
  
  # expect_equal(mean((
  #   EstimateCate(tl, cate_problem$feat_te) - cate_problem$tau_te
  # ) ^ 2),
  # 70.08503,
  # tolerance = 1)

})
