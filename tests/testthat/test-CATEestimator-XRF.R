library(testthat)
test_that("Tests that XhRF is working correctly", {
  context('X-RF')
  set.seed(1423614230)

  feat <- iris[, -1]
  tr <- rbinom(nrow(iris), 1, .5)
  yobs <- iris[, 1]

  expect_output(xl <- X_RF(feat = feat,
                           tr = tr,
                           yobs = yobs,
                           verbose = TRUE))
  EstimateCate(xl, feat)[1]
  expect_equal(EstimateCate(xl, feat)[1], 0.07867132, tolerance = 1)

  set.seed(432)
  cate_problem <-
    simulate_causal_experiment(
      ntrain = 400,
      ntest = 100,
      dim = 20,
      alpha = .1,
      feat_distribution = "normal",
      testseed = 543,
      trainseed = 234)

  expect_output(xl <- X_RF(feat = cate_problem$feat_tr,
                           yobs = cate_problem$Yobs_tr,
                           tr = cate_problem$W_tr, 
                           verbose = TRUE))
  
  expect_equal(mean((
    EstimateCate(xl, cate_problem$feat_te) - cate_problem$tau_te
  ) ^ 2),
  10.20051,
  tolerance = 1)

})

