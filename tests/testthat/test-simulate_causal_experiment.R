library(testthat)
test_that("Test causal experiment creator", {
  context('causal experiment creator')
  
  
  ce <- simulate_causal_experiment(
    ntrain = 20,
    ntest = 20,
    dim = 3,
    alpha = .1,
    feat_distribution = "normal",
    given_features = NULL,
    pscore = "rct5",
    mu0 = "sparseLinearStrong",
    tau = "sparseLinearWeak",
    testseed = 4972609,
    trainseed = 1184332
  ) 
  
  # expect_equal(ce$Yobs_tr[1:10], 
  #              c(-24.939932, -10.492052, -16.162362, 9.945661, -6.344967, 
  #                -104.458783, 37.289647, 27.835309, 41.925713, 112.604875))
  
})
