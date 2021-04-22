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
  #              c(-67.037264, -48.915089, -17.727470 , 72.952805, 29.012947,
  #                -57.851446, -3.209296, -64.463066, -5.855422, 104.064621))
  
})
