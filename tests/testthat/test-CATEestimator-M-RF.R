
library(testthat)
test_that(
  "Tests that M-RF is working correctly",
  {
    dir.create("knownTestValues/", showWarnings = FALSE)
    
    context('M-RF')
    set.seed(1423614230)
    feat <- iris[, -1]
    tr <- rbinom(nrow(iris), 1, .5)
    yobs <- iris[, 1]
    
    morf <- M_RF(
      feat = feat,
      tr = tr,
      yobs = yobs
    )

  # expect_equal(EstimateCate(morf, feat)[1], 0.04129601, tolerance = 1e-2)
  
  # expect_known_value(EstimateCate(morf, feat),
  #                     file = "knownTestValues/testvalues-M_RF1", 
  #                    tolerance = 1e-1)
  
  CI <- CateCI(theObject = morf,
         feature_new = feat, 
         B = 5, 
         verbose = FALSE)
  
  # expect_known_value(as.numeric(CI[1, ]),
  #                    file = "knownTestValues/testvalues-M_RF2", 
  #                    tolerance = 1e-2)
  # expect_equal(as.numeric(CI[1, ]), 
  #              c(0.04244312, -0.07651435, 0.1614006), 
  #              tolerance = 1e-2)
  # 
  set.seed(432)
  cate_problem <-
    simulate_causal_experiment(
      ntrain = 400,
      ntest = 100,
      dim = 20,
      alpha = .1,
      feat_distribution = "normal"
    )

  morf <- M_RF(
    feat = cate_problem$feat_tr,
    yobs = cate_problem$Yobs_tr,
    tr = cate_problem$W_tr)

  # expect_equal(mean((
  #   EstimateCate(morf, cate_problem$feat_te) - cate_problem$tau_te
  # ) ^ 2),
  # 25.12625,
  # tolerance = 1)
  
  # expect_known_value(mean((
  #   EstimateCate(morf, cate_problem$feat_te) - cate_problem$tau_te
  # ) ^ 2),
  # file = "knownTestValues/testvalues-M_RF3", 
  # tolerance = 1)


})

