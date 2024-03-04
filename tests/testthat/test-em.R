test_that("EM cov and prec", {

  library(psych)
  data(bfi)
  cov.res <- em.cov(bfi[,1:25])
  
  prec.res <- em.prec(bfi[,1:25])
  
  expect_snapshot(cov.res$S)
  expect_snapshot(cov.res$mu)
  
  expect_snapshot(prec.res$S)
  expect_snapshot(prec.res$K)
  expect_snapshot(prec.res$mu)
  
  expect_equal(cov.res$S, prec.res$S)
  
})
