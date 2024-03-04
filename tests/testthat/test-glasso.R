test_that("glasso", {

  library(psych)
  data(bfi)

  reg.15 <- em.prec(bfi[,1:25], glassoversion = "glasso", rho=.15)
  reg1 <- em.prec(bfi[,1:25], glassoversion = "glasso", rho=1)
  
  expect_snapshot(reg.15$K)
  expect_snapshot(reg1$K)
  
  # others are slow... omit from CRAN testing
  skip_on_cran()
  rho <- seq(.001,.5,length.out = 25)
  ebic1 <- EMggm(bfi[,1:25], rho = rho, glassoversion = "glasso",
                    rhoselect = "ebic")
     
  kfold1 <- EMggm(bfi[,1:25], rho = rho, glassoversion = "glasso",
                  seed = 1234,
                  rhoselect = "kfold")
  
  ebic2 <- EMggm(bfi[,1:25], rho = rho, glassoversion = "glassoFast",
                  rhoselect = "ebic")
  
  expect_snapshot(ebic1$graph)
  expect_snapshot(kfold1$graph)
  expect_snapshot(ebic2$graph)
  
})
