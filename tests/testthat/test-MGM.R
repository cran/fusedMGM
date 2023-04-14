test_that("Making MGM works", {
  library(fusedMGM)
  options(bigmemory.allow.dimnames=TRUE)
  
  data("data_all")
  data("ind_disc")
  
  X <- data_all[,-ind_disc, drop=FALSE] ;
  Y <- data_all[,ind_disc, drop=FALSE] ;
  g <- rep(c(1,2), each=250) ;
  
  expect_s3_class(MGM(X,Y,g), "MGM") ;
})