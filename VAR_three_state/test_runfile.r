library(MASS, quietly=T)
library(mvtnorm, quietly=T)
library(LaplacesDemon, quietly=T)
library(Matrix, quietly=T)
library(RcppArmadillo)
library(RcppDist)
library(Matrix, quietly=T)
library(rbenchmark)
library(foreach, quietly=T)
library(doParallel, quietly=T)
library(bayesSurv)

library(Rcpp)
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")

# load('test_post.rda')
# load('test_list.rda')
Rcpp::sourceCpp("likelihood_fnc_arm.cpp")
set.seed(5)



