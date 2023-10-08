library(RcppArmadillo)
library(RcppDist)
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")

Rcpp::sourceCpp("likelihood_fnc_arm.cpp")
set.seed(5)


test_fnc()
