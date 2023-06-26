library(RcppArmadillo)
library(RcppDist)
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")

# load('test_post.rda')
# load('test_list.rda')
Rcpp::sourceCpp("likelihood_fnc_arm.cpp")
set.seed(5)


load('Data/data_format_new.rda')

Y = data_format[, c('EID','hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
EIDs = unique(data_format[,'EID'])

otype = !is.na(Y[, c('hemo','hr','map','lactate')])
colnames(otype) = c('hemo','hr','map','lactate')

test_fnc(EIDs, Y, otype)
