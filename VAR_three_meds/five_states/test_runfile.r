library(RcppArmadillo)
library(RcppDist)
library(Rcpp)
# Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
# Sys.setenv("PKG_LIBS" = "-fopenmp")

Rcpp::sourceCpp("test_run.cpp")
set.seed(5)


adjacency_matrix <- matrix(c(1, 1, 1,
                             1, 1, 1,
                             1, 0, 1), nrow = 3, byrow = TRUE)

s = 2

a = list()
b = list()
c = list()

for(i in 1:nrow(adjacency_matrix)) {
    a_temp = do.call(rbind, getAllStateSequences_forward(i-1, adjacency_matrix, s+1))
    a_temp = a_temp[,-1]
    c_temp = do.call(rbind, getAllStateSequences_backward(i-1, adjacency_matrix, s+1))
    c_temp = c_temp[,-(s+1)]
    
    a[[i]] = a_temp
    c[[i]] = c_temp
    
    b[[i]] = list()
    for(j in 1:ncol(adjacency_matrix)) {
        b_temp = do.call(rbind, getAllStateSequences_both(i-1, j-1, adjacency_matrix, s+2))
        b_temp = b_temp[,-c(1, s+2)]
        b[[i]][[j]] = b_temp
    }
}
