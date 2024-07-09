library(RcppArmadillo)
library(RcppDist)
library(Rcpp)

Rcpp::sourceCpp("likelihood_fnc_arm.cpp")
set.seed(5)

get_Omega_list_r <- function(adj_mat, s) {
    a = list()
    b = list()
    c = list()
    
    for(i in 1:nrow(adj_mat)) {
        a_temp_list = getAllStateSequences_forward(i-1, adj_mat, s+1)
        if(length(a_temp_list) > 0) {
            a_temp = do.call(rbind, a_temp_list)
            a_temp = a_temp[,-1,drop=F]   
        } else {
            a_temp = matrix(-2, nrow = 1, ncol = s)
        }
        
        c_temp_list = getAllStateSequences_backward(i-1, adj_mat, s+1)
        if(length(c_temp_list) > 0) {
            c_temp = do.call(rbind, c_temp_list)
            c_temp = c_temp[,-(s+1),drop=F]
            # c_temp = c_temp[order(c_temp[,s]), ]
        } else {
            c_temp = matrix(-2, nrow = 1, ncol = s)
        }
        
        a[[i]] = matrix(c(a_temp) + 1, ncol = s)
        c[[i]] = matrix(c(c_temp) + 1, ncol = s)
        
        
        b[[i]] = list()
        for(j in 1:ncol(adj_mat)) {
            b_temp_list = getAllStateSequences_both(i-1, j-1, adj_mat, s+2)
            if(length(b_temp_list) > 0) {
                b_temp = do.call(rbind, b_temp_list)
                b_temp = b_temp[,-c(1, s+2),drop=F]
                b[[i]][[j]] = matrix(c(b_temp) + 1, ncol = s)
            } else {
                b[[i]][[j]] = matrix(-1, nrow = 1, ncol = s)
            }
        }
    }
    
    Omega_List = list()
    Omega_List[[1]] = c
    Omega_List[[2]] = b
    Omega_List[[3]] = a
    
    return(Omega_List)
}

adjacency_matrix <- matrix(c(1, 1, 0, 1, 0,
                             0, 1, 1, 1, 0,
                             1, 1, 1, 1, 0,
                             0, 1, 0, 1, 1,
                             1, 1, 0, 1, 1), nrow = 5, byrow = TRUE)

adjacency_matrix_sub <- matrix(c(1, 0, 0, 1, 0,
                                 0, 1, 0, 0, 0,
                                 1, 0, 1, 1, 0,
                                 0, 0, 0, 1, 1,
                                 1, 0, 0, 1, 1), nrow = 5, byrow = TRUE)

t_pt_length = 2
Omega_List_GLOBAL_sub_multi = get_Omega_list(adj_mat = adjacency_matrix_sub, 
                                             s = t_pt_length)
Omega_List_GLOBAL_multi     = get_Omega_list(adj_mat = adjacency_matrix, 
                                             s = t_pt_length)

N = ncol(adjacency_matrix)

# print("Case (c) Full")
# for(w in 1:N) {
#     print(paste0("() -> () -> ",w))
#     print(paste0(nrow(Omega_List_GLOBAL_multi[[1]][[w]]), " combos"))
#     print(Omega_List_GLOBAL_multi[[1]][[w]][order(Omega_List_GLOBAL_multi[[1]][[w]][,2]), ])
# }
# 
# print("Case (b) Full")
# for(i in 1:N) {
#     for(j in 1:N) {
#         print(paste0(i, "-->", j))
#         print(paste0(nrow(Omega_List_GLOBAL_multi[[2]][[i]][[j]]), " combos"))
#         print(Omega_List_GLOBAL_multi[[2]][[i]][[j]])
#     }
# }
# 
# print("Case (a) Full")
# for(w in 1:N) {
#     print(paste0(w, " -> () -> ()"))
#     print(paste0(nrow(Omega_List_GLOBAL_multi[[3]][[w]]), " combos"))
#     print(Omega_List_GLOBAL_multi[[3]][[w]])
# }

test_fnc()
