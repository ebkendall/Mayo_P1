library(RcppArmadillo)
library(RcppDist)
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")

# load('test_post.rda')
# load('test_list.rda')
Rcpp::sourceCpp("likelihood_fnc_arm.cpp")
set.seed(5)

simulation = T
data_format = NULL
if(simulation) {
    load('Data/use_data1_1.rda')
    data_format = use_data
    trialNum = 1
} else {
    load('Data/data_format_new.rda')
    pace_id = c(53475, 110750, 125025, 260625, 273425, 296500, 310100, 384925,
                417300, 448075, 538075, 616025, 660075, 665850, 666750, 677225,
                732525, 758025, 763050, 843000)
    data_format = data_format[!(data_format[,'EID'] %in% pace_id), ]
    trialNum = 7 # CHANGE THIS EVERY TIME **********************
}

Y = data_format[, c('EID','hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
EIDs = as.character(unique(data_format[,'EID']))

x = data_format[,c('n_RBC_admin'), drop=F]
p = ncol(x)

z = cbind(1, data_format[,c('RBC_ordered'), drop=F])
m = ncol(z)

otype = !is.na(Y[, c('hemo','hr','map','lactate')])
colnames(otype) = c('hemo','hr','map','lactate')

B = list()
for(i in EIDs){
    
    if(simulation) {
        B[[i]] = data_format[data_format[,'EID']==as.numeric(i), "b_true", drop=F]
    } else {
        b_temp = matrix( 1, sum(Y[,'EID']==as.numeric(i)), 1)
        # b_length = nrow(b_temp)
        # b_temp[(b_length-5):b_length, ] = 1
        B[[i]] = b_temp
    }
}

# load('test_post.rda')
load('A_pre.rda')
load('A_post.rda')
load('prop_R_test.rda')
load('true_alpha_i.rda')
A_true = list()
for(i in 1:length(prop_R_test$EIDs)) {
    A_true[[i]] = t(true_alpha_i[i,,drop=F])
}

proposal_R <- function(nu_R, psi_R, Y, Dn, Xn, A, par, par_index, EIDs){
    
    eids = Y[,1]
    vec_A_total = par[par_index$vec_A]
    A_all_state = matrix(vec_A_total, nrow = 4, ncol = 3)
    vec_A = A_all_state[,1] 
    A_1 = diag(vec_A)
    little_a = c(A_1)
    
    psi_prop_R_interm = matrix(0, nrow = 4, ncol = 4)
    
    for(i in 1:length(EIDs)) {
        
        Y_temp = Y[eids == EIDs[i], ]
        Y_i = Y_temp[,2:5]
        Y_i = t(Y_i)
        vec_Y_i = c(Y_i)
        
        vec_alpha_i = A[[i]]
        vec_beta = par[par_index$vec_beta]
        Xn_i = Xn[[i]]
        Dn_i = Dn[[i]]
        
        
        script_N_full = Dn_i %*% vec_alpha_i + Xn_i %*% vec_beta;
        bold_Z = Y_i[,-ncol(Y_i)]
        I_4 = diag(4)
        
        script_Z = t(bold_Z) %x% I_4
        script_N = script_N_full[-(1:4)]
        script_Y = vec_Y_i[-(1:4)]
        
        vec_M = script_N + script_Z %*% little_a
        M = matrix(vec_M, nrow = 4)
        hold = matrix(script_Y, nrow = 4) - M
        psi_prop_R_interm = psi_prop_R_interm + hold %*% t(hold)
    }
    
    psi_prop_R = psi_prop_R_interm + psi_R
    nu_prop_R = nrow(Y) - length(EIDs) + nu_R
    
    return(list(psi_prop_R = psi_prop_R, nu_prop_R = nu_prop_R))
}

nu_R = prop_R_test$nu_R
psi_R = prop_R_test$psi_R
Y = prop_R_test$Y
Dn = prop_R_test$Dn
Xn = prop_R_test$Xn
par = prop_R_test$par
par_index = prop_R_test$par_index
EIDs = prop_R_test$EIDs

test_R = proposal_R(nu_R, psi_R, Y, Dn, Xn, A_true, par, par_index, EIDs)

print("Mean")
print(test_R$psi_prop_R / (test_R$nu_prop_R - 5))

diff_check = matrix(ncol = 3, nrow = length(EIDs))
colnames(diff_check) = c("mean_sample", "mean_true", "sample_true")
for(i in 1:length(EIDs)) {
    diff_check[i,] = c(sum((A_pre[[i]] - A_post[[i]])^2),
                       sum((A_pre[[i]] - A_true[[i]])^2),
                       sum((A_post[[i]] - A_true[[i]])^2))
}


s_time = Sys.time()
curr_psi_nu = proposal_R(nu_R, psi_R, Y, test_post$Dn, test_post$Xn, test_post$A,
                         test_post$par, test_post$par_index, EIDs)
e_time = Sys.time() - s_time; print(e_time)
print(curr_psi_nu)

s_time = Sys.time()
curr_psi_nu_c = proposal_R_cpp(nu_R, psi_R, Y, test_post$Dn, test_post$Xn, test_post$A,
                         test_post$par, test_post$par_index, as.numeric(EIDs))
e_time = Sys.time() - s_time; print(e_time)
print(curr_psi_nu_c)

