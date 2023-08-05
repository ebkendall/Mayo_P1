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

load('test_post.rda')

proposal_R <- function(nu_R, psi_R, Y, Dn, Xn, A, par, par_index, EIDs){
    
    eids = Y[,1]
    vec_A_total = par[par_index$vec_logit_A]
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
        bold_Z_vec = vec_Y_i - script_N_full;
        bold_Z = matrix(bold_Z_vec, nrow = 4, ncol = ncol(Y_i))
        bold_Z = bold_Z[,-ncol(bold_Z)]
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

nu_R = 6
psi_R = matrix(c(1.6, -0.8,  0.8, -0.8,
                 -0.8,   16, -0.8,  0.8,
                 0.8, -0.8,   16, -0.8,
                 -0.8,  0.8, -0.8,  1.6), nrow = 4, byrow = T)

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

