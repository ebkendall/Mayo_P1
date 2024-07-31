library(RcppArmadillo)
library(RcppDist)
library(Rcpp)

sim_dat_num = 2
load(paste0('Data_sim/use_data1_', sim_dat_num, '.rda'))
data_format = use_data

Y = data_format[, c('EID','hemo', 'hr', 'map', 'lactate', 
                    'RBC_rule', 'clinic_rule')] 
EIDs = as.character(unique(data_format[,'EID']))

x = data_format[,c('n_RBC_admin'), drop=F]
p = ncol(x)

z = cbind(1, data_format[,c('RBC_ordered'), drop=F])
m = ncol(z)

par_index = list()
par_index$vec_beta = 1:4
par_index$vec_alpha_tilde = 5:24
par_index$vec_sigma_upsilon = 25:424
par_index$vec_A = 425:444
par_index$vec_R = 445:460
par_index$vec_zeta = 461:484
par_index$vec_init = 485:488
par_index$omega_tilde = 489:576
par_index$vec_upsilon_omega = 577:664

load(paste0('Data_sim/true_pars_', sim_dat_num, '.rda'))
load(paste0('Data_sim/alpha_i_mat_', sim_dat_num, '.rda'))
load(paste0('Data_sim/omega_i_mat_', sim_dat_num, '.rda'))
load(paste0('Data_sim/Dn_omega_sim_', sim_dat_num, '.rda'))
load(paste0('Data_sim/bleed_indicator_sim_', sim_dat_num,'.rda'))

par = true_pars
Dn_omega = Dn_omega_sim
rm(Dn_omega_sim)

# Artificially increase the noise of the VAR process
par[par_index$vec_R] = c(diag(c(4.58, 98.2, 101.3, 7.6)))

A = list()
W = list()
B = list()

for(i in EIDs){
    A[[i]] = alpha_i_mat[[which(EIDs == i)]]
    W[[i]] = omega_i_mat[[which(EIDs == i)]]
    B[[i]] = matrix(1, nrow = sum(Y[,"EID"] == i), ncol = 1)
}

Rcpp::sourceCpp("likelihood_fnc_arm.cpp")

# -----------------------------------------------------------------------------
# Function that initializes the state space based on the "highest" likelihood -
# -----------------------------------------------------------------------------
Xn_initial = update_Dn_Xn_cpp( as.numeric(EIDs), B, Y, par, par_index, x, 10)
Xn = Xn_initial[[2]]

B = list()
B = initialize_b_i(as.numeric(EIDs), par, par_index, A, Y, z, Xn, Dn_omega, W, 10)


# How close is this to the true state sequences
true_B = data_format[,"b_true"]
initial_B = c(do.call( rbind, B))
side_by_side_state = cbind(data_format[,"EID"], cbind(true_B, initial_B))
not_correct = side_by_side_state[(side_by_side_state[,2] != side_by_side_state[,3]), ]

# -----------------------------------------------------------------------------
# Focusing on subject 72450 ---------------------------------------------------
# -----------------------------------------------------------------------------

i = 72450
ii = which(EIDs == i)
t_pts = c(28:54)
A = alpha_i_mat[[ii]]
Y_temp = Y[Y[,"EID"] == i, ]
z_temp = z[Y[,"EID"] == i, ]
Dn_omega_temp = Dn_omega[[ii]]
W_temp = W[[ii]]

sub_interest = data_format[data_format[,"EID"] == i, ]


b_tester = function(b) {
    
    B[[ii]] = matrix(b, ncol = 1)
    
    Dn_Xn = update_Dn_Xn_cpp( as.numeric(EIDs), B, Y, par, par_index, x, 10)
    Dn = Dn_Xn[[1]]; names(Dn) = EIDs
    Xn = Dn_Xn[[2]]
    
    Dn_temp = Dn[[ii]]
    Xn_temp = Xn[[ii]]
    
    like_val = log_f_i_cpp(i, ii, t_pts, par, par_index, A, B[[ii]], Y_temp, z_temp, Dn_temp,
                           Xn_temp, Dn_omega_temp, W_temp)
    
    return(like_val)
}

# Initial state sequence
b_tester(rep(1, 54))
# True state sequence
b_tester(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
           2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3))
# Correct identification of change
b_tester(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
           1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2))






