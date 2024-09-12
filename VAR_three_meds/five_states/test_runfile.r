library(RcppArmadillo)
library(RcppDist)
library(Rcpp)

Rcpp::sourceCpp("likelihood_fnc_arm.cpp")

sim_dat_num = 4
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

A = list()
W = list()
B = list()

for(i in EIDs){
    A[[i]] = alpha_i_mat[[which(EIDs == i)]]
    W[[i]] = omega_i_mat[[which(EIDs == i)]]
    # B[[i]] = matrix(1, nrow = sum(Y[,"EID"] == i), ncol = 1)
    B[[i]] = matrix(data_format[data_format[,"EID"] == i, "b_true"], ncol = 1)
}

# -----------------------------------------------------------------------------
# Isolate one subject and test efficacy of state sampler and see:
# (a) actual posterior distribution of states
# (b) compare to distribution of state proposals
# (c) compare to "stationary" distribution (marginal probability of states at each time)
# Also, write a sanity check code to ensure that the proposal probabilities are what we expect
# Then, dig into MH ratio and see how we can incorporate the response values
# -----------------------------------------------------------------------------
# i = 502250
# ii = which(EIDs == i)
# EIDs_temp = i
# par_temp = true_pars
# Y_temp = Y[Y[,"EID"] == i, ]
# z_temp = z[Y[,"EID"] == i, ]
# n_i = nrow(Y_temp)
# b_i_true = data_format[data_format[,"EID"] == i,"b_true"]



# -----------------------------------------------------------------------------
# Focusing on a subject -------------------------------------------------------
# -----------------------------------------------------------------------------
# Data set 3
# i = 166350
# i = 259825
# i = 425400
# Data set 4
# i = 288775
# i = 259825
# i = 194350
# i = 234375
# i = 747775

# -----------------------------------------------------------------------------
# Focusing on a few subjects --------------------------------------------------
# -----------------------------------------------------------------------------
id_list = c(194350, 234375, 259825, 288775, 747775)
EIDs_temp = id_list
par_temp = true_pars
ii_s = which(EIDs %in% EIDs_temp)
Y_temp_big = Y[Y[,"EID"] %in% EIDs_temp, ]
z_temp_big = z[Y[,"EID"] %in% EIDs_temp, ]
b_i_true = data_format[data_format[,"EID"] %in% EIDs_temp, "b_true"]
# Initialize a state sequence -------------------------------------------------
set.seed(1)
for(i in EIDs_temp) {
    Y_temp = Y_temp_big[Y_temp_big[,"EID"] == i, ]
    z_temp = z_temp_big[Y_temp_big[,"EID"] == i, ]
    
    n_i = nrow(Y_temp)
    dat_temp = data_format[data_format[,"EID"] == i,]
    
    init_logit = true_pars[par_index$vec_init]
    init_logit = c(0, init_logit)
    P_i = exp(init_logit) / sum(exp(init_logit))
    zeta = matrix(true_pars[par_index$vec_zeta], nrow = 2)
    colnames(zeta) = c('(1) 1->2', '(2) 1->4','(3) 2->3', '(4) 2->4', '(5) 3->1', 
                       '(6) 3->2', '(7) 3->4','(8) 4->2', '(9) 4->5', '(10) 5->1', 
                       '(11) 5->2', '(12) 5->4')
    b_i = NULL
    stat_dist = matrix(ncol = 5, nrow = n_i)
    for(k in 1:n_i){
        if(k==1){
            b_i = sample(1:5, size=1, prob=P_i)
            stat_dist[k,] = P_i
        } else{
            q1   = exp(z_temp[k,, drop=F] %*% zeta[,  1, drop=F]) 
            q2   = exp(z_temp[k,, drop=F] %*% zeta[,  2, drop=F])
            q3   = exp(z_temp[k,, drop=F] %*% zeta[,  3, drop=F])
            q4   = exp(z_temp[k,, drop=F] %*% zeta[,  4, drop=F])
            q5   = exp(z_temp[k,, drop=F] %*% zeta[,  5, drop=F]) 
            q6   = exp(z_temp[k,, drop=F] %*% zeta[,  6, drop=F])
            q7   = exp(z_temp[k,, drop=F] %*% zeta[,  7, drop=F])
            q8   = exp(z_temp[k,, drop=F] %*% zeta[,  8, drop=F])
            q9   = exp(z_temp[k,, drop=F] %*% zeta[,  9, drop=F]) 
            q10  = exp(z_temp[k,, drop=F] %*% zeta[,  10, drop=F])
            q11  = exp(z_temp[k,, drop=F] %*% zeta[,  11, drop=F])
            q12  = exp(z_temp[k,, drop=F] %*% zeta[,  12, drop=F])
            
            Q = matrix(c(   1,   q1,  0,  q2,  0,
                            0,    1, q3,  q4,  0,
                            q5,   q6,  1,  q7,  0,
                            0,   q8,  0,   1, q9,
                            q10,  q11,  0, q12,  1), ncol=5, byrow=T)
            
            P_i = Q / rowSums(Q)
            # Sample the latent state sequence
            b_i = c( b_i, sample(1:5, size=1, prob=P_i[tail(b_i,1),]))
            
            stat_dist[k,] = stat_dist[k-1,,drop=F] %*% P_i
        }
        
    }
    
    B[[which(EIDs == i)]] = matrix(b_i, ncol = 1)
    # B[[ii]] = matrix(1, nrow = length(b_i), ncol = 1)
}


Dn_Xn = update_Dn_Xn_cpp( as.numeric(EIDs), B, Y, true_pars, par_index, x, 10)
Dn = Dn_Xn[[1]]; names(Dn) = EIDs
Xn = Dn_Xn[[2]]

A_temp = B_temp = Dn_temp = Xn_temp = Dn_omega_temp = W_temp = list(); 

for(i in 1:length(EIDs_temp)) {
    A_temp[[i]] = alpha_i_mat[[which(EIDs == EIDs_temp[i])]]
    B_temp[[i]] = B[[which(EIDs == EIDs_temp[i])]]
    Dn_temp[[i]] = Dn[[which(EIDs == EIDs_temp[i])]]
    Xn_temp[[i]] = Xn[[which(EIDs == EIDs_temp[i])]]
    Dn_omega_temp[[i]] = Dn_omega[[which(EIDs == EIDs_temp[i])]]
    W_temp[[i]] = W[[which(EIDs == EIDs_temp[i])]]
}

bleed_indicator_temp = bleed_indicator[Y[,"EID"] %in% EIDs_temp]
n_cores = 10
#  ----------------------------------------------------------------------------

# Sampling type ---------------------------------------------------------------
# samp_type = 1: random selection MH
# samp_type = 2: transition prob. MH
# samp_type = 3: sub. likelihood  MH
# samp_type = 4: Gibbs update

testing_combo_list = cbind(rep(1:4, 4), rep(2:5, each=4))

# Rcpp::sourceCpp("likelihood_fnc_arm.cpp")
# args = commandArgs(TRUE)
# samp_and_p = as.numeric(args[1]) # a number between 1 and 16
samp_and_p = 1

samp_type = testing_combo_list[samp_and_p, 1]
t_pt_length = testing_combo_list[samp_and_p, 2]

it_length = 2000
post_prob_b = matrix(nrow = it_length, ncol = nrow(Y_temp_big))
accur_b = matrix(0, nrow = it_length, ncol = length(EIDs_temp) + 1)
colnames(accur_b) = c('t_diff', paste0('percent_', EIDs_temp))

start_t = Sys.time()
algorithm_s_time = proc.time()
for(it in 1:it_length) {
    if(it %% 10 == 0) {print(it)}
    
    if(samp_type == 1) {
        B_Dn = update_b_i_MH(EIDs_temp, par_temp, par_index, A_temp, B_temp, Y_temp_big, 
                             z_temp_big, Dn_temp, Xn_temp, Dn_omega_temp, W_temp, 
                             bleed_indicator_temp, n_cores, t_pt_length, 1)
    } else if(samp_type == 2) {
        B_Dn = update_b_i_MH(EIDs_temp, par_temp, par_index, A_temp, B_temp, Y_temp_big, 
                             z_temp_big, Dn_temp, Xn_temp, Dn_omega_temp, W_temp, 
                             bleed_indicator_temp, n_cores, t_pt_length, 2)
    } else if(samp_type == 3) {
        B_Dn = update_b_i_MH(EIDs_temp, par_temp, par_index, A_temp, B_temp, Y_temp_big, 
                             z_temp_big, Dn_temp, Xn_temp, Dn_omega_temp, W_temp, 
                             bleed_indicator_temp, n_cores, t_pt_length, 3)
    } else {
        B_Dn = update_b_i_gibbs(EIDs_temp, par_temp, par_index, A_temp, B_temp, Y_temp_big, 
                                z_temp_big, Dn_temp, Xn_temp, Dn_omega_temp, W_temp, 
                                bleed_indicator_temp, n_cores, t_pt_length)
    }
    
    B_temp = B_Dn[[1]]
    Dn_temp = B_Dn[[2]]
    
    # Posterior probability distribution
    post_prob_b[it, ] = do.call( 'c', B_temp)
    
    # Accuracy and timeliness of the sampling routine
    algorithm_e_time = proc.time();
    accur_b[it, 1] = algorithm_e_time[3] - algorithm_s_time[3]
    
    for(sub in 2:ncol(accur_b)) {
        b_ind = which(Y_temp_big[,"EID"] == EIDs_temp[sub-1])
        accur_b[it, sub] = mean(post_prob_b[it, b_ind] == b_i_true[b_ind])
    }
    
    if(algorithm_e_time[3] - algorithm_s_time[3] > 36000) {
        # Max time limit is 10 hours
        print(paste0("Reached maximum time with ", it, " iterations"))
        
        post_prob_b = post_prob_b[1:it, ]
        accur_b = accur_b[1:it,]
        break;
    }
}
end_t = Sys.time()
print(paste0("Total MCMC time: ", end_t - start_t))

sampling_out = list(post_prob_b = post_prob_b, accur_b = accur_b)
save(sampling_out, 
     file = paste0('Model_out/sampling_out_',samp_type,'_',t_pt_length,'.rda'))


# Plot the results
pdf(paste0('Plots/sampling_out_',samp_type,'_',t_pt_length,'.pdf'))
par(mfrow = c(2,1))
for(i in 1:length(EIDs_temp)) {
    b_ind = which(Y_temp_big[,"EID"] == EIDs_temp[i])
    pb = barplot(rbind(colMeans(post_prob_b[1:it_length, b_ind] == 1),
                       colMeans(post_prob_b[1:it_length, b_ind] == 2),
                       colMeans(post_prob_b[1:it_length, b_ind] == 3),
                       colMeans(post_prob_b[1:it_length, b_ind] == 4),
                       colMeans(post_prob_b[1:it_length, b_ind] == 5)), 
                 col=c( 'dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray'), 
                 xlab='time', space=0, col.main='black', border=NA, axes = F, plot = F)
    
    barplot(rbind(colMeans(post_prob_b[1:it_length, b_ind] == 1),
                  colMeans(post_prob_b[1:it_length, b_ind] == 2),
                  colMeans(post_prob_b[1:it_length, b_ind] == 3),
                  colMeans(post_prob_b[1:it_length, b_ind] == 4),
                  colMeans(post_prob_b[1:it_length, b_ind] == 5)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray'), 
            xlab='true states', space=0, col.main='black', border=NA,
            xlim=range(pb) + c(-0.5,0.5), main = EIDs_temp[i]) 
    grid( nx=20, NULL, col='white')
    axis( side=1, at=1:length(b_ind)-0.5, labels = b_i_true[b_ind], cex.axis=0.3)
}

plot(x = accur_b[,1], y = accur_b[,2], type = 'l', ylim=c(0,1), lwd=2,
     xlab = 'compute time (sec.)', ylab = 'state accuracy (at each it.)')
for(j in 2:length(EIDs_temp)) {
    lines(x = accur_b[,1], y = accur_b[,j+1], lwd = 2, col = j)
}
dev.off()

# post_prob_b_no = matrix(nrow = 1000, ncol = n_i)
# for(it in 1:1000) {
#     print(it)
#     prev_prob = 1
#     for(k in 1:(n_i - t_pt_length + 1)) {
#         if (k == 1) {
#             # () -> () -> 1-5
#             omega_set = omega_list_4[[1]][[b_i[t_pt_length+1]]]
#             omega_set_no = omega_list_4[[1]][[b_i_no[t_pt_length+1]]]
#         } else if (k <= n_i - t_pt_length) {
#             # 1-5 -> () -> () -> 1-5
#             omega_set = omega_list_4[[2]][b_i[k - 1], b_i[k + t_pt_length]][[1]]
#             omega_set_no = omega_list_4[[2]][b_i_no[k - 1], b_i_no[k + t_pt_length]][[1]]
#         } else if (k == n_i - t_pt_length + 1) {
#             # 1-5 -> () -> ()
#             omega_set = omega_list_4[[3]][[b_i[n_i - t_pt_length]]]
#             omega_set_no = omega_list_4[[3]][[b_i_no[n_i - t_pt_length]]]
#         }
#         
#         # Informed sampling
#         ss_prob = state_prob_gibbs(k, n_i, t_pt_length, true_pars, par_index, z_i, b_i, omega_set, prev_prob)
#         # row_ind = sample(x = 1:nrow(omega_set), size = 1, prob = c(ss_prob[,1]))
#         row_ind = csample_num(x = 1:nrow(omega_set), size = 1, prob = c(ss_prob[,1]), replace = F)
#         prev_prob = ss_prob[row_ind, 2]
#         b_i[k:(k+t_pt_length-1)] = c(omega_set[row_ind,])
#         
#         # Uninformed sampling
#         row_ind_no = sample(x = 1:nrow(omega_set_no), size = 1)
#         b_i_no[k:(k+t_pt_length-1)] = c(omega_set_no[row_ind_no,])
#     }
#     
#     post_prob_b[it, ] = b_i
#     post_prob_b_no[it, ] = b_i_no
# }


# barplot(rbind(colMeans(post_prob_b_no[, 1:n_i] == 1),
#               colMeans(post_prob_b_no[, 1:n_i] == 2),
#               colMeans(post_prob_b_no[, 1:n_i] == 3),
#               colMeans(post_prob_b_no[, 1:n_i] == 4),
#               colMeans(post_prob_b_no[, 1:n_i] == 5)), 
#         col=c( 'dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray'), 
#         xlab='time', space=0, col.main='green', border=NA,
#         xlim=range(pb) + c(-0.5,0.5)) 
# grid( nx=20, NULL, col='white')
# axis( side=1, at=1:n_i-0.5, col.axis='green', labels = 1:n_i)


# -----------------------------------------------------------------------------
# Focusing on subject 72450 ---------------------------------------------------
# -----------------------------------------------------------------------------
# Rcpp::sourceCpp("likelihood_fnc_arm.cpp")
# 
# i = 72450
# ii = which(EIDs == i)
# t_pts = c(28:54)
# A = alpha_i_mat[[ii]]
# Y_temp = Y[Y[,"EID"] == i, ]
# z_temp = z[Y[,"EID"] == i, ]
# Dn_omega_temp = Dn_omega[[ii]]
# W_temp = W[[ii]]
# 
# sub_interest = data_format[data_format[,"EID"] == i, ]
# 
# 
# b_tester = function(b) {
#     
#     B[[ii]] = matrix(b, ncol = 1)
#     
#     Dn_Xn = update_Dn_Xn_cpp( as.numeric(EIDs), B, Y, par, par_index, x, 10)
#     Dn = Dn_Xn[[1]]; names(Dn) = EIDs
#     Xn = Dn_Xn[[2]]
#     
#     Dn_temp = Dn[[ii]]
#     Xn_temp = Xn[[ii]]
#     
#     like_val = log_f_i_cpp(i, ii, t_pts, par, par_index, A, B[[ii]], Y_temp, z_temp, Dn_temp,
#                            Xn_temp, Dn_omega_temp, W_temp)
#     
#     return(like_val)
# }
# 
# # Initial state sequence
# b_tester(rep(1, 54))
# # True state sequence
# b_tester(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
#            2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3))
# # Correct identification of change
# b_tester(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
#            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2))


# # -----------------------------------------------------------------------------
# # Function that initializes the state space based on the "highest" likelihood -
# # -----------------------------------------------------------------------------
# Xn_initial = update_Dn_Xn_cpp( as.numeric(EIDs), B, Y, par, par_index, x, 10)
# Xn = Xn_initial[[2]]
# 
# B = list()
# B = initialize_b_i(as.numeric(EIDs), par, par_index, A, Y, z, Xn, Dn_omega, W, 10)

# # How close is this to the true state sequences
# true_B = data_format[,"b_true"]
# initial_B = c(do.call( rbind, B))
# side_by_side_state = cbind(data_format[,"EID"], cbind(true_B, initial_B))
# not_correct = side_by_side_state[(side_by_side_state[,2] != side_by_side_state[,3]), ]
# 
# print(paste0("proportion incorrect: ", nrow(not_correct) / nrow(side_by_side_state)))
# # -----------------------------------------------------------------------------
# # -----------------------------------------------------------------------------