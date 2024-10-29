source('mcmc_routine_arma.r')

args = commandArgs(TRUE)
sampling_num = as.numeric(args[1])
seed_num = NULL
if(sampling_num <= 4) {
    seed_num = 1
} else if(sampling_num > 4 & sampling_num <= 8) {
    seed_num = 2
    sampling_num = sampling_num - 4
} else {
    seed_num = 3
    sampling_num = sampling_num - 8
} 

# df_num_list = rep(1:10, each = 3)
# df_num = df_num_list[seed_num]
df_num = 4

set.seed(seed_num)
ind = seed_num
# ind_list = rep(1:3, 10)
# ind = ind_list[seed_num]
# print(ind)

simulation = T
data_format = NULL

if(simulation) {
    steps  = 20000
    burnin =  5000
    sim_dat_num = 5
    
    load(paste0('Data_sim/use_data1_', sim_dat_num, '.rda'))
    data_format = use_data
    trialNum = 3
    max_ind = 5

    print(paste0('seed ', seed_num, ' samp ', sampling_num, ' trial ', trialNum))
} else {
    steps  = 50000
    burnin = 5000
    
    data_name = paste0('Data_updates/data_format_', df_num, '.rda')
    load(data_name)
    
    trialNum = 1
    max_ind = 5
    if(max_ind > 5) burnin = 0
}

Y = data_format[, c('EID','hemo', 'hr', 'map', 'lactate', 
                    'RBC_rule', 'clinic_rule')] 
EIDs = as.character(unique(data_format[,'EID']))

x = data_format[,c('n_RBC_admin'), drop=F]
p = ncol(x)

z = cbind(1, data_format[,c('RBC_ordered'), drop=F])
m = ncol(z)

# Parameters ------------------------------------------------------------------
beta = c(0.25, -1, 2, -0.25)

# columns: hemo, hr, map, lactate
alpha_tilde = c( 9.57729783,          -1,        0.1, 0, 0,
                88.69780576,  5.04150472,         -4, 0, 0,
                79.74903940, -5.04150472,          4, 0, 0,
                  5.2113319,   0.5360813, -0.6866748, 0, 0)

sigma_upsilon = c(diag(c(  4, 2, 2, 25, 25, 
                         100, 9, 9, 25, 25, 
                         100, 9, 9, 25, 25, 
                           4, 2, 2, 25, 25)))

vec_A = rep(1.5, 20)

# columns: hemo, hr, map, lactate
R = diag(c(4.58, 98.2, 101.3, 7.6))

# transitions: 1->2, 1->4, 2->3, 2->4, 3->1, 3->2, 3->4, 4->2, 4->5, 5->1, 5->2, 5->4
# zeta = c(-7.2405, 1.5, -5.2152,   1, -2.6473,  -1, -5.1475,  -1, 
#          -9.4459,  -1, -7.2404,   2, -5.2151,   1, -7.1778, 1.5, 
#          -2.6523,   0, -9.4459,  -1, -7.2404, 1.5, -5.2151,   1)
zeta = c(-6.2405, 3.5, -5.2152,   1, -3.6473,  -2, -5.1475,  -2, 
         -7.4459,  -1, -6.2404,   2, -5.2151,   1, -5.1778, 2.5, 
         -2.6523,   0, -6.4459,  -1, -7.2404, 3.5, -4.2151,   1)

omega = c(-1, -1,  1, -1,  1,  1, -1, -1, -1,  1,  1,  1,  1,
          -1,  1,  1,  1, -1,  1, -1,  1, -1, -1, -1, -1, -1,
          -1, -1, -1, -1,  1, -1,  1, -1, -1,  1,  1, -1,  1,
          -1, -1,  1,  1, -1, -1, -1, -1, -1, -1,  1, -1,  1,
           1, -1,  1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,
          -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1, -1,
          -1, -1,  1, -1, -1, -1, -1, -1, -1,  1)

omega = 3 * omega
upsilon_omega = rep(1, length(omega))

init_logit = c(-1, -0.5, -0.1, -0.3)
init_logit = exp(init_logit)

par = c(beta, c(alpha_tilde), c(sigma_upsilon), c(vec_A), c(R), c(zeta), 
        log(init_logit), omega, log(upsilon_omega))
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
# -----------------------------------------------------------------------------

if(simulation) {
    load(paste0('Data_sim/true_pars_', sim_dat_num, '.rda'))
    load(paste0('Data_sim/alpha_i_mat_', sim_dat_num, '.rda'))
    load(paste0('Data_sim/omega_i_mat_', sim_dat_num, '.rda'))
    load(paste0('Data_sim/Dn_omega_sim_', sim_dat_num, '.rda'))
    load(paste0('Data_sim/bleed_indicator_sim_', sim_dat_num,'.rda'))

    par = true_pars
    Dn_omega = Dn_omega_sim
    
    b_chain = data_format[, "b_true"]
} else {
    # ----------------------------------------------------------------------
    if(max_ind > 5) {
        prev_file = paste0('Model_out/mcmc_out_interm_', ind, '_', trialNum, 
                           'it', max_ind-5, '_df', df_num, '.rda')
        load(prev_file)
        
        par_temp = mcmc_out_temp$chain[nrow(mcmc_out_temp$chain), ]
        rownames(par_temp) = NULL
        par = par_temp
        
        b_chain = c(mcmc_out_temp$B_chain[nrow(mcmc_out_temp$B_chain), ])
        
        rm(mcmc_out_temp)
    }
    # ----------------------------------------------------------------------
    
    load('Data_updates/Dn_omega1.rda')
    load('Data_updates/all_EIDs.rda')
    Dn_omega_big = Dn_omega
    eid_index = which(all_EIDs %in% EIDs)
    Dn_omega = vector(mode = 'list', length = length(eid_index))
    for(jjj in 1:length(eid_index)) {
        Dn_omega[[jjj]] = Dn_omega_big[[eid_index[jjj]]]
    }
    rm(Dn_omega_big)
    
    bleed_indicator = b_ind_fnc(data_format)
}
# -----------------------------------------------------------------------------
A = list()
W = list()
B = list()

for(i in EIDs){
    if(simulation) {
        A[[i]] = alpha_i_mat[[which(EIDs == i)]]
        W[[i]] = omega_i_mat[[which(EIDs == i)]]
    } else {
        A[[i]] = matrix(par[par_index$vec_alpha_tilde], ncol =1)
        W[[i]] = matrix(par[par_index$omega_tilde], ncol =1)
    }
    
    if(max_ind > 5) {
        b_temp = b_chain[data_format[,"EID"] == as.numeric(i)]   
    } else {
        init_logit = par[par_index$vec_init]
        init_logit = c(0, init_logit)
        P_i = exp(init_logit) / sum(exp(init_logit))
        zeta = matrix(par[par_index$vec_zeta], nrow = 2)
        colnames(zeta) = c('(1) 1->2', '(2) 1->4','(3) 2->3', '(4) 2->4', '(5) 3->1', 
                           '(6) 3->2', '(7) 3->4','(8) 4->2', '(9) 4->5', '(10) 5->1', 
                           '(11) 5->2', '(12) 5->4')
        z_temp = z[Y[,"EID"] == i, ]
        n_i = nrow(z_temp)
        
        b_temp = NULL
        for(k in 1:n_i){
            if(k==1){
                b_temp = sample(1:5, size=1, prob=P_i)
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
                b_temp = c( b_temp, sample(1:5, size=1, prob=P_i[tail(b_temp,1),]))
            }
            
        }
    }
    
    B[[i]] = matrix(b_temp, ncol = 1)
}
# save(B, file = 'Data_sim/init_B.rda')
# load('Data_sim/init_B.rda')
# -----------------------------------------------------------------------------

print("Starting values for the chain")
print("alpha_tilde")
print(round(par[par_index$vec_alpha_tilde], 3))

print("A")
vec_A_t_logit = par[par_index$vec_A]
vec_A_t = exp(vec_A_t_logit) / (1 + exp(vec_A_t_logit))
mat_A_t = matrix(vec_A_t, nrow = 4)
print(mat_A_t)

print("R")
R_t = matrix(par[par_index$vec_R], ncol = 4)
print(R_t)

print("zeta")
zed = matrix(par[par_index$vec_zeta], nrow = 2)
colnames(zed) = c('(1) 1->2', '(2) 1->4','(3) 2->3', '(4) 2->4', '(5) 3->1', 
                   '(6) 3->2', '(7) 3->4','(8) 4->2', '(9) 4->5', '(10) 5->1', 
                   '(11) 5->2', '(12) 5->4')
print(zed)

vec_A1 = par[par_index$vec_A]
scale_A1 = exp(vec_A1) / (1+exp(vec_A1))

diag_gamma = c(R_t[1,1] / (scale_A1[1]^2), R_t[2,2] / (scale_A1[2]^2),
               R_t[3,3] / (scale_A1[3]^2), R_t[4,4] / (scale_A1[4]^2),
               R_t[1,1] / (scale_A1[5]^2), R_t[2,2] / (scale_A1[6]^2),
               R_t[3,3] / (scale_A1[7]^2), R_t[4,4] / (scale_A1[8]^2),
               R_t[1,1] / (scale_A1[9]^2), R_t[2,2] / (scale_A1[10]^2),
               R_t[3,3] / (scale_A1[11]^2), R_t[4,4] / (scale_A1[12]^2),
               R_t[1,1] / (scale_A1[13]^2), R_t[2,2] / (scale_A1[14]^2),
               R_t[3,3] / (scale_A1[15]^2), R_t[4,4] / (scale_A1[16]^2),
               R_t[1,1] / (scale_A1[17]^2), R_t[2,2] / (scale_A1[18]^2),
               R_t[3,3] / (scale_A1[19]^2), R_t[4,4] / (scale_A1[20]^2))
print(round(sqrt(diag_gamma), 3))

s_time = Sys.time()
mcmc_out = mcmc_routine( par, par_index, A, W, B, Y, x, z, steps, burnin, ind, 
                         trialNum, Dn_omega, simulation, bleed_indicator, 
                         max_ind, df_num, sampling_num)
e_time = Sys.time() - s_time; print(e_time)


# Summary statistics for the transition rate parameter estimates
# transitions: 1->2, 1->4, 2->3, 2->4, 3->1, 3->2, 3->4, 4->2, 4->5, 5->1, 5->2, 5->4
# load('Data_updates/data_format.rda')
# EID_big = unique(data_format[,"EID"])
# est_bleeds = unique(data_format[data_format[,"RBC_rule"] != 0 | data_format[,"clinic_rule"] == 1, "EID"])
# 
# los_patients = matrix(0, ncol = 2, nrow = length(EID_big))
# colnames(los_patients) = c('EID', 'los')
# for(i in 1:length(EID_big)){
#     los_patients[i,1] = EID_big[i]
#     los_patients[i,2] = sum(data_format[,"EID"] == EID_big[i])
# }
# # Mean and median length of stay (los)
# mean(los_patients[,2])
# median(los_patients[,2])

