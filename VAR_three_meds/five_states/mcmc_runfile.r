source('mcmc_routine_arma.r')

args = commandArgs(TRUE)
ind = as.numeric(args[1])

set.seed(ind)
print(ind)

simulation = F
data_format = NULL

if(simulation) {
    steps  = 50000
    burnin =  5000
    sim_dat_num = 1
    
    load(paste0('Data_sim/use_data1_', sim_dat_num, '.rda'))
    data_format = use_data
    trialNum = 1
    
    max_ind = 5
} else {
    steps  = 50000
    burnin = 5000
    real_dat_num = 3
    
    # Loading the old data to keep the same IDs -----------------
    load(paste0('../Data/data_format_new', real_dat_num, '.rda'))
    old_EIDs = unique(data_format[,"EID"])
    old_time = data_format[,"time"]
    old_ids = data_format[,"EID"]
    
    load('Data_updates/data_format.rda')
    new_EIDs = unique(data_format[,'EID'])

    data_format = data_format[data_format[,"EID"] %in% old_EIDs, ]
    # ID 642
    
    # trial 2: starting seed was 3
    # trial 3: starting seed was 1
    # trial 4: continuation of trial 2
    # trial 5: new EIDs
    # trial 6: new EIDs and new initial values
    # trial 7 (it 1 - 5): new EIDs, trial 6 initial values, and R update has greater variance
    # trial 8: new EIDs, old initial values and old priors
    # trial 9: new EIDs, continuation of trial 7, more diffuse R prior, tighter prior on beta
    # trial 10: strong prior on R
    # trial 11: diffuse prior on R, new R proposal scheme
    # trial 12: reintroduce hemo,lact > 0, change priors for zeta and alpha tilde
    trialNum = 12
    max_ind = 5
}

Y = data_format[, c('EID','hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
EIDs = as.character(unique(data_format[,'EID']))

x = data_format[,c('n_RBC_admin'), drop=F]
p = ncol(x)

z = cbind(1, data_format[,c('RBC_ordered'), drop=F])
m = ncol(z)

# Parameters ------------------------------------------------------------------
beta = c(0.25, -1, 2, -0.25)

# columns: hemo, hr, map, lactate
alpha_tilde = c(9.57729783,          -1,        0.1, 0, 0,
                88.69780576,  5.04150472,         -4, 0, 0,
                79.74903940, -5.04150472,          4, 0, 0,
                5.2113319,   0.5360813, -0.6866748, 0, 0)

sigma_upsilon = c(diag(c(9, 2, 2, 2, 2, 
                         400, 9, 9, 9, 9, 
                         100, 9, 9, 9, 9, 
                         9, 2, 2, 2, 2)))

vec_A = rep(0, 20)

# columns: hemo, hr, map, lactate
R = diag(c(4.58, 98.2, 101.3, 7.6))

# transitions: 1->2, 1->4, 2->3, 2->4, 3->1, 3->2, 3->4, 4->2, 4->5, 5->1, 5->2, 5->4
zeta = c(-7.2405, 1.5, -5.2152,   1, -2.6473,  -1, -5.1475,  -1, 
         -9.4459,  -1, -7.2404,   2, -5.2151,   1, -7.1778, 1.5, 
         -2.6523,   0, -9.4459,  -1, -7.2404, 1.5, -5.2151,   1)

omega = c(-1, -1,  1, -1,  1,  1, -1, -1, -1,  1,  1,  1,  1,
          -1,  1,  1,  1, -1,  1, -1,  1, -1, -1, -1, -1, -1,
          -1, -1, -1, -1,  1, -1,  1, -1, -1,  1,  1, -1,  1,
          -1, -1,  1,  1, -1, -1, -1, -1, -1, -1,  1, -1,  1,
           1, -1,  1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,
          -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1, -1,
          -1, -1,  1, -1, -1, -1, -1, -1, -1,  1)

omega = 3 * omega
upsilon_omega = rep(1, length(omega))

init_logit = c(-1, -1, 0.5, -1)
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
} else {
    # ----------------------------------------------------------------------
    prev_file = 'Model_out/mcmc_out_interm_3_1it2.rda'
    load(prev_file)
    # ----------------------------------------------------------------------
    
    load('Data_updates/Dn_omega1.rda')
    Dn_omega_big = Dn_omega
    eid_index = which(new_EIDs %in% old_EIDs)
    Dn_omega = vector(mode = 'list', length = length(eid_index))
    for(jjj in 1:length(eid_index)) {
        Dn_omega[[jjj]] = Dn_omega_big[[eid_index[jjj]]]
    }
    rm(Dn_omega_big)
    
    bleed_indicator = b_ind_fnc(data_format)
    
    par_temp = colMeans(mcmc_out_temp$chain[800:nrow(mcmc_out_temp$chain), ])
    rownames(par_temp) = NULL
    par[c(par_index$vec_A, par_index$vec_R)] = par_temp[c(par_index$vec_A, par_index$vec_R)]
    
    b_chain = c(mcmc_out_temp$B_chain[nrow(mcmc_out_temp$B_chain), ])

    print("initial state sequence based on:")
    print(prev_file)
    
    rm(mcmc_out_temp)
}
# -----------------------------------------------------------------------------
A = list()
W = list()
B = list()

for(i in EIDs){
    if(simulation) {
        A[[i]] = alpha_i_mat[[which(EIDs == i)]]
        B[[i]] = data_format[data_format[,'EID']==as.numeric(i), "b_true", drop=F]
        W[[i]] = omega_i_mat[[which(EIDs == i)]]
    } else {
        # # ----------------------------------------------------------------------
        # old_t = old_time[old_ids == as.numeric(i)]
        # curr_t = data_format[data_format[,"EID"] == as.numeric(i), "time"]
        # if(sum(floor(old_t) %in% floor(curr_t)) == length(old_t)) {
        #     b_temp = b_chain[old_ids==as.numeric(i)]
        # } else {
        #     b_index = which(floor(old_t) %in% floor(curr_t))
        # 
        #     b_temp_init = b_chain[old_ids==as.numeric(i)]
        #     b_temp = b_temp_init[b_index]
        # }
        # 
        # # temporary fix to there existing state 2 in the clinic rule = -1
        # if(unique(Y[Y[,'EID']==as.numeric(i), 'clinic_rule']) < 0) {
        #     print(paste0("Clinic rule -1: ", i))
        #     b_temp = rep(1, length(b_temp))
        # }
        # # ----------------------------------------------------------------------
        
        # b_temp = b_chain[data_format[,"EID"] == as.numeric(i)]
        b_temp = rep(1, sum(data_format[,"EID"] == as.numeric(i)))
        
        B[[i]] = matrix(b_temp, ncol = 1)
        A[[i]] = matrix(par[par_index$vec_alpha_tilde], ncol =1)
        W[[i]] = matrix(par[par_index$omega_tilde], ncol =1)
    }
}
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
print(zed)

s_time = Sys.time()
mcmc_out = mcmc_routine( par, par_index, A, W, B, Y, x, z, steps, burnin, ind, 
                         trialNum, Dn_omega, simulation, bleed_indicator, max_ind)
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

