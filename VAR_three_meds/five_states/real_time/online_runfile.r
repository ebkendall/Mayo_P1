source('online_routine.r')

args = commandArgs(TRUE)
seed_num = as.numeric(args[1])
set.seed(seed_num)

df_num = 1

ind_list = rep(1:3, 10)
ind = ind_list[seed_num]
print(ind)

data_format = NULL

steps  = 50000
burnin = 0

data_name = paste0('../Data_updates/data_format_', df_num, '.rda')
load(data_name)
print(data_name)

# Subset to only include a few patients which will update over time ------------
pred_id = c(208100, 208575, 209675, 327375, 498400, 
            509800, 510100, 544200, 544450, 820775)
df_sub = data_format[data_format[,"EID"] %in% pred_id, ]

# Start with an hour of data ---------------------------------------------------
df_sub_start = NULL
for(j in unique(df_sub[,"EID"])) {
    df_j = df_sub[df_sub[,'EID'] == j, ]
    df_sub_start = rbind(df_sub_start, df_j[1:4, ])
}

trialNum = 1
max_ind = 10

# Y = data_format[, c('EID','hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
Y = df_sub_start[, c('EID','hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
EIDs = as.character(unique(Y[,'EID']))

x = df_sub_start[,c('n_RBC_admin'), drop=F]
p = ncol(x)

z = cbind(1, df_sub_start[,c('RBC_ordered'), drop=F])
m = ncol(z)

# Parameters ------------------------------------------------------------------
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
prev_file = paste0('../Model_out/mcmc_out_interm_', 3, '_', 15, 'it', 5, '_df', df_num, '.rda')
load(prev_file)
# ----------------------------------------------------------------------

load('../Data_updates/Dn_omega1.rda')
load('../Data_updates/all_EIDs.rda')
Dn_omega_big = Dn_omega
eid_index = which(all_EIDs %in% EIDs)
Dn_omega = vector(mode = 'list', length = length(eid_index))
for(jjj in 1:length(eid_index)) {
    Dn_omega[[jjj]] = Dn_omega_big[[eid_index[jjj]]]
}
rm(Dn_omega_big)

par_temp = mcmc_out_temp$chain[nrow(mcmc_out_temp$chain), ]
rownames(par_temp) = NULL
par = par_temp
par[par_index$vec_zeta] = c(-7.2405, 1.5, -5.2152,   1, -2.6473,  -1, -5.1475,  -1, 
                            -9.4459,  -1, -7.2404,   2, -5.2151,   1, -7.1778, 1.5, 
                            -2.6523,   0, -9.4459,  -1, -7.2404, 1.5, -5.2151,   1)

b_chain = c(mcmc_out_temp$B_chain[nrow(mcmc_out_temp$B_chain), ])

print("initial state sequence based on:")
print(prev_file)

pcov = mcmc_out_temp$pcov
pscale = mcmc_out_temp$pscale

rm(mcmc_out_temp)
# -----------------------------------------------------------------------------
A = list()
W = list()
B = list()

for(i in EIDs){

    A[[i]] = matrix(par[par_index$vec_alpha_tilde], ncol =1)
    W[[i]] = matrix(par[par_index$omega_tilde], ncol =1)
    
    b_temp = rep(1, sum(Y[,"EID"] == as.numeric(i)))
    if(sum(z[Y[,"EID"] == as.numeric(i), 2]) >= 3) {
        sub_z = z[Y[,"EID"] == as.numeric(i), 2]
        s_bleed = min(which(sub_z != 0))
        if(s_bleed == 1) {
            b_temp = c(2,2,3,3)
        } else if(s_bleed == 2) {
            b_temp = c(1,2,2,3)
        } else if(s_bleed == 3) {
            b_temp = c(1,1,2,2)
        } else {
            b_temp = c(1,1,1,2)
        }
    }
    
    B[[i]] = matrix(b_temp, ncol = 1)
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
colnames(zed) = c('1->2', '1->4', '2->3', '2->4', '3->1', '3->2', '3->4', '4->2', 
                  '4->5', '5->1', '5->2', '5->4')
print(zed)



n_cores = 1#strtoi(Sys.getenv(c("LSB_DJOB_NUMPROC")))
print(paste0("Number of cores: ", n_cores))

# Start with no RBC or clinic rule
Y[,"RBC_rule"] = 0

# Initial otype
otype = !is.na(Y[, c('hemo','hr','map','lactate')])
colnames(otype) = c('hemo','hr','map','lactate')

# Change inputs every minute and save the resulting mcmc_out output ------------
it_num = 1
while(nrow(Y) > 0) {
    print(paste0("iteration = ", it_num))
    # Initialize Y -------------------------------------------------------------
    if(it_num == 1) {
        Y = initialize_Y(Y, EIDs)
    } else {
        na_hemo = which(is.na(Y[,"hemo"]))
        if(length(na_hemo) > 0) Y[na_hemo,"hemo"] = Y[na_hemo-1,"hemo"]
        
        na_hr = which(is.na(Y[,"hr"]))
        if(length(na_hr) > 0) Y[na_hr,"hr"] = Y[na_hr-1,"hr"]
        
        na_map = which(is.na(Y[,"map"]))
        if(length(na_map) > 0) Y[na_map,"map"] = Y[na_map-1,"map"]
        
        na_lact = which(is.na(Y[,"lactate"]))
        if(length(na_lact) > 0) Y[na_lact,"lactate"] = Y[na_lact-1,"lactate"]
    }
    
    # Run procedure ------------------------------------------------------------
    s_time = proc.time();
    mcmc_out = mcmc_online( par, par_index, A, W, B, Y, x, z, steps, burnin, ind, 
                            trialNum, Dn_omega, simulation, max_ind, df_num, 
                            n_cores, otype, pcov, pscale, s_time)
    
    # Plot output / update chart plots -----------------------------------------
    B_chain   = mcmc_out$B_chain[1:mcmc_out$ttt, ]
    Hr_chain  = mcmc_out$hr_chain[1:mcmc_out$ttt, ]
    Map_chain = mcmc_out$bp_chain[1:mcmc_out$ttt, ]
    Hc_chain  = mcmc_out$hc_chain[1:mcmc_out$ttt, ]
    La_chain  = mcmc_out$la_chain[1:mcmc_out$ttt, ]
    if(it_num == 1) {
        prev_state_avg = rbind( colMeans(B_chain == 1), colMeans(B_chain == 2),
                                colMeans(B_chain == 3), colMeans(B_chain == 4),
                                colMeans(B_chain == 5))
        prev_state_avg = rbind(Y[,'EID'], prev_state_avg)
    }
    
    if(it_num == 1 | it_num %% 5 == 0) {
        plotting_fnc(df_sub_start, EIDs, B_chain, Hc_chain, Hr_chain, Map_chain, 
                     La_chain, it_num, prev_state_avg)   
    }
    
    # Keep track of previous posterior probabilities ---------------------------
    if(it_num > 1) {
        prev_state_avg = rbind( colMeans(B_chain == 1), colMeans(B_chain == 2),
                                colMeans(B_chain == 3), colMeans(B_chain == 4),
                                colMeans(B_chain == 5))
        prev_state_avg = rbind(Y[,'EID'], prev_state_avg)
    }
    
    # Continue to add data to Y ------------------------------------------------
    Y_curr_hemo = colMeans(Hc_chain)
    Y_curr_hr = colMeans(Hr_chain)
    Y_curr_map = colMeans(Map_chain)
    Y_curr_lact = colMeans(La_chain)
    Y_curr_EIDs = Y[,'EID']
    b_curr = B_chain[nrow(B_chain), ]
    
    # Reinitialize all of the other variables ----------------------------------
    new_data = add_data(it_num, pred_id, data_format, b_curr, Y_curr_EIDs)
    Y = new_data$Y
    if(nrow(Y) == 0) break; # No more updates
    
    EIDs = new_data$EIDs
    x = new_data$x
    z = new_data$z
    otype = new_data$otype
    
    A = new_data$A
    W = new_data$W
    B = new_data$B
    df_sub_start = new_data$df_sub_start
    
    
    # Initialize missing Y values ----------------------------------------------
    for(e in EIDs) {
        row_e = which(Y[,"EID"] == e)
        Y[row_e[-length(row_e)], "hemo"] = Y_curr_hemo[Y_curr_EIDs==e]
        Y[row_e[-length(row_e)], "hr"] = Y_curr_hr[Y_curr_EIDs==e]
        Y[row_e[-length(row_e)], "map"] = Y_curr_map[Y_curr_EIDs==e]
        Y[row_e[-length(row_e)], "lactate"] = Y_curr_lact[Y_curr_EIDs==e]
    }
    
    it_num = it_num + 1
}

