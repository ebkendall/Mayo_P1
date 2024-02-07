source('mcmc_routine_arma.r')

args = commandArgs(TRUE)

ind = as.numeric(args[1])
set.seed(ind)
print(ind)

# Simulation trials:
# trial 3: fix everything but VAR parameters
# trial 4: run with state space fixed
# trial 5: fix everything but VAR pars and state space

# trial 1: Full sim with updated data (start at correct states, no state update)
# trial 2: Full sim with updated data (start at correct states)

# trial 7 start from same values, trial 8 - on are continuations of the previous 
# iteration

simulation = T
sim_dat_num = 6
real_dat_num = 3

data_format = NULL
if(simulation) {
    steps  = 50000
    burnin =  5000
    
    load(paste0('Data/use_data1_', sim_dat_num, '.rda'))
    data_format = use_data
    trialNum = 2
} else {
    steps  = 50000
    burnin =  0
    
    load(paste0('Data/data_format_new', real_dat_num, '.rda'))
    trialNum = 7
}

Y = data_format[, c('EID','hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
EIDs = as.character(unique(data_format[,'EID']))

x = data_format[,c('n_RBC_admin'), drop=F]
p = ncol(x)

z = cbind(1, data_format[,c('RBC_ordered'), drop=F])
m = ncol(z)

# Parameters ------------------------------------------------------------------
beta = rep(0, 4)

# columns: hemo, hr, map
alpha_tilde = matrix( c( 9.57729783, 88.69780576, 79.74903940,  5.2113319,
                                 -1,  9.04150472, -7.42458547,  0.5360813,
					            0.1,          -4,           4, -0.6866748), 
					            ncol=4, byrow=T)

sigma_upsilon = diag(c(4, 0.25, 0.25, 36, 1, 1, 64, 1, 1, 4, 0.25, 0.25))
Lambda = diag(12)

vec_A = rep(0, 12)

# columns: hemo, hr, map, lactate
R = diag(4)

# transitions: 1->2, 2->3, 3->1, 3->2
zeta = matrix(c(-5.236006, -3.078241,        -4,     -5.23,
                 2.006518, -1.688983, -0.056713,  2.044297), nrow = 2,byrow = T)

omega = c( 1,  1, -1,  1, -1, -1,  1,  1, -1, -1,  1,  1,  1, -1,  1,
          -1, -1, -1, -1, -1, -1, -1,  1,  1, -1, -1,  1, -1, -1,  1, 
          -1,  1, -1, -1, -1, -1,  1, -1,  1, -1,  1, -1, -1, -1, -1,  
           1, -1, -1,  1, -1,  1, -1,  1,  1, -1, -1, -1, -1,  1, -1, 
          -1, -1,  1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1, -1,
          -1, -1, -1, -1, -1,  1,  1, -1, -1, -1, -1, -1, -1)

omega = 6 * omega
upsilon_omega = rep(1, length(omega))

init_logit = c(-5,-5)
init_logit = exp(init_logit)

par = c(beta, c(alpha_tilde), c(sigma_upsilon), c(vec_A), c(R), c(zeta), 
        log(init_logit), omega, log(upsilon_omega))
par_index = list()
par_index$vec_beta = 1:4
par_index$vec_alpha_tilde = 5:16
par_index$vec_sigma_upsilon = 17:160
par_index$vec_A = 161:172
par_index$vec_R = 173:188
par_index$vec_zeta = 189:196
par_index$vec_init = 197:198
par_index$omega_tilde = 199:286
par_index$vec_upsilon_omega = 287:374
# -----------------------------------------------------------------------------

if(simulation) {
    load(paste0('Data/true_pars_', sim_dat_num, '.rda'))
    load(paste0('Data/alpha_i_mat_', sim_dat_num, '.rda'))
    load(paste0('Data/omega_i_mat_', sim_dat_num, '.rda'))
    load(paste0('Data/Dn_omega_sim_', sim_dat_num, '.rda'))
    load(paste0('Data/bleed_indicator_sim_', sim_dat_num,'.rda'))
    
    par = true_pars
    Dn_omega = Dn_omega_sim
} else {
    load('Model_out/mcmc_out_interm_2_6it4.rda')
    load(paste0('Data/Dn_omega', real_dat_num, '.rda'))
    bleed_indicator = b_ind_fnc(data_format)
    par_temp = mcmc_out_temp$chain[1001, ]
    rownames(par_temp) = NULL

    par = par_temp
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
      b_temp = mcmc_out_temp$B_chain[1001, Y[,'EID']==as.numeric(i)]
    #   b_temp = rep( 1, sum(Y[,'EID']==as.numeric(i)))
      B[[i]] = matrix(b_temp, ncol = 1)
      A[[i]] = matrix(par[par_index$vec_alpha_tilde], ncol =1)
      W[[i]] = matrix(par[par_index$omega_tilde], ncol =1)
  }
}
if(!simulation) rm(mcmc_out_temp)
# -----------------------------------------------------------------------------

print("Starting values for the chain")
print("alpha_tilde")
print(round(par[par_index$vec_alpha_tilde], 3))

print("diag of Sigma_Upsilon")
Sigma_t = matrix(par[par_index$vec_sigma_upsilon], ncol = 12)
print(round(diag(Sigma_t), 3))

print("A1")
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

print("omega_tilde")
print(par[par_index$omega_tilde])

print("log upsilon omega")
print(par[par_index$vec_upsilon_omega])

s_time = Sys.time()
mcmc_out = mcmc_routine( par, par_index, A, W, B, Y, x, z, steps, burnin, ind, 
                         trialNum, Dn_omega, simulation, bleed_indicator)
e_time = Sys.time() - s_time; print(e_time)
