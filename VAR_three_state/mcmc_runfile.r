source('mcmc_routine_arma.r')

args = commandArgs(TRUE)

ind = as.numeric(args[1])
set.seed(ind)
print(ind)

steps  = 20000
burnin =  5000

load('Data/data_format_new.rda')
# load('Data/Dn_omega.rda')

Y = data_format[, c('EID','hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
EIDs = as.character(unique(data_format[,'EID']))

x = data_format[,c('n_RBC_admin'), drop=F]
p = ncol(x)

z = cbind(1, data_format[,c('RBC_ordered'), drop=F])
m = ncol(z)

# Parameters ------------------------------------------------------------------
beta = c(0.6261, -1.3286, 1.6741, -0.1)

# columns: hemo, hr, map
alpha_tilde = matrix( c( 9.57729783, 88.69780576, 79.74903940,  5.2113319,
                                 -1,  9.04150472, -7.42458547,  0.5360813,
					                      0.1,          -4,           4, -0.6866748), ncol=4, byrow=T)

sigma_upsilon = diag(12)
# Lambda = diag(c(   2,.1,.1,   3,.1,.1,   4,.25,.25,  2,.1,.1))
Lambda = diag(rep(1, 12))
Upsilon = Lambda %*% sigma_upsilon %*% Lambda

# columns correspond to the different states
# Each column corresponds to a different state
# vec_A = matrix( 0 , nrow = 4, ncol = 3) 
vec_A = rep(0, 4)

# columns: hemo, hr, map, lactate
# diagonal elements only (because R is a diagonal matrix)
R = rep(2, 4)

# transitions: 1->2, 2->3, 3->1, 3->2
zeta = matrix(c(-5.236006, -3.078241,        -4,     -5.23,
                 2.006518, -1.688983, -0.056713,  2.044297), nrow = 2, byrow = T)

omega = c(3, -3,   -3, 3,   3, -3,   -3, 3)
upsilon_omega = c(diag(8))

init_logit = c(-5,0.5)
init_logit = exp(init_logit)

par = c(beta, c(alpha_tilde), c(sigma_upsilon), c(vec_A), log(R), c(zeta), 
        log(init_logit), log(diag(Lambda)), omega, upsilon_omega)
par_index = list()
par_index$vec_beta = 1:4
par_index$vec_alpha_tilde = 5:16
par_index$vec_sigma_upsilon = 17:160
par_index$vec_logit_A = 161:164
par_index$vec_R = 165:168
par_index$vec_zeta = 169:176
par_index$vec_init = 177:178
par_index$log_lambda = 179:190
par_index$omega_tilde = 191:198
par_index$vec_upsilon_omega = 199:262
# -----------------------------------------------------------------------------

A = list()
W = list()
B = list()
Dn_omega = list()
for(i in EIDs){
  A[[i]] = c(alpha_tilde)
  W[[i]] = rep(0, length(omega))
  Dn_omega[[i]] = diag(4)
  
  temp = data_format[data_format[,'EID']==as.numeric(i), ]
  b_temp = matrix( 1, sum(Y[,'EID']==as.numeric(i)), 1)
  
  # b_length = nrow(b_temp)
  # b_temp[(b_length-5):b_length, ] = 1
  
  B[[i]] = b_temp
}

trialNum = 1 # CHANGE THIS EVERY TIME **********************

# -----------------------------------------------------------------------------
# load('Model_out/mcmc_out_interm_1_5it7.rda')
# par_temp = colMeans(mcmc_out_temp$chain)
# rownames(par_temp) = NULL
# par = par_temp
# par[par_index$log_theta] = 9
# par[par_index$vec_R] = par[par_index$vec_R] * 25
# rm(mcmc_out_temp)
# -----------------------------------------------------------------------------

s_time = Sys.time()
mcmc_out = mcmc_routine( par, par_index, A, W, B, Y, x, z, steps, burnin, ind, trialNum, Dn_omega)
e_time = Sys.time() - s_time; print(e_time)