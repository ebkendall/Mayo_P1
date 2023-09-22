source('mcmc_routine_arma.r')

args = commandArgs(TRUE)

ind = as.numeric(args[1])
set.seed(ind)
print(ind)

simulation = F
data_num = 3

steps  = 50000
burnin =  5000

data_format = NULL
if(simulation) {
  load(paste0('Data/use_data1_', data_num, '.rda'))
  data_format = use_data
  trialNum = 9
} else {
  load('Data/data_format_new.rda')
  pace_id = c(53475, 110750, 125025, 260625, 273425, 296500, 310100, 384925,
              417300, 448075, 538075, 616025, 660075, 665850, 666750, 677225,
              732525, 758025, 763050, 843000)
  data_format = data_format[!(data_format[,'EID'] %in% pace_id), ]
  trialNum = 2 # CHANGE THIS EVERY TIME **********************
}

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
Lambda = diag(c(   2,.1,.1,   3,.1,.1,   4,.25,.25,  2,.1,.1))
Upsilon = Lambda %*% sigma_upsilon %*% Lambda

# columns correspond to the different states
# Each column corresponds to a different state
vec_A = c(matrix(c(2.3, 0.5, 1.3,
                   1.9,  -1, 0.5,
                   1.9,  -1, 0.5,
                   2.3, 0.5, 1.3), ncol = 3, byrow = T))

# columns: hemo, hr, map, lactate
# diagonal elements only (because R is a diagonal matrix)
R = diag(1, 4)

# transitions: 1->2, 2->3, 3->1, 3->2
zeta = matrix(c(-5.236006, -3.078241,        -4,     -5.23,
                 2.006518, -1.688983, -0.056713,  2.044297), nrow = 2, byrow = T)

omega = c(3, -3,   -3, 3,   3, -3,   -3, 3)
upsilon_omega = c(diag(8))

init_logit = c(-5,-5)
init_logit = exp(init_logit)

par = c(beta, c(alpha_tilde), c(sigma_upsilon), c(vec_A), c(R), c(zeta), 
        log(init_logit), log(diag(Lambda)), omega, upsilon_omega)
par_index = list()
par_index$vec_beta = 1:4
par_index$vec_alpha_tilde = 5:16
par_index$vec_sigma_upsilon = 17:160
par_index$vec_A = 161:172
par_index$vec_R = 173:188
par_index$vec_zeta = 189:196
par_index$vec_init = 197:198
par_index$log_lambda = 199:210
par_index$omega_tilde = 211:218
par_index$vec_upsilon_omega = 219:282
# -----------------------------------------------------------------------------

A = list()
W = list()
B = list()
Dn_omega = list()

load(paste0('Data/true_pars_', data_num, '.rda'))
load(paste0('Data/alpha_i_mat_', data_num, '.rda'))

for(i in EIDs){
  W[[i]] = rep(0, length(omega))
  Dn_omega[[i]] = diag(4)
  
  if(simulation) {
      A[[i]] = alpha_i_mat[[which(EIDs == i)]]
      B[[i]] = data_format[data_format[,'EID']==as.numeric(i), "b_true", drop=F]
  } else {
      b_temp = matrix( 1, sum(Y[,'EID']==as.numeric(i)), 1)
      B[[i]] = b_temp
      A[[i]] = matrix(par[par_index$vec_alpha_tilde], ncol =1)
  }
}

# -----------------------------------------------------------------------------
if(simulation) {
  load(paste0('Data/true_pars_', data_num, '.rda'))
  par[1:210] = true_pars

  su = matrix(par[par_index$vec_sigma_upsilon], ncol = 12)
  la = diag(exp(par[par_index$log_lambda]))
  up_true = la %*% su %*% la

  par[par_index$vec_sigma_upsilon] = c(up_true)
  par[par_index$log_lambda] = 0
} else {
  # load('Model_out/mcmc_out_interm_5_6it3.rda')
  # par_temp = colMeans(mcmc_out_temp$chain)
  # rownames(par_temp) = NULL
  # par[1:172] = par_temp[1:172]
  # par[189:282] = par_temp[177:270]
  # rm(mcmc_out_temp) 
}
# -----------------------------------------------------------------------------

print("Starting values for the chain")
print("alpha_tilde")
print(round(par[par_index$vec_alpha_tilde], 3))

print("mean alpha_i")
A_stacked = do.call( cbind, A)
print(apply(A_stacked, 1, mean))

print("single alpha_i")
a_ind = sample(size = 1, 1:180)
print(a_ind)
print(c(A[[a_ind]]))

print("diag of Sigma_Upsilon")
Sigma_t = matrix(par[par_index$vec_sigma_upsilon], ncol = 12)
# Lambda_t = diag(exp(par[par_index$log_lambda]))
# Upsilon_t = Lambda_t %*% Sigma_t %*% Lambda_t
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

s_time = Sys.time()
mcmc_out = mcmc_routine( par, par_index, A, W, B, Y, x, z, steps, burnin, ind, trialNum, Dn_omega, simulation)
e_time = Sys.time() - s_time; print(e_time)
