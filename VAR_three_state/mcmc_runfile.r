source('mcmc_routine_arma.r')

args = commandArgs(TRUE)

ind = as.numeric(args[1])
set.seed(ind)
print(ind)

simulation = T

steps  = 30000
burnin =  5000

data_format = NULL
if(simulation) {
  load('Data/use_data1_1.rda')
  data_format = use_data
  trialNum = 3
} else {
  load('Data/data_format_new.rda')
  pace_id = c(53475, 110750, 125025, 260625, 273425, 296500, 310100, 384925,
              417300, 448075, 538075, 616025, 660075, 665850, 666750, 677225,
              732525, 758025, 763050, 843000)
  data_format = data_format[!(data_format[,'EID'] %in% pace_id), ]
  trialNum = 7 # CHANGE THIS EVERY TIME **********************
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
vec_A = matrix( 2.3 , nrow = 4, ncol = 3) 

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
for(i in EIDs){
  load('Data/true_pars_1.rda')
  A[[i]] = matrix(true_pars[par_index$vec_alpha_tilde], ncol =1)
  W[[i]] = rep(0, length(omega))
  Dn_omega[[i]] = diag(4)
  
  if(simulation) {
      B[[i]] = data_format[data_format[,'EID']==as.numeric(i), "b_true", drop=F]
  } else {
      b_temp = matrix( 1, sum(Y[,'EID']==as.numeric(i)), 1)
      # b_length = nrow(b_temp)
      # b_temp[(b_length-5):b_length, ] = 1
      B[[i]] = b_temp
  }
}

# -----------------------------------------------------------------------------
if(simulation) {
  load('Data/true_pars_1.rda')
  par[1:210] = true_pars
} else {
  load('Model_out/mcmc_out_interm_5_6it3.rda')
  par_temp = colMeans(mcmc_out_temp$chain)
  rownames(par_temp) = NULL
  par[1:172] = par_temp[1:172]
  par[189:282] = par_temp[177:270]
  rm(mcmc_out_temp) 
}
# -----------------------------------------------------------------------------

s_time = Sys.time()
mcmc_out = mcmc_routine( par, par_index, A, W, B, Y, x, z, steps, burnin, ind, trialNum, Dn_omega, simulation)
e_time = Sys.time() - s_time; print(e_time)
