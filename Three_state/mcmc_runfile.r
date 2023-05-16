source('mcmc_routine_arma.r')

args = commandArgs(TRUE)

ind = as.numeric(args[1])
set.seed(ind)
print(ind)

steps  = 10000
burnin =  5000

load("Data/data_format_FULL_48hr_update_RBC_sub.rda")
load('Data/Dn_omega.rda')

# Removing pacing patients
pace_id = c(18075, 108825, 110750, 125025, 173750, 260100, 304700, 307225, 310100,
            382450, 429375, 516150, 533075, 666750, 677225, 732525, 763050, 767500, 
            769025, 777175, 794900, 799125, 819225)
data_format = data_format[!(data_format[,'EID'] %in% pace_id), ]

# Manually perturbing instances to see if things change
# 100950
# data_format[which(data_format[,"EID"] == 100950)[15:18], 'map'] = 51
# data_format[which(data_format[,"EID"] == 100950)[18], 'hemo'] = 6

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

# theta = exp(par_means[par_index$log_theta])
theta = 1

# columns: hemo, hr, map, lactate
R = matrix( c(    .2, -.1,  .1, -.1,
                 -.1,   2, -.1,  .1,
                  .1, -.1,   2, -.1,
                 -.1,  .1, -.1,  .2), ncol=4, byrow=TRUE)

# transitions: 1->2, 2->3, 3->1, 3->2
zeta = matrix(c(-5.236006, -3.078241,        -4,     -5.23,
                 2.006518, -1.688983, -0.056713,  2.044297), nrow = 2, byrow = T)

omega = c(3, -3,   -3, 3,   3, -3,   -3, 3)
upsilon_omega = c(diag(8))

# init_logit = c(0, par_means[par_index$vec_init])
init_logit = c(0,-5,0.5)
init_logit = exp(init_logit)

par = c(beta, c(alpha_tilde), c(sigma_upsilon), c(log(theta)), c(R), c(zeta), 
        c(-5,0.5), log(diag(Lambda)), omega, upsilon_omega)
par_index = list()
par_index$vec_beta = 1:4
par_index$vec_alpha_tilde = 5:16
par_index$vec_sigma_upsilon = 17:160
par_index$log_theta = 161
par_index$vec_R = 162:177
par_index$vec_zeta = 178:185
par_index$vec_init = 186:187
par_index$log_lambda = 188:199
par_index$omega_tilde = 200:207
par_index$vec_upsilon_omega = 208:271
# -----------------------------------------------------------------------------

A = list()
W = list()
B = list()
for(i in EIDs){
  A[[i]] = c(alpha_tilde)
  W[[i]] = rep(0, length(omega))
  
  temp = data_format[data_format[,'EID']==as.numeric(i), ]
  b_temp = matrix( 3, sum(Y[,'EID']==as.numeric(i)), 1)
  
  b_length = nrow(b_temp)
  b_temp[(b_length-5):b_length, ] = 1
  
  B[[i]] = b_temp
}

trialNum = 6 # CHANGE THIS EVERY TIME **********************

# -----------------------------------------------------------------------------
load('Model_out/mcmc_out_interm_4_4it5.rda')
par_temp = colMeans(mcmc_out_temp$chain)
rownames(par_temp) = NULL
par = par_temp

rm(mcmc_out_temp)
# -----------------------------------------------------------------------------

print(par)

s_time = Sys.time()
mcmc_out = mcmc_routine( par, par_index, A, W, B, Y, x, z, steps, burnin, ind, trialNum, Dn_omega)
e_time = Sys.time() - s_time; print(e_time)

# save( mcmc_out, file=paste0('Model_out/post_mcmc_out_dev',ind,'_', trialNum, '.rda'))
