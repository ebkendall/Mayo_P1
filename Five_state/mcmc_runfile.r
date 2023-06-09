source('mcmc_routine_arma.r')

args = commandArgs(TRUE)

ind = as.numeric(args[1])
set.seed(ind)
print(ind)

n_cores = 20

steps  = 100000
burnin =  5000

load("Data/data_format_FULL_48hr_update_RBC_sub.rda")

# Removing pacing patients
pace_id = c(18075, 108825, 110750, 125025, 173750, 260100, 304700, 307225, 310100,
            382450, 429375, 516150, 533075, 666750, 677225, 732525, 763050, 767500, 
            769025, 777175, 794900, 799125, 819225)
data_format = data_format[!(data_format[,'EID'] %in% pace_id), ]

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
					            0.1,          -2,           2, -0.6866748,
					              0,           0,           0,          0,
					              0,           0,           0,          0), ncol=4, byrow=T)

sigma_upsilon = diag(20)
Lambda = diag(c(   2,.1,.1,.1,.1,   3,.1,.1,.1,.1,   4,.25,.25,.25,.25,  2,.1,.1,.1,.1))
Upsilon = Lambda %*% sigma_upsilon %*% Lambda

theta = 1

# columns: hemo, hr, map, lactate
R = matrix( c(    .2, -.1,  .1, -.1,
                 -.1,   2, -.1,  .1,
                  .1, -.1,   2, -.1,
                 -.1,  .1, -.1,  .2), ncol=4, byrow=TRUE)

# transitions: 1->2,1->4, 2->3, 2->4, 3->1, 3->2, 3->4, 4->2, 4->5, 5->1, 5->2, 5->4
zeta = matrix( c(-6.1708, -6.1708, -2.652, -2.652, -3, -3, -3, -2.652, -2.652, -3, -3, -3,
                       2,       0,     -2,     -2,  0,  2,  0,      2,      0,  0,  2,  0),
               ncol=12, byrow=T) 


# init_logit = c(0, par_means[par_index$vec_init])
init_logit = c(0,-5,0.5,-5,-5)
init_logit = exp(init_logit)

par = c(beta, c(alpha_tilde), c(sigma_upsilon), c(log(theta)), c(R), c(zeta), c(-5,0.5,-5,-5), log(diag(Lambda)))
par_index = list()
par_index$vec_beta = 1:4
par_index$vec_alpha_tilde = 5:24
par_index$vec_sigma_upsilon = 25:424
par_index$log_theta = 425
par_index$vec_R = 426:441
par_index$vec_zeta = 442:465
par_index$vec_init = 466:469
par_index$log_lambda = 470:489
# -----------------------------------------------------------------------------

A = list()
B = list()
for(i in EIDs){
  A[[i]] = c(alpha_tilde)
  
  temp = data_format[data_format[,'EID']==as.numeric(i), ]
  b_temp = matrix( 1, sum(Y[,'EID']==as.numeric(i)), 1)
  
  b_length = nrow(b_temp)
  b_temp[(b_length-5):b_length, ] = 1
  
  B[[i]] = b_temp
}

trialNum = 1 # CHANGE THIS EVERY TIME **********************

# -----------------------------------------------------------------------------
# index_post = 8000:10000
# load(paste0('Model_out/mcmc_out_interm_', 5, '_', trialNum - 1,'it5.rda'))
# par_temp = colMeans(mcmc_out_temp$chain[index_post,])
# par_temp = colMeans(mcmc_out_temp$chain)
# rownames(par_temp) = NULL
# par = par_temp
# 
# rm(mcmc_out_temp)
# -----------------------------------------------------------------------------

print(par)

s_time = Sys.time()
mcmc_out = mcmc_routine( par, par_index, A, B, Y, x, z, steps, burnin, n_cores, ind, trialNum)
e_time = Sys.time() - s_time; print(e_time)

# save( mcmc_out, file=paste0('Model_out/post_mcmc_out_dev',ind,'_', trialNum, '.rda'))
