source('mcmc_routine_arma.r')

args = commandArgs(TRUE)

ind = as.numeric(args[1])
set.seed(ind)
print(ind)

simulation = T
data_num = 2

steps  = 50000
burnin =  5000

data_format = NULL
if(simulation) {
  load(paste0('Data/use_data', ind,'_', data_num, '.rda'))
  print(paste0('Data/use_data', ind,'_', data_num, '.rda'))
  data_format = use_data
  trialNum = 2
} else {
  load('Data/data_format_new2.rda')
  trialNum = 1 # CHANGE THIS EVERY TIME **********************
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
					                      0.1,          -4,           4, -0.6866748), ncol=4, byrow=T)

sigma_upsilon = diag(c(4, 0.25, 0.25, 36, 1, 1, 64, 1, 1, 4, 0.25, 0.25))

vec_A = rep(0, 4)

# columns: hemo, hr, map, lactate
R_base = diag(4)

R_bleed = diag(4)


# transitions: 1->2, 2->3, 3->1, 3->2
zeta = matrix(c(-5.236006, -3.078241,        -4,     -5.23,
                 2.006518, -1.688983, -0.056713,  2.044297), nrow = 2, byrow = T)

omega = c(3, -3,   -3, 3,   3, -3,   -3, 3)
upsilon_omega = c(diag(8))

init_logit = c(-5,-5)
init_logit = exp(init_logit)

par = c(beta, c(alpha_tilde), c(sigma_upsilon), c(vec_A), c(R_base), c(zeta), 
        log(init_logit), c(R_bleed), omega, upsilon_omega)
par_index = list()
par_index$vec_beta = 1:4
par_index$vec_alpha_tilde = 5:16
par_index$vec_sigma_upsilon = 17:160
par_index$vec_A = 161:164
par_index$vec_R = 165:180
par_index$vec_zeta = 181:188
par_index$vec_init = 189:190
par_index$vec_R_bleed = 191:206
par_index$omega_tilde = 207:214
par_index$vec_upsilon_omega = 215:278
# -----------------------------------------------------------------------------
load(paste0('Data/true_pars_', data_num, '.rda'))
# load(paste0('Data/alpha_i_mat_', data_num, '.rda'))

if(simulation) {
    par[1:206] = true_pars
} else {
    # par_temp = colMeans(mcmc_out_temp$chain)
    # rownames(par_temp) = NULL
    # par = par_temp
}
# -----------------------------------------------------------------------------
A = list()
W = list()
B = list()
Dn_omega = list()

for(i in EIDs){
  W[[i]] = rep(0, length(omega))
  Dn_omega[[i]] = diag(4)
  
  if(simulation) {
      A[[i]] = matrix(par[par_index$vec_alpha_tilde], ncol = 1)
      B[[i]] = data_format[data_format[,'EID']==as.numeric(i), "b_true", drop=F]
  } else {
      b_temp = rep( 1, sum(Y[,'EID']==as.numeric(i)))
      # b_temp = mcmc_out_temp$B_chain[1001, Y[,'EID']==as.numeric(i)]
      B[[i]] = matrix(b_temp, ncol = 1)
      A[[i]] = matrix(par[par_index$vec_alpha_tilde], ncol = 1)
  }
}
# rm(mcmc_out_temp)
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

print("R_base")
R_t = matrix(par[par_index$vec_R], ncol = 4)
print(R_t)

print("R_bleed")
R_b = matrix(par[par_index$vec_R_bleed], ncol = 4)
print(R_b)

print("zeta")
zed = matrix(par[par_index$vec_zeta], nrow = 2)
print(zed)

s_time = Sys.time()
mcmc_out = mcmc_routine( par, par_index, A, W, B, Y, x, z, steps, burnin, ind, trialNum, Dn_omega, simulation)
e_time = Sys.time() - s_time; print(e_time)
