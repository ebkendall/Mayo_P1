library(mvtnorm, quietly=T)

# Load in the existing data and save the covariate combinations
real_dat_num = 3
load(paste0('Data/data_format_new', real_dat_num, '.rda'))

it_num = 6
set.seed(2018)
N = length(unique(data_format[,"EID"]))

# Making an indicator variable about the first RBC to indicate bleed event
bleed_pat = unique(data_format[data_format[,"RBC_rule"] != 0, "EID"])
bleed_indicator = rep(0, nrow(data_format))
for(i in 1:length(bleed_pat)) {
    sub_dat = data_format[data_format[,"EID"] == bleed_pat[i], ]
    
    # Check in any 12 hour period
    max_time = tail(sub_dat[,"time"], 1)
    when_rbc = c(1, which(diff(sub_dat[,"n_RBC_admin"]) != 0))
    
    for(j in 1:length(when_rbc)) {
        s_time = sub_dat[when_rbc[j], "time"]
        e_time_12 = s_time + 720
        e_time_24 = s_time + 1440
        RBC_diff_12 = RBC_diff_24 = 0
        
        if (e_time_12 <= max_time) {
            s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
            ind_12 = order(abs(sub_dat[,"time"] - e_time_12))[1]
            RBC_diff_12 = sub_dat[ind_12, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
        } else {
            s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
            e_ind = order(abs(sub_dat[,"time"] - max_time))[1]
            RBC_diff_12 = sub_dat[e_ind, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
        }
        if (e_time_24 <= max_time) {
            s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
            ind_24 = order(abs(sub_dat[,"time"] - e_time_24))[1]
            RBC_diff_24 = sub_dat[ind_24, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
        } else {
            s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
            e_ind = order(abs(sub_dat[,"time"] - max_time))[1]
            RBC_diff_24 = sub_dat[e_ind, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
        }
        
        if(RBC_diff_12 >=3 | RBC_diff_24 >= 6) {
            admin_times = sub_dat[sub_dat[,"RBC_admin"] != 0, "time"]
            if(RBC_diff_12 >=3) {
                a_t = which(admin_times >= s_time & admin_times < e_time_12)
                first_time = admin_times[a_t[1]]
                order_times = sub_dat[sub_dat[,"RBC_ordered"] != 0, "time"]
                if(sum(order_times <= first_time) == 0) {
                    print(paste0(i, ", ", sub_dat[1,"EID"]))
                    first_order_time = first_time
                } else {
                    first_order_time = max(order_times[order_times <= first_time])   
                }
            } else if (RBC_diff_24 >= 6) {
                a_t = which(admin_times >= s_time & admin_times < e_time_24)
                first_time = admin_times[a_t[1]]  
                order_times = sub_dat[sub_dat[,"RBC_ordered"] != 0, "time"]
                if(sum(order_times <= first_time) == 0) {
                    print(paste0(i, ", ", sub_dat[1,"EID"]))
                    first_order_time = first_time
                } else {
                    first_order_time = max(order_times[order_times <= first_time])   
                }
            }

            bleed_indicator[data_format[,"EID"] == bleed_pat[i] & 
                                data_format[,"time"] == first_order_time] = 1
            break
        }
        
    }
    
}
# save(bleed_indicator, file = 'Data/bleed_indicator_real.rda')
# -----------------------------------------------------------------------------

data_format = cbind(data_format, bleed_indicator)

load(paste0('Data/Dn_omega', real_dat_num, '.rda'))

Y = data_format[, c('EID','time','hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
EIDs = as.character(unique(data_format[,'EID']))

x = data_format[,c('n_RBC_admin'), drop=F]
p = ncol(x)

z = cbind(1, data_format[,c('RBC_ordered'), drop=F])
m = ncol(z)

# Loading the parameter values to base this off of
load('Model_out/mcmc_out_interm_2_6it4.rda')
pars_mean = colMeans(mcmc_out_temp$chain[900:1001,])

# Initializing par_index
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

save(par_index, file = paste0('Data/true_par_index_', it_num, '.rda'))

# Parameters ------------------------------------------------------------------
# beta = c(0.6261, -1.3286, 1.6741, -0.1)
beta = pars_mean[par_index$vec_beta]

# alpha_tilde = matrix( c( 9.57729783, 88.69780576, 79.74903940, 5.2113319,
#                                  -1,  5.04150472, -5.42458547,         1,
#                                   1, 	      -4, 		    4,        -1), ncol=4, byrow=T)
alpha_tilde = matrix(pars_mean[par_index$vec_alpha_tilde], ncol = 4)

# sigma_upsilon = Upsilon = diag(c(4, 0.25, 0.25, 36, 1, 1, 36, 1, 1, 4, 0.25, 0.25))
sigma_upsilon = Upsilon = matrix(pars_mean[par_index$vec_sigma_upsilon], ncol = 12)

# A_mat = matrix(c(2, -2, 0,
#                  2, -2, 0,
#                  2, -2, 0,
#                  2, -2, 0), ncol = 3, byrow = T)
A_mat = matrix(pars_mean[par_index$vec_A], ncol = 3)
vec_A = c(A_mat)
correct_scale_A = exp(vec_A) / (1 + exp(vec_A))
A_mat_scale = matrix(correct_scale_A, nrow = 4)

# columns: hemo, hr, map, lactate
# R = diag(4)
R = matrix(pars_mean[par_index$vec_R], ncol = 4)

# transitions: 1->2, 2->3, 3->1, 3->2
# zeta = matrix(c(      -5, -4.078241, -7.000000, -7.230000,
#                 3.006518,      -1.2,   -1.6713,  3.544297), ncol = 4, byrow=T)
zeta = matrix(pars_mean[par_index$vec_zeta], ncol = 4)

# init_logit = c(0,-5,-2)
# init_logit = exp(init_logit)
init_logit = pars_mean[par_index$vec_init]
init_logit = c(0, init_logit)

# omega = c( 1,  1, -1,  1, -1, -1,  1,  1, -1, -1,  1,  1,  1, -1,  1,
#           -1, -1, -1, -1, -1, -1, -1,  1,  1, -1, -1,  1, -1, -1,  1, 
#           -1,  1, -1, -1, -1, -1,  1, -1,  1, -1,  1, -1, -1, -1, -1,  
#            1, -1, -1,  1, -1,  1, -1,  1,  1, -1, -1, -1, -1,  1, -1, 
#           -1, -1,  1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1, -1,
#           -1, -1, -1, -1, -1,  1,  1, -1, -1, -1, -1, -1, -1)
# omega = 6 * omega
omega = pars_mean[par_index$omega_tilde]


# Changing some medication effects
load(paste0('Data/med_select_FINAL', real_dat_num, '.rda'))
load(paste0('Data/Dn_omega_names', real_dat_num, '.rda'))
load(paste0('Data/hr_map_names', real_dat_num, '.rda'))

Dn_omega_sim = vector(mode = 'list', length = N)

hr_max_ind = max(which(hr_map_names == 'hr_disc'))
map_max_ind = max(which(hr_map_names == 'map_disc'))

hr_med = med_select_FINAL$med_name_admin[med_select_FINAL$hr != 0]
map_med = med_select_FINAL$med_name_admin[med_select_FINAL$map != 0]
hr_med_priority = sort(c(table(hr_med)))
map_med_priority = sort(c(table(map_med)))

hr_med_order = sample(1:hr_max_ind, 24)
hr_med_0 = names(hr_med_priority)[hr_med_order[1:12]]
hr_med_1 = names(hr_med_priority)[hr_med_order[13:24]]

map_med_order = sample(1:(map_max_ind - hr_max_ind), 36)
map_med_0 = names(map_med_priority)[map_med_order[1:18]]
map_med_1 = names(map_med_priority)[map_med_order[19:36]]

zero_ind = c(which(Dn_omega_names[1:hr_max_ind] %in% hr_med_0), 
             hr_max_ind + which(Dn_omega_names[(hr_max_ind+1):map_max_ind] %in% map_med_0))
one_ind = c(which(Dn_omega_names[1:hr_max_ind] %in% hr_med_1), 
            hr_max_ind + which(Dn_omega_names[(hr_max_ind+1):map_max_ind] %in% map_med_1))

# omega[zero_ind] = omega[zero_ind] / 6
# omega[one_ind] = omega[one_ind] / 2

# upsilon_omega = runif(length(omega))
upsilon_omega = exp(pars_mean[par_index$vec_upsilon_omega])

# true_pars = c(beta, c(alpha_tilde), c(sigma_upsilon), c(vec_A), c(R), c(zeta), 
#               init_logit[2:3], omega, log(upsilon_omega))
true_pars = pars_mean
save(true_pars, file = paste0('Data/true_pars_', it_num, '.rda'))
# -----------------------------------------------------------------------------

alpha_i_mat = vector(mode = "list", length = N)
omega_i_mat = vector(mode = "list", length = N)
bleed_indicator_update = NULL

for (www in 1:1) {
    
    Dir = 'Data/'
    
    use_data = NULL
    
    rbc_bleed_correct = NULL
    for(i in 1:N){
        
        id_num = EIDs[i]
        ind_i = which(EIDs == id_num)
        Dn_omega_sim[[i]] = Dn_omega[[ind_i]]
        D_i_omega = Dn_omega_sim[[i]]
        
        print(paste0(i, ", ", id_num))
        
        rbc_rule = as.logical(head(Y[Y[,'EID']==as.numeric(id_num),"RBC_rule"], 1))
        correct_bleed = T
        if(rbc_rule) correct_bleed = F
        
        n_i = sum(Y[,'EID']==as.numeric(id_num))
        
        x_i = x[ Y[,'EID']==as.numeric(id_num),, drop=F]
        z_i = z[ Y[,'EID']==as.numeric(id_num),, drop=F]
        bleed_ind_i = bleed_indicator[Y[,'EID']==as.numeric(id_num)]
        
        bleed_indicator_update = c(bleed_indicator_update, bleed_ind_i)
        
        # Generate realizations of latent bleeding process ---------------------
        D_i = vector(mode = 'list', length = n_i)
        X_i = vector(mode = 'list', length = n_i)
        
        if(length(D_i_omega) != n_i) {
            print(paste0("issue n_i: ", i))
        }
        
        P_i = exp(init_logit) / sum(exp(init_logit))
        for(k in 1:n_i){
            if(k==1){
                b_i = sample(1:3, size=1, prob=P_i)
            } else{
                q1   = exp(z_i[k,, drop=F] %*% zeta[,  1, drop=F]) 
                q2   = exp(z_i[k,, drop=F] %*% zeta[,  2, drop=F])
                q3   = exp(z_i[k,, drop=F] %*% zeta[,  3, drop=F])
                q4   = exp(z_i[k,, drop=F] %*% zeta[,  4, drop=F])
                
                # transitions: 1->2, 2->3, 3->1, 3->2
                Q = matrix(c(  1,  q1,  0,
                               0,   1, q2,
                               q3,  q4,  1), ncol=3, byrow=T)
                P_i = Q / rowSums(Q)
                # Sample the latent state sequence
                b_i = c( b_i, sample(1:3, size=1, prob=P_i[tail(b_i,1),]))
            }
            
            D_i_temp = matrix(c( 1, sum(b_i[1:k]==2), sum(b_i[1:k]==3)), nrow = 1, ncol = 3)
            D_i[[k]] = diag(4) %x% D_i_temp
            
            x_i_temp = matrix(c(x_i[k,]), ncol = 1)
            X_i[[k]] = diag(4) %x% x_i[k,]
        }
        
        if(rbc_rule) {
            rbc_bleed_correct = c(rbc_bleed_correct, -1)
            # print(paste0("bleed check: ", i, ", ", id_num))
            if(2 %in% b_i) {
                first_bleed_ind = which(bleed_ind_i == 1)
                sim_bleed_ind = which(b_i == 2)
                if((sum(sim_bleed_ind <= first_bleed_ind) > 0) & 
                   (2 %in% b_i[c(first_bleed_ind, first_bleed_ind - 1)])) {
                    correct_bleed = T
                    rbc_bleed_correct[length(rbc_bleed_correct)] = 1
                } 
            } 
        }
        # ---------------------------------------------------------------------------
        
        # Generate realizations of hc, hr, and bp -----------------------------------
        Y_i = matrix(nrow = n_i, ncol = 4)
        vec_alpha_i = rmvnorm( n=1, mean=c(alpha_tilde), sigma=Upsilon)
        vec_omega_i = rmvnorm( n=1, mean=c(omega), sigma=diag(upsilon_omega))

        alpha_i_mat[[i]] = matrix(vec_alpha_i, ncol = 1)
        omega_i_mat[[i]] = matrix(vec_omega_i, ncol = 1)
        
        for(k in 1:n_i) {
            if(k==1)  {
                A_state_k = A_mat_scale[,b_i[k]]
                Gamma = matrix(c(R[1,1] / (1 - A_state_k[1] * A_state_k[1]), 
                                 R[1,2] / (1 - A_state_k[1] * A_state_k[2]),
                                 R[1,3] / (1 - A_state_k[1] * A_state_k[3]), 
                                 R[1,4] / (1 - A_state_k[1] * A_state_k[4]),
                                 R[2,1] / (1 - A_state_k[2] * A_state_k[1]), 
                                 R[2,2] / (1 - A_state_k[2] * A_state_k[2]),
                                 R[2,3] / (1 - A_state_k[2] * A_state_k[3]), 
                                 R[2,4] / (1 - A_state_k[2] * A_state_k[4]),
                                 R[3,1] / (1 - A_state_k[3] * A_state_k[1]), 
                                 R[3,2] / (1 - A_state_k[3] * A_state_k[2]),
                                 R[3,3] / (1 - A_state_k[3] * A_state_k[3]), 
                                 R[3,4] / (1 - A_state_k[3] * A_state_k[4]),
                                 R[4,1] / (1 - A_state_k[4] * A_state_k[1]), 
                                 R[4,2] / (1 - A_state_k[4] * A_state_k[2]),
                                 R[4,3] / (1 - A_state_k[4] * A_state_k[3]), 
                                 R[4,4] / (1 - A_state_k[4] * A_state_k[4])), 
                               ncol = 4, byrow = T)
                mean_vecY_i_k = D_i[[k]]%*%matrix(vec_alpha_i,ncol=1) + 
                                X_i[[k]]%*%matrix(beta,ncol=1) + 
                                D_i_omega[[k]]%*%matrix(vec_omega_i,ncol=1)
                Y_i[k,] = rmvnorm(n=1, mean = mean_vecY_i_k, sigma = Gamma)
            } else {
                A_state_k = A_mat_scale[,b_i[k]]
                A_1 = diag(A_state_k)
                
                nu_k = D_i[[k]]%*%matrix(vec_alpha_i,ncol=1) + 
                       X_i[[k]]%*%matrix(beta,ncol=1) +
                       D_i_omega[[k]]%*%matrix(vec_omega_i,ncol=1)
                nu_k_1 = D_i[[k-1]]%*%matrix(vec_alpha_i,ncol=1) + 
                         X_i[[k-1]]%*%matrix(beta,ncol=1) +
                         D_i_omega[[k-1]]%*%matrix(vec_omega_i,ncol=1)
                diff_vec = c(Y_i[k-1,] - nu_k_1)
                
                mean_vecY_i_k = nu_k + A_1 %*% matrix(diff_vec,ncol=1)
                
                Y_i[k,] = rmvnorm(n=1, mean = mean_vecY_i_k, sigma = R)
            }
        }
        # ---------------------------------------------------------------------------
        
        t_i = Y[ Y[,'EID']==as.numeric(id_num), 'time', drop=F]
        
        # Removing clinic_rule temporarily
        rules = Y[ Y[,'EID']==as.numeric(id_num), c('RBC_rule'), drop=F]
        rules = cbind(rules, 0)

        if((1 %in% rules[,1]) & !(correct_bleed)) {
            rules[,1] = 0
            print("Bleed rule is changed to 0")
        }
        
        use_data = rbind( use_data, cbind( i, t_i, Y_i, b_i, 
                                           z_i[,2],
                                           x_i[,1], 
                                           rules))
    }
    use_data = matrix(as.numeric(use_data), ncol = ncol(use_data))
    colnames(use_data) = c( 'EID', 'time', 'hemo', 'hr', 'map','lactate', 
                            'b_true', 'RBC_ordered', 
                            'n_RBC_admin', 'RBC_rule', 'clinic_rule')
    
    save( use_data, file=paste0(Dir,'use_data',www,'_', it_num, '.rda'))
    
    cat('\n','Proption of occurances in each state:','\n')
    print(table(use_data[,'b_true'])/dim(use_data)[1])
    cat('\n')
    
    print(paste0("RBC rule is found in ", length(rbc_bleed_correct), " patients"))
    print(paste0(sum(rbc_bleed_correct == 1), " were correct with the bleed event"))
}

bleed_indicator = bleed_indicator_update

save(alpha_i_mat, file = paste0('Data/alpha_i_mat_', it_num, '.rda'))
save(omega_i_mat, file = paste0('Data/omega_i_mat_', it_num, '.rda'))
save(bleed_indicator, file = paste0('Data/bleed_indicator_sim_', it_num, '.rda'))
save(Dn_omega_sim, file = paste0('Data/Dn_omega_sim_', it_num, '.rda'))

# Visualize the noise --------------------------------------------------------
load(paste0('Data/use_data', 1, '_', it_num, '.rda'))

EIDs = unique(use_data[,'EID'])
simulation = T

# New patients ---------------------------------------------------------------
pdf(paste0('Plots/initial_charts', it_num, '.pdf'))
panels = c(4, 1)
par(mfrow=panels, mar=c(2,4,2,4), bg='black', fg='green')
for(i in EIDs){
    # print(which(EIDs == i))
    indices_i = (use_data[,'EID']==i)
    n_i = sum(indices_i)
    t_grid = seq( 0, n_i, by=5)[-1]
    rbc_times = which(use_data[indices_i, 'RBC_ordered'] != 0)
    rbc_admin_times = which(diff(use_data[indices_i, 'n_RBC_admin']) > 0) + 1

    if(simulation){
        b_i = use_data[ indices_i,'b_true']
        to_s1 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==1]
        to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
        to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]
    }
    # HEART RATE --------------------------------------------------------------
    if(sum(!is.na(use_data[indices_i, 'hr']))==0){
        plot.new()
    } else{
        plot(use_data[indices_i, 'hr'], main=paste0('heart rate: ', i, ', RBC Rule = ', mean(use_data[indices_i, 'RBC_rule'])),
             xlab='time', ylab=NA, col.main='green', col.axis='green', pch=20, cex=1)
        grid( nx=20, NULL, col='white')
        axis( side=1, at=t_grid, col.axis='green', labels=t_grid / 4)
        abline(v = rbc_times, col = 'darkorchid2', lwd = 1)
        abline(v = rbc_admin_times, col = 'grey', lwd = 1)
    }
    if(simulation){
        abline( v=to_s1, col='dodgerblue', lwd=2)
        abline( v=to_s2, col='firebrick1', lwd=2)
        abline( v=to_s3, col='yellow2', lwd=2)
        col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
        abline( v= 1, col = col_choice[b_i[1]], lwd = 2)
    }

    # MAP --------------------------------------------------------------
    if(sum(!is.na(use_data[indices_i, 'map']))==0){
        plot.new()
    } else{
        plot(use_data[indices_i, 'map'], main=paste0('map: ', i),
             xlab='time', ylab=NA, col.main='green', col.axis='green', pch=20, cex=1)
        grid( nx=20, NULL, col='white')
        axis( side=1, at=t_grid, col.axis='green', labels=t_grid / 4)
        abline(v = rbc_times, col = 'darkorchid2', lwd = 1)
        abline(v = rbc_admin_times, col = 'grey', lwd = 1)
    }
    if(simulation){
        abline( v=to_s1, col='dodgerblue', lwd=2)
        abline( v=to_s2, col='firebrick1', lwd=2)
        abline( v=to_s3, col='yellow2', lwd=2)
        col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
        abline( v= 1, col = col_choice[b_i[1]], lwd = 2)
    }

    # HEMO --------------------------------------------------------------
    if(sum(!is.na(use_data[indices_i, 'hemo']))==0){
        plot.new()
    } else{
        plot(use_data[indices_i, 'hemo'], main=paste0('hemo: ', i),
             xlab='time', ylab=NA, col.main='green', col.axis='green', pch=20, cex=1)
        grid( nx=20, NULL, col='white')
        axis( side=1, at=t_grid, col.axis='green', labels=t_grid / 4)
        abline(v = rbc_times, col = 'darkorchid2', lwd = 1)
        abline(v = rbc_admin_times, col = 'grey', lwd = 1)
    }
    if(simulation){
        abline( v=to_s1, col='dodgerblue', lwd=2)
        abline( v=to_s2, col='firebrick1', lwd=2)
        abline( v=to_s3, col='yellow2', lwd=2)
        col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
        abline( v= 1, col = col_choice[b_i[1]], lwd = 2)
    }

    # LACTATE --------------------------------------------------------------
    if(sum(!is.na(use_data[indices_i, 'lactate']))==0){
        plot.new()
    } else{
        plot(use_data[indices_i, 'lactate'], main=paste0('lactate: ', i),
             xlab='time', ylab=NA, col.main='green', col.axis='green', pch=20, cex=1)
        grid( nx=20, NULL, col='white')
        axis( side=1, at=t_grid, col.axis='green', labels=t_grid / 4)
        abline(v = rbc_times, col = 'darkorchid2', lwd = 1)
        abline(v = rbc_admin_times, col = 'grey', lwd = 1)
    }
    if(simulation){
        abline( v=to_s1, col='dodgerblue', lwd=2)
        abline( v=to_s2, col='firebrick1', lwd=2)
        abline( v=to_s3, col='yellow2', lwd=2)
        col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
        abline( v= 1, col = col_choice[b_i[1]], lwd = 2)
    }

}
dev.off()