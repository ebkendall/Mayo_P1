# library(MASS, quietly=T)
library(mvtnorm, quietly=T)
library(bayesSurv, quietly=T)
library(expm, quietly=T)

it_num = 5
N = 400

# Load in the existing data and save the covariate combinations
load('Data/data_format_new2.rda')

Y = data_format[, c('EID','time','hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
EIDs = as.character(unique(data_format[,'EID']))

x = data_format[,c('n_RBC_admin'), drop=F]
p = ncol(x)

z = cbind(1, data_format[,c('RBC_ordered'), drop=F])
m = ncol(z)

# Parameters ------------------------------------------------------------------
beta = c(0.6261, -1.3286, 1.6741, -0.1)

# columns: hemo, hr, map
alpha_tilde = matrix( c( 9.57729783, 88.69780576, 79.74903940, 5.2113319,
                                 -1,  5.04150472, -5.42458547,         1,
                                  1, 	      -4, 		    4,        -1), ncol=4, byrow=T)

sigma_upsilon = Upsilon = diag(c(4, 0.25, 0.25, 36, 1, 1, 36, 1, 1, 4, 0.25, 0.25))

A_mat = matrix(c(1.5,  -1, 0,
                 1.5,  -1, 0,
                 1.5,  -1, 0,
                 1.5,  -1, 0), ncol = 3, byrow = T)
vec_A = c(A_mat)
correct_scale_A = exp(vec_A) / (1 + exp(vec_A))
A_mat_scale = matrix(correct_scale_A, nrow = 4)

# columns: hemo, hr, map, lactate
R = diag(4)


# transitions: 1->2, 2->3, 3->1, 3->2
zeta = matrix(c(      -4, -2.578241, -5.000000, -5.230000,
                2.006518,      -1.2,   -1.6713,  1.044297), ncol = 4, byrow=T)


init_logit = c(0,-5,-2)
init_logit = exp(init_logit)

true_pars = c(beta, c(alpha_tilde), c(sigma_upsilon), c(vec_A), c(R), c(zeta), log(init_logit)[2:3])
par_index = list()
par_index$vec_beta = 1:4
par_index$vec_alpha_tilde = 5:16
par_index$vec_sigma_upsilon = 17:160
par_index$vec_A = 161:172
par_index$vec_R = 173:188
par_index$vec_zeta = 189:196
par_index$vec_init = 197:198
save(par_index, file = paste0('Data/true_par_index_', it_num, '.rda'))
save(true_pars, file = paste0('Data/true_pars_', it_num, '.rda'))
# -----------------------------------------------------------------------------

alpha_i_mat = vector(mode = "list", length = N)

for (www in 1:1) {
    set.seed(2023)
    
    Dir = 'Data/'
    
    use_data = NULL
    for(i in 1:N){

        id_num = EIDs[i]
        print(paste0(i, ", ", id_num))
        
        n_i = sum(Y[,'EID']==as.numeric(id_num))
        
        x_i = x[ Y[,'EID']==as.numeric(id_num),, drop=F]
        z_i = z[ Y[,'EID']==as.numeric(id_num),, drop=F]
        
        # Generate realizations of latent bleeding process ---------------------
        D_i = vector(mode = 'list', length = n_i)
        X_i = vector(mode = 'list', length = n_i)
        
        P_i = init_logit / sum(init_logit)
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
        # ---------------------------------------------------------------------------
        
        # Generate realizations of hc, hr, and bp -----------------------------------
        Y_i = matrix(nrow = n_i, ncol = 4)
        vec_alpha_i = rmvnorm( n=1, mean=c(alpha_tilde), sigma=Upsilon)

        alpha_i_mat[[i]] = matrix(vec_alpha_i, ncol = 1)
        
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
                mean_vecY_i_k = D_i[[k]]%*%matrix(vec_alpha_i,ncol=1) + X_i[[k]]%*%matrix(beta,ncol=1)
                Y_i[k,] = rmvnorm(n=1, mean = mean_vecY_i_k, sigma = Gamma)
            } else {
                A_state_k = A_mat_scale[,b_i[k]]
                A_1 = diag(A_state_k)
                
                nu_k = D_i[[k]]%*%matrix(vec_alpha_i,ncol=1) + X_i[[k]]%*%matrix(beta,ncol=1)
                nu_k_1 = D_i[[k-1]]%*%matrix(vec_alpha_i,ncol=1) + X_i[[k-1]]%*%matrix(beta,ncol=1)
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

        if((1 %in% rules[,1]) & !(2 %in% b_i)) {
            rules[,1] = 0
            print("changed rule")
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
}

save(alpha_i_mat, file = paste0('Data/alpha_i_mat_', it_num, '.rda'))

# Visualize the noise --------------------------------------------------------
load(paste0('Data/use_data', 1, '_', it_num, '.rda'))

EIDs = unique(use_data[,'EID'])
simulation = T

# New patients ---------------------------------------------------------------
pdf(paste0('Plots/initial_charts', it_num, '.pdf'))
# mar=c(b,l,t,r) oma=c(b,l,t,r) 
panels = c(4, 1)
par(mfrow=panels, mar=c(2,4,2,4), bg='black', fg='green')
for(i in EIDs){
    print(which(EIDs == i))
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



# # Investigation into each plot -----------------------------------------------
# load('Model_out/final_debug1_it1.rda')
# barplot( rbind( colMeans(final_debug[[2]][[1]][1:4999, 136:163] == 1),
#                 colMeans(final_debug[[2]][[1]][1:4999, 136:163] == 2),
#                 colMeans(final_debug[[2]][[1]][1:4999, 136:163] == 3)), 
#          col=c( 'dodgerblue', 'firebrick1', 'yellow2'), 
#          xlab='time', xaxt='n', space=0, 
#          col.main='green', border='gray') 
# final_debug[[2]][[2]][[1]][,136:162]

# load('Model_out/final_debug1_it1_b.rda')
# barplot( rbind( colMeans(final_debug[[2]][[1]][1:4999, 136:163] == 1),
#                 colMeans(final_debug[[2]][[1]][1:4999, 136:163] == 2),
#                 colMeans(final_debug[[2]][[1]][1:4999, 136:163] == 3)), 
#          col=c( 'dodgerblue', 'firebrick1', 'yellow2'), 
#          xlab='time', xaxt='n', space=0, 
#          col.main='green', border='gray') 
# final_debug[[2]][[2]][[1]][,136:162]
