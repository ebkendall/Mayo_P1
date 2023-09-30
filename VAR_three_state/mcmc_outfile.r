library(matrixStats)
library(mvtnorm)
library(MASS)

dir = 'Model_out/' # Change this everytime!!!! ****************

# Size of posterior sample from mcmc chains
n_post = 1000
# Step number at 3ich the adaptive tuning scheme was frozen
burnin = 0
# Total number of steps the mcmc algorithm is computed for
steps = 1000
# Matrix row indices for the posterior sample to use for GFF computation
index_post = (steps - burnin - n_post + 1):(steps - burnin)
# par_index$vec_beta = 1:4
# par_index$vec_alpha_tilde = 5:16
# par_index$vec_sigma_upsilon = 17:160
# par_index$vec_logit_A = 161:164
# par_index$vec_R = 165:168
# par_index$vec_zeta = 169:176
# par_index$vec_init = 177:178
# par_index$log_lambda = 179:190
# par_index$omega_tilde = 191:198
# par_index$vec_upsilon_omega = 199:262

simulation = T
data_num = 4

labels = c("beta (n_RBC_admin): hemo", "beta (n_RBC_admin): hr", 
           "beta (n_RBC_admin): map", "beta (n_RBC_admin): lact",
           "intercept (hemo)", "slope bleeding (hemo)", "slope recovery (hemo)",
           "intercept (hr)", "slope bleeding (hr)", "slope recovery (hr)",
           "intercept (map)", "slope bleeding (map)", "slope recovery (map)",
           "intercept (lact)", "slope bleeding (lact)", "slope recovery (lact)",
           paste0("Upsilon (", 1:12, ", ", rep(1:12, each = 12), ")"), 
           "A1 (baseline)", "A2 (baseline)", "A3 (baseline)", "A4 (baseline)",
           "A1 (bleed)", "A2 (bleed)", "A3 (bleed)", "A4 (bleed)",
           "A1 (recovery)", "A2 (recovery)", "A3 (recovery)", "A4 (recovery)",
           "Var(hemo)", "Cov(hemo, hr)", "Cov(hemo, map)", "Cov(hemo, lact)", 
           "Cov(hr, hemo)", "Var(hr)", "Cov(hr, map)", "Cov(hr, lact)",
           "Cov(map, hemo)", "Cov(map, hr)", "Var(map)", "Cov(map, lact)",
           "Cov(lact, hemo)", "Cov(lact, hr)", "Cov(lact, map)", "Var(lact)",
           "intercept: S1 --> S2", "RBC_order: S1 --> S2",  "intercept: S2 --> S3", "RBC_order: S2 --> S3", 
           "intercept: S3 --> S1", "RBC_order: S3 --> S1",  "intercept: S3 --> S2", "RBC_order: S3 --> S2",
           "logit Pr(init S2)", "logit Pr(init S3)",
           "log(lambda): intercept (hemo)", "log(lambda): slope bleeding (hemo)", "log(lambda): slope recovery (hemo)",
           "log(lambda): intercept (hr)", "log(lambda): slope bleeding (hr)", "log(lambda): slope recovery (hr)",
           "log(lambda): intercept (map)", "log(lambda): slope bleeding (map)", "log(lambda): slope recovery (map)",
           "log(lambda): intercept (lact)", "log(lambda): slope bleeding (lact)", "log(lambda): slope recovery (lact)"
        #    "omega_tilde(1,1)", "omega_tilde(1,2)", "omega_tilde(2,1)", "omega_tilde(2,2)",
        #    "omega_tilde(3,1)", "omega_tilde(3,2)", "omega_tilde(4,1)", "omega_tilde(5,2)",
        #    paste0("Upsilon_omega (", 1:8, ", ", rep(1:8, each = 8), ")")
           ) 
additional_labels = c("Gamma(1,1) stable", "Gamma(2,2) stable", "Gamma(3,3) stable", "Gamma(4,4) stable",
                      "Gamma(1,1) bleed", "Gamma(2,2) bleed", "Gamma(3,3) bleed", "Gamma(4,4) bleed",
                      "Gamma(1,1) recov", "Gamma(2,2) recov", "Gamma(3,3) recov", "Gamma(4,4) recov")

if(simulation) {
    index_seeds = c(1:5)
    trialNum = 3
    itNum = 3
} else {
    index_seeds = c(1:2,4:5)
    trialNum = 6 # Change this everytime!!!! ****************
    itNum = 5
}
# load('Model_out/mcmc_out_interm_3_13it10.rda')
# par_temp = colMeans(mcmc_out_temp$chain)
# rownames(par_temp) = NULL
if(simulation) {
    load(paste0('Data/true_pars_', data_num, '.rda'))
    load(paste0('Data/true_par_index_', data_num, '.rda'))
    # true_pars[par_index$vec_sigma_upsilon] = diag(exp(true_pars[par_index$log_lambda])) %*% 
    #                                             matrix(true_pars[par_index$vec_sigma_upsilon], ncol=12) %*% 
    #                                                 diag(exp(true_pars[par_index$log_lambda])) 
    true_par = true_pars                                                
} else {
    true_par = NULL
}

par_index = NULL
accept_rat = rep(NA, length(index_seeds))

data_format = NULL
if(simulation) {
  load(paste0('Data/use_data1_', data_num, '.rda'))
  data_format = use_data
} else {
  load('Data/data_format_new2.rda')
#   pace_id = c(53475, 110750, 125025, 260625, 273425, 296500, 310100, 384925,
#               417300, 448075, 538075, 616025, 660075, 665850, 666750, 677225,
#               732525, 758025, 763050, 843000, 117525)
#   data_format = data_format[!(data_format[,'EID'] %in% pace_id), ]
}

n_subjects = length(unique(data_format[,'EID']))

chain_base = chain_bleed = chain_recov = 
    chain_NBE = chain_l_recov = vector(mode = 'list', length = 4)

chain_base[[1]] = chain_base[[2]] = chain_base[[3]] = chain_base[[4]] =     
    chain_bleed[[1]] = chain_bleed[[2]] = chain_bleed[[3]] = chain_bleed[[4]] =
        chain_recov[[1]] = chain_recov[[2]] = chain_recov[[3]] = chain_recov[[4]] =
                    rep(NA, n_subjects * length(index_seeds))

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
chain_list = vector(mode = "list", length = length(index_seeds))

ind = 0

for(seed in index_seeds){
    if(simulation) {
        file_name = paste0(dir,'mcmc_out_interm_',toString(seed),'_', trialNum,'it', itNum, '_sim.rda') 
    } else {
        file_name = paste0(dir,'mcmc_out_interm_',toString(seed),'_', trialNum,'it', itNum, '.rda')
    }
    if (file.exists(file_name)) {
        load(file_name)
        ind = ind + 1
        print(paste0(ind, ": ", file_name))
        print(mcmc_out_temp$accept)
        
        par_index = mcmc_out_temp$par_index

        chain_list[[ind]] = mcmc_out_temp$chain[index_post,]

        # Adding the A_chain components
        for (s in 1:4) {
            print(s)
            main_ind = min(which(is.na(chain_bleed[[s]])))
            for(i in 1:n_subjects) {
                chain_base[[s]][main_ind] = mcmc_out_temp$A_chain[[5]][[i]][3*s-2]
                chain_bleed[[s]][main_ind] = mcmc_out_temp$A_chain[[5]][[i]][3*s-1]
                chain_recov[[s]][main_ind] = mcmc_out_temp$A_chain[[5]][[i]][3*s]
                main_ind = main_ind + 1
            }
        }

        rm(mcmc_out_temp)
    }
}

stacked_chains = do.call( rbind, chain_list)

# # Re-calculating the Upsilon matrix
if(simulation) {
    true_R = matrix(true_par[par_index$vec_R], ncol = 4)
    true_A = true_par[par_index$vec_A]
    true_A_scale = exp(true_A) / (1+exp(true_A))
    true_gamma = c(true_R[1,1] / (true_A_scale[1]^2), true_R[2,2] / (true_A_scale[2]^2),
                   true_R[3,3] / (true_A_scale[3]^2), true_R[4,4] / (true_A_scale[4]^2),
                   true_R[1,1] / (true_A_scale[5]^2), true_R[2,2] / (true_A_scale[6]^2),
                   true_R[3,3] / (true_A_scale[7]^2), true_R[4,4] / (true_A_scale[8]^2),
                   true_R[1,1] / (true_A_scale[9]^2), true_R[2,2] / (true_A_scale[10]^2),
                   true_R[3,3] / (true_A_scale[11]^2), true_R[4,4] / (true_A_scale[12]^2))
} else {
    true_gamma = NULL
}

gamma_chain = matrix(nrow = nrow(stacked_chains), ncol = 12)
for(i in 1:nrow(stacked_chains)) {
    R = matrix(stacked_chains[i, par_index$vec_R], ncol = 4)
    vec_A1 = stacked_chains[i, par_index$vec_A]
    scale_A1 = exp(vec_A1) / (1+exp(vec_A1))
    
    diag_gamma = c(R[1,1] / (scale_A1[1]^2), R[2,2] / (scale_A1[2]^2),
                   R[3,3] / (scale_A1[3]^2), R[4,4] / (scale_A1[4]^2),
                   R[1,1] / (scale_A1[5]^2), R[2,2] / (scale_A1[6]^2),
                   R[3,3] / (scale_A1[7]^2), R[4,4] / (scale_A1[8]^2),
                   R[1,1] / (scale_A1[9]^2), R[2,2] / (scale_A1[10]^2),
                   R[3,3] / (scale_A1[11]^2), R[4,4] / (scale_A1[12]^2))

    gamma_chain[i, ] = diag_gamma
}

pdf_title = NULL
if(simulation) {
    pdf_title = paste0('Plots/trace_plot_', trialNum, '_it', itNum, '_sim.pdf')
} else {
    pdf_title = paste0('Plots/trace_plot_', trialNum, '_it', itNum, '.pdf')
}
pdf(pdf_title)
par(mfrow=c(3, 2))
lab_ind = 0
for(s in names(par_index)){
    if(which(s == names(par_index)) < 9) {
        temp_par = par_index[[s]]
        if (s == names(par_index)[3]) {
            temp_par = temp_par[c(1 ,  14,  27,  40,  53, 66, 79,
                                  92, 105, 118, 131, 144)]
        }
        
        for(r in temp_par){
            # lab_ind = lab_ind + 1
            lab_ind = r
            parMean = round( mean(stacked_chains[,r]), 4)
            parMedian = round( median(stacked_chains[,r]), 4)
            upper = quantile( stacked_chains[,r], prob=.975)
            lower = quantile( stacked_chains[,r], prob=.025)
            
            y_limit = range(stacked_chains[,r])
            plot( NULL, ylab=NA, main=labels[lab_ind], xlim=c(1,length(index_post)),
                  ylim=y_limit, xlab = paste0("95% CI: [", round(lower, 4),
                                              ", ", round(upper, 4), "]"))
            
            for(seed in 1:length(chain_list)) lines( chain_list[[seed]][,r], type='l', col=seed)
            
            if (simulation) {
                x_label = paste0('Mean =',toString(parMean),
                                 ' Median =',toString(parMedian),
                                 ' True =', round(true_par[r], 3))
            } else {
                x_label = paste0('Mean =',toString(parMean),' Median =',toString(parMedian))
            }
            hist( stacked_chains[,r], breaks=sqrt(nrow(stacked_chains)), ylab=NA, main=NA, freq=FALSE,
                  xlab=x_label)
            abline( v=upper, col='red', lwd=2, lty=2)
            abline( v=true_par[r], col='green', lwd=2, lty=2)
            abline( v=lower, col='purple', lwd=2, lty=2)
        }   
    }
}

chain_list_gamma = vector(mode = 'list', length = nrow(stacked_chains) / 1000)
for(i in 1:length(chain_list_gamma)) {
    max_ind = i * 1000
    chain_list_gamma[[i]] = gamma_chain[(max_ind - 999):max_ind, ]
}

for(rr in 1:ncol(gamma_chain)){
    # lab_ind = lab_ind + 1
    lab_ind = rr
    parMean = round( mean(gamma_chain[,rr]), 4)
    parMedian = round( median(gamma_chain[,rr]), 4)
    upper = quantile( gamma_chain[,rr], prob=.975)
    lower = quantile( gamma_chain[,rr], prob=.025)
    
    y_limit = range(gamma_chain[,rr])
    
    plot( NULL, ylab=NA, main=additional_labels[lab_ind], xlim=c(1,length(index_post)),
          ylim=y_limit, xlab = paste0("95% CI: [", round(lower, 4),
                                      ", ", round(upper, 4), "]"))
    
    for(seed in 1:length(chain_list_gamma)) lines( chain_list_gamma[[seed]][,rr], type='l', col=seed)
    
    if (simulation) {
        x_label = paste0('Mean =',toString(parMean),
                         ' Median =',toString(parMedian),
                         ' True =', round(true_gamma[rr], 3))
    } else {
        x_label = paste0('Mean =',toString(parMean),' Median =',toString(parMedian))
    }
    hist( gamma_chain[,rr], breaks=sqrt(nrow(gamma_chain)), ylab=NA, main=NA, freq=FALSE,
          xlab=x_label)
    abline( v=upper, col='red', lwd=2, lty=2)
    abline( v=true_gamma[rr], col='green', lwd=2, lty=2)
    abline( v=lower, col='purple', lwd=2, lty=2)
}

hist_names = c("alpha_i slopes for hemo",
               "alpha_i slopes for hr",
               "alpha_i slopes for map",
               "alpha_i slopes for lactate")
base_names = c("alpha_i baseline for hemo",
               "alpha_i baseline for hr",
               "alpha_i baseline for map",
               "alpha_i baseline for lactate")

for (s in 1:4) {

    xlim_calc = c(min(c(chain_bleed[[s]], chain_recov[[s]])), max(c(chain_bleed[[s]], chain_recov[[s]])))

    hist(chain_base[[s]], main = base_names[s], col = "darkolivegreen4", breaks = sqrt(length(chain_base[[s]])),
         xlab = "baseline")

    hist(chain_recov[[s]], main = hist_names[s], col = "yellow2", breaks = sqrt(length(chain_recov[[s]])),
            xlim = xlim_calc, 
            xlab = paste0("S2 -> (", round(mean(chain_bleed[[s]]), 3), ", ", round(sd(chain_bleed[[s]]), 3),
                          "), S3 -> (", round(mean(chain_recov[[s]]), 3), ", ", round(sd(chain_recov[[s]]), 3), ")"))
    hist(chain_bleed[[s]], col = "firebrick1", breaks = sqrt(length(chain_bleed[[s]])),
         add = T) 
    hist(chain_recov[[s]], col = "yellow2", breaks = sqrt(length(chain_recov[[s]])), add = T)
    hist(chain_bleed[[s]], col=rgb(1,0,0,0.5), breaks = sqrt(length(chain_bleed[[s]])), add=T)
    
    temp_df = data.frame("bleed" = chain_bleed[[s]],
                         "recov_bleed" = chain_recov[[s]])
                        #  "NBE" = chain_NBE[[s]],
                        #  "recov_NBE" = chain_l_recov[[s]]
    # boxplot(temp_df, col = c('firebrick1', 'yellow2', 'green', 'darkgray'),
    #         main = hist_names[s])

    print(base_names[s])
    print(summary(temp_df$bleed)); print(summary(temp_df$recov_bleed))

}

dev.off()

# Investigation into each plot -----------------------------------------------
# load('Model_out/final_debug1_6it1i.rda')
# barplot( rbind( colMeans(final_debug[["l_1"]][[1]][1:4999, 1:187] == 1),
#                 colMeans(final_debug[["l_1"]][[1]][1:4999, 1:187] == 2),
#                 colMeans(final_debug[["l_1"]][[1]][1:4999, 1:187] == 3)),
#          col=c( 'dodgerblue', 'firebrick1', 'yellow2'),
#          xlab='time', xaxt='n', space=0,
#          col.main='green', border='gray')
# final_debug[["l_1"]][[2]][[1]][,136:162]

# load('Model_out/final_debug1_6it1g.rda')
# barplot( rbind( colMeans(final_debug[["l_2"]][[1]][1:4999, 1:117] == 1),
#                 colMeans(final_debug[["l_2"]][[1]][1:4999, 1:117] == 2),
#                 colMeans(final_debug[["l_2"]][[1]][1:4999, 1:117] == 3)),
#          col=c( 'dodgerblue', 'firebrick1', 'yellow2'),
#          xlab='time', xaxt='n', space=0,
#          col.main='green', border='gray')
# final_debug[["l_2"]][[2]][[1]][,90:116]
