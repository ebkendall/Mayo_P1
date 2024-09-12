# Size of posterior sample from mcmc chains
steps = 1001

simulation = T

# load("../Data/Dn_omega_names3.rda")
# load('../Data/hr_map_names3.rda')
load('Data_updates/Dn_omega_names1.rda')
load('Data_updates/hr_map_names1.rda')

labels = c("beta (n_RBC_admin): hemo", "beta (n_RBC_admin): hr", 
           "beta (n_RBC_admin): map", "beta (n_RBC_admin): lact",
           "intercept (hemo)", "slope bleeding (hemo)", "slope recovery (hemo)", "slope state 4 (hemo)", "slope state 5 (hemo)",
           "intercept (hr)", "slope bleeding (hr)", "slope recovery (hr)", "slope state 4 (hr)", "slope state 5 (hr)",
           "intercept (map)", "slope bleeding (map)", "slope recovery (map)", "slope state 4 (map)", "slope state 5 (map)",
           "intercept (lact)", "slope bleeding (lact)", "slope recovery (lact)", "slope state 4 (lact)", "slope state 5 (lact)",
           paste0("Upsilon (", 1:20, ", ", rep(1:20, each = 20), ")"), 
           "A1 (baseline)", "A2 (baseline)", "A3 (baseline)", "A4 (baseline)",
           "A1 (bleed)", "A2 (bleed)", "A3 (bleed)", "A4 (bleed)",
           "A1 (recovery)", "A2 (recovery)", "A3 (recovery)", "A4 (recovery)",
           "A1 (state 4)", "A2 (state 4)", "A3 (state 4)", "A4 (state 4)",
           "A1 (state 5)", "A2 (state 5)", "A3 (state 5)", "A4 (state 5)",
           "Var(hemo)", "Cov(hemo, hr)", "Cov(hemo, map)", "Cov(hemo, lact)", 
           "Cov(hr, hemo)", "Var(hr)", "Cov(hr, map)", "Cov(hr, lact)",
           "Cov(map, hemo)", "Cov(map, hr)", "Var(map)", "Cov(map, lact)",
           "Cov(lact, hemo)", "Cov(lact, hr)", "Cov(lact, map)", "Var(lact)",
           "intercept: S1 --> S2", "RBC_order: S1 --> S2",  "intercept: S1 --> S4", "RBC_order: S1 --> S4",
           "intercept: S2 --> S3", "RBC_order: S2 --> S3",  "intercept: S2 --> S4", "RBC_order: S2 --> S4", 
           "intercept: S3 --> S1", "RBC_order: S3 --> S1",  "intercept: S3 --> S2", "RBC_order: S3 --> S2",
           "intercept: S3 --> S4", "RBC_order: S3 --> S4",  "intercept: S4 --> S2", "RBC_order: S4 --> S2",
           "intercept: S4 --> S5", "RBC_order: S4 --> S5",  "intercept: S5 --> S1", "RBC_order: S5 --> S1",
           "intercept: S5 --> S2", "RBC_order: S5 --> S2",  "intercept: S5 --> S4", "RBC_order: S5 --> S4",
           "logit Pr(init S2)", "logit Pr(init S3)","logit Pr(init S4)", "logit Pr(init S5)",
           paste0("mean (hr): ", Dn_omega_names[1:35]), paste0("mean (map): ", Dn_omega_names[36:88]), 
           paste0("log Upsilon (hr): ", Dn_omega_names[1:35]), paste0("log Upsilon (map): ", Dn_omega_names[36:88])) 
additional_labels = c("Gamma(1,1) stable", "Gamma(2,2) stable", "Gamma(3,3) stable", "Gamma(4,4) stable",
                      "Gamma(1,1) bleed", "Gamma(2,2) bleed", "Gamma(3,3) bleed", "Gamma(4,4) bleed",
                      "Gamma(1,1) recov", "Gamma(2,2) recov", "Gamma(3,3) recov", "Gamma(4,4) recov",
                      "Gamma(1,1) state 4", "Gamma(2,2) state 4", "Gamma(3,3) state 4", "Gamma(4,4) state 4",
                      "Gamma(1,1) state 5", "Gamma(2,2) state 5", "Gamma(3,3) state 5", "Gamma(4,4) state 5")


dir = 'Model_out/'

if(simulation) {
    index_seeds = c(1:1)
    trialNum = 1
    itNum = 4
    long_chain = T
    
    data_num = 4
    load(paste0('Data_sim/true_pars_', data_num, '.rda'))
    true_par = true_pars     
} else {
    index_seeds = c(1:3)
    trialNum = 1
    itNum = 1
    long_chain = T
    df_num = 1
    
    true_par = NULL
}

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
chain_list = vector(mode = "list", length = length(index_seeds))
a_chain_list = vector(mode = 'list', length = length(index_seeds))

for(a in 1:length(a_chain_list)) {
    a_chain_list[[a]] = vector(mode = 'list', length = 10)
}

ind = 0

for(seed in index_seeds){
    ind = ind + 1
    if(long_chain) {
        it_seq = 1:itNum
    } else {
        it_seq = itNum
    }
    
    for(it in it_seq) {
        if(simulation) {
            file_name = paste0(dir,'mcmc_out_interm_',toString(seed),'_', 
                                trialNum,'it', it, '_sim.rda') 
        } else {
            file_name = paste0(dir,'mcmc_out_interm_',toString(seed),'_', 
                               trialNum,'it', it, '_df', df_num, '.rda')
        }

        load(file_name)
        print(paste0(ind, ": ", file_name))
        print("accept")
        print(mcmc_out_temp$accept)
        print("pscale")
        print(mcmc_out_temp$pscale)

        par_index = mcmc_out_temp$par_index

        if(it == 1) {
            chain_list[[ind]] = mcmc_out_temp$chain
            for(a in 1:length(a_chain_list[[ind]])) {
                a_chain_list[[ind]][[a]] = mcmc_out_temp$A_chain[[a]]
            }
        } else {
            chain_list[[ind]] = rbind(chain_list[[ind]], mcmc_out_temp$chain)
            for(a in 1:length(a_chain_list[[ind]])) {
                a_chain_list[[ind]][[a]] = cbind(a_chain_list[[ind]][[a]], mcmc_out_temp$A_chain[[a]])
            }
        }
        
        rm(mcmc_out_temp)
    }
}

stacked_chains = do.call( rbind, chain_list)

# Re-calculating the Upsilon matrix
true_gamma = NULL

gamma_chain = matrix(nrow = nrow(stacked_chains), ncol = 20)
for(i in 1:nrow(stacked_chains)) {
    R = matrix(stacked_chains[i, par_index$vec_R], ncol = 4)
    vec_A1 = stacked_chains[i, par_index$vec_A]
    scale_A1 = exp(vec_A1) / (1+exp(vec_A1))
    
    diag_gamma = c(R[1,1] / (scale_A1[1]^2), R[2,2] / (scale_A1[2]^2),
                   R[3,3] / (scale_A1[3]^2), R[4,4] / (scale_A1[4]^2),
                   R[1,1] / (scale_A1[5]^2), R[2,2] / (scale_A1[6]^2),
                   R[3,3] / (scale_A1[7]^2), R[4,4] / (scale_A1[8]^2),
                   R[1,1] / (scale_A1[9]^2), R[2,2] / (scale_A1[10]^2),
                   R[3,3] / (scale_A1[11]^2), R[4,4] / (scale_A1[12]^2),
                   R[1,1] / (scale_A1[13]^2), R[2,2] / (scale_A1[14]^2),
                   R[3,3] / (scale_A1[15]^2), R[4,4] / (scale_A1[16]^2),
                   R[1,1] / (scale_A1[17]^2), R[2,2] / (scale_A1[18]^2),
                   R[3,3] / (scale_A1[19]^2), R[4,4] / (scale_A1[20]^2))

    gamma_chain[i, ] = diag_gamma
}

pdf_title = NULL
if(simulation) {
    pdf_title = paste0('Plots/trace_plot_', trialNum, '_it', itNum, '_sim.pdf')
} else {
    pdf_title = paste0('Plots/trace_plot_', trialNum, '_it', itNum, 
                       '_df', df_num, '.pdf')
}
pdf(pdf_title)
par(mfrow=c(3, 2))
lab_ind = 0
for(s in names(par_index)){
    temp_par = par_index[[s]]
    for(r in temp_par){
        # lab_ind = lab_ind + 1
        lab_ind = r
        parMean = round( mean(stacked_chains[,r]), 4)
        parMedian = round( median(stacked_chains[,r]), 4)
        upper = quantile( stacked_chains[,r], prob=.975)
        lower = quantile( stacked_chains[,r], prob=.025)
        
        title_color = "black"
        if(s == names(par_index)[8]) {
            if(0 < lower) {
                title_color = "red"
            }
            if(0 > upper) {
                title_color = "red"
            }
        }
        
        y_limit = range(stacked_chains[,r])
        plot( NULL, ylab=NA, main=labels[lab_ind], xlim=c(1,nrow(chain_list[[1]])),
              ylim=y_limit, xlab = paste0("95% CI: [", round(lower, 4),
                                          ", ", round(upper, 4), "]"),
              col.main = title_color)
        
        for(seed in 1:length(chain_list)) {
            lines( chain_list[[seed]][,r], type='l', col=seed)
            abline(h = chain_list[[seed]][1,r], col=seed)
        }
        
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
        abline( v=lower, col='purple', lwd=2, lty=2)
        abline( v=true_par[r], col='green', lwd=2, lty=2)
    }   
}

if(long_chain) {
    chain_list_gamma = vector(mode = 'list', length = nrow(stacked_chains) / (itNum*steps))
} else {
    chain_list_gamma = vector(mode = 'list', length = nrow(stacked_chains) / steps)   
}
for(i in 1:length(chain_list_gamma)) {
    if(long_chain) {
        max_ind = i * (itNum*steps)
        chain_list_gamma[[i]] = gamma_chain[(max_ind - (itNum*steps - 1)):max_ind, ]
    } else {
        max_ind = i * steps   
        chain_list_gamma[[i]] = gamma_chain[(max_ind - (steps - 1)):max_ind, ]
    }
}

for(rr in 1:ncol(gamma_chain)){

    lab_ind = rr
    parMean = round( mean(gamma_chain[,rr]), 4)
    parMedian = round( median(gamma_chain[,rr]), 4)
    upper = quantile( gamma_chain[,rr], prob=.975)
    lower = quantile( gamma_chain[,rr], prob=.025)

    y_limit = range(gamma_chain[,rr])

    plot( NULL, ylab=NA, main=additional_labels[lab_ind], xlim=c(1,nrow(chain_list[[1]])),
          ylim=y_limit, xlab = paste0("95% CI: [", round(lower, 4),
                                      ", ", round(upper, 4), "]"))

    for(seed in 1:length(chain_list_gamma)) {
        lines( chain_list_gamma[[seed]][,rr], type='l', col=seed)
        abline(h = chain_list_gamma[[seed]][1,rr], col=seed)
    }

    x_label = paste0('Mean =',toString(parMean),' Median =',toString(parMedian))
    
    hist( gamma_chain[,rr], breaks=sqrt(nrow(gamma_chain)), ylab=NA, main=NA, freq=FALSE,
          xlab=x_label)
    abline( v=upper, col='red', lwd=2, lty=2)
    abline( v=lower, col='purple', lwd=2, lty=2)
    abline( v=true_gamma[rr], col='green', lwd=2, lty=2)
}

# Plotting the sampled alpha_i
a_chain_id = c(3, 86, 163, 237, 427, 521, 632, 646, 692, 713)
hist_a_chain_list = vector(mode = 'list', length = length(a_chain_id))
for(i in 1:length(a_chain_id)) {
    for(j in 1:length(index_seeds)) {
        if(j==1) {
            hist_a_chain_list[[i]] = a_chain_list[[j]][[i]]
        } else {
            hist_a_chain_list[[i]] = cbind(hist_a_chain_list[[i]], a_chain_list[[j]][[i]])
        }
    }
}

hist_names = c("alpha_i baseline for hemo", "alpha_i slopes for hemo",
               "alpha_i baseline for hr", "alpha_i slopes for hr",
               "alpha_i baseline for map", "alpha_i slopes for map",
               "alpha_i baseline for lactate", "alpha_i slopes for lactate")

for (s in 1:length(a_chain_id)) {

    patient_a = hist_a_chain_list[[s]]
    
    for(k in 1:4) {
        base_ind = 5*k - 4
        b_ind = 5*k - 3
        r_ind = 5*k - 2
        s4_ind = 5*k - 1
        s5_ind = 5*k
        
        temp_df = data.frame("bleed" = patient_a[b_ind,],
                             "recovery" = patient_a[r_ind,],
                             "NBE" = patient_a[s4_ind, ],
                             "recov_NBE" = patient_a[s5_ind, ])
        
        hist(patient_a[base_ind,], main = paste0(a_chain_id[s], ": ", hist_names[2*k-1]), 
             col = "darkolivegreen4", breaks = floor(sqrt(ncol(patient_a))),
             xlab = "baseline")
        
        boxplot(temp_df, col = c('firebrick1', 'yellow2', 'green', 'darkgray'),
                main = hist_names[2*k], outline=FALSE)
        abline(h = 0, col = 'blue')
    }

}

dev.off()
