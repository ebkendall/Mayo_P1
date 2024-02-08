# Size of posterior sample from mcmc chains
n_post = 1001
burnin = 0
steps = 1001
index_post = (steps - burnin - n_post + 1):(steps - burnin)

simulation = F
load("Data/Dn_omega_names3.rda")
load('Data/hr_map_names3.rda')

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
           paste0("mean (hr): ", Dn_omega_names[1:36]), paste0("mean (map): ", Dn_omega_names[37:88]), 
           paste0("log Upsilon (hr): ", Dn_omega_names[1:36]), paste0("log Upsilon (map): ", Dn_omega_names[37:88])) 
additional_labels = c("Gamma(1,1) stable", "Gamma(2,2) stable", "Gamma(3,3) stable", "Gamma(4,4) stable",
                      "Gamma(1,1) bleed", "Gamma(2,2) bleed", "Gamma(3,3) bleed", "Gamma(4,4) bleed",
                      "Gamma(1,1) recov", "Gamma(2,2) recov", "Gamma(3,3) recov", "Gamma(4,4) recov")

dir = 'Model_out/'

if(simulation) {
    index_seeds = c(1:3)
    trialNum = 2
    itNum = 5
    long_chain = F
    
    data_num = 6
    load(paste0('Data/true_pars_', data_num, '.rda'))
    true_par = true_pars     
} else {
    index_seeds = c(1:3)
    trialNum = 8
    itNum = 1
    long_chain = T
    
    true_par = NULL
}

data_format = NULL
if(simulation) {
  load(paste0('Data/use_data1_', data_num, '.rda'))
  data_format = use_data
} else {
  load('Data/data_format_new3.rda')
}

n_subjects = length(unique(data_format[,'EID']))

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
                                trialNum,'it', it, '.rda')
        }

        load(file_name)
        print(paste0(ind, ": ", file_name))
        print(mcmc_out_temp$accept)
        
        par_index = mcmc_out_temp$par_index

        if(it == 1) {
            chain_list[[ind]] = mcmc_out_temp$chain[index_post,]
            for(a in 1:length(a_chain_list[[ind]])) {
                a_chain_list[[ind]][[a]] = mcmc_out_temp$A_chain[[a]]
            }
        } else {
            chain_list[[ind]] = rbind(chain_list[[ind]], mcmc_out_temp$chain[index_post,])
            for(a in 1:length(a_chain_list[[ind]])) {
                a_chain_list[[ind]][[a]] = cbind(a_chain_list[[ind]][[a]], mcmc_out_temp$A_chain[[a]])
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
red_par = matrix(0, nrow=1,ncol=2)
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
                red_par = rbind(red_par, c(r-198, 1))
            }
            if(0 > upper) {
                title_color = "red"
                red_par = rbind(red_par, c(r-198, -1))
            }
        }
        
        y_limit = range(stacked_chains[,r])
        plot( NULL, ylab=NA, main=labels[lab_ind], xlim=c(1,nrow(chain_list[[1]])),
              ylim=y_limit, xlab = paste0("95% CI: [", round(lower, 4),
                                          ", ", round(upper, 4), "]"),
              col.main = title_color)
        
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
        abline( v=lower, col='purple', lwd=2, lty=2)
        abline( v=true_par[r], col='green', lwd=2, lty=2)
    }   
}
red_par = red_par[-1, ]
red_par = cbind(red_par, 0, "hr", 0)
red_par[as.numeric(red_par[,1]) > 36, 4] = "map"
red_par[,5] = Dn_omega_names[as.numeric(red_par[,1])]
red_par = as.data.frame(red_par)
red_par[,1] = as.numeric(red_par[,1])
red_par[,2] = as.numeric(red_par[,2])
red_par[,3] = as.numeric(red_par[,3])
colnames(red_par) = c('ind', 'fit_up_down', 'true_up_down', 'vital', 'name')

if(simulation) {
    print(paste0(sum(true_par[par_index$omega_tilde] != 0), " out of ", length(par_index$omega_tilde), " medications have nonzero effect"))
}

load('Data/med_select_FINAL3.rda')
for(i in 1:nrow(red_par)) {
    if(red_par$vital[i] == "hr") {
        val = as.numeric(unique(med_select_FINAL$hr[med_select_FINAL$med_name_admin == red_par$name[i]]))
        red_par$true_up_down[i] = val
    } else {
        val = as.numeric(unique(med_select_FINAL$map[med_select_FINAL$med_name_admin == red_par$name[i]]))
        red_par$true_up_down[i] = val
    }
}

red_par_diff = red_par[red_par[,'fit_up_down'] != red_par[,'true_up_down'], ,drop = F]
if(nrow(red_par_diff) == 0) {
    print(paste0("All ", nrow(red_par), " match the truth"))
} else {
    print(paste0(nrow(red_par_diff), " out of ", nrow(red_par), " differ"))
    print(red_par_diff)
}

chain_list_gamma = vector(mode = 'list', length = nrow(stacked_chains) / steps)
for(i in 1:length(chain_list_gamma)) {
    max_ind = i * steps
    chain_list_gamma[[i]] = gamma_chain[(max_ind - (steps - 1)):max_ind, ]
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
        base_ind = 3*k - 2
        b_ind = 3*k - 1
        r_ind = 3*k
        xlim_calc = c(min(c(patient_a[b_ind,], patient_a[r_ind,])), 
                      max(c(patient_a[b_ind,], patient_a[r_ind,])))
        
        hist(patient_a[base_ind,], main = paste0(a_chain_id[s], ": ", hist_names[2*k-1]), 
             col = "darkolivegreen4", breaks = floor(sqrt(ncol(patient_a))),
             xlab = "baseline")
        
        hist(patient_a[r_ind,], main = paste0(a_chain_id[s], ": ", hist_names[2*k]), 
             col = "yellow2", breaks = floor(sqrt(ncol(patient_a))), xlim = xlim_calc, 
             xlab = paste0("S2 -> (", 
                           round(mean(patient_a[b_ind,]), 3), ", ", round(sd(patient_a[b_ind,]), 3),
                           "), S3 -> (", 
                           round(mean(patient_a[r_ind,]), 3), ", ", round(sd(patient_a[r_ind,]), 3), ")"))
        
        hist(patient_a[b_ind,], col = "firebrick1", 
             breaks = floor(sqrt(ncol(patient_a))),
             add = T) 
        hist(patient_a[r_ind,], col = "yellow2", 
             breaks = floor(sqrt(ncol(patient_a))), 
             add = T)
        hist(patient_a[b_ind,], col=rgb(1,0,0,0.5), 
             breaks = floor(sqrt(ncol(patient_a))), 
             add=T)   
    }

}

dev.off()
