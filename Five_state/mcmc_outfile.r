library(matrixStats)
library(mvtnorm)
library(MASS)

dir = 'Model_out/' # Change this everytime!!!! ****************

# Size of posterior sample from mcmc chains
n_post = 700
# Step number at 3ich the adaptive tuning scheme was frozen
burnin = 0
# Total number of steps the mcmc algorithm is computed for
steps = 1000
# Matrix row indices for the posterior sample to use for GFF computation
index_post = (steps - burnin - n_post + 1):(steps - burnin)

labels = c("beta (n_RBC_admin): hemo", "beta (n_RBC_admin): hr", 
           "beta (n_RBC_admin): map", "beta (n_RBC_admin): lact",
           "intercept (hemo)", "slope bleeding (hemo)", "slope recovery (hemo)",
           "slope NBE (hemo)", "slope NBE recov (hemo)",
           "intercept (hr)", "slope bleeding (hr)", "slope recovery (hr)",
           "slope NBE (hr)", "slope NBE recov (hr)",
           "intercept (map)", "slope bleeding (map)", "slope recovery (map)",
           "slope NBE (map)", "slope NBE recov (map)",
           "intercept (lact)", "slope bleeding (lact)", "slope recovery (lact)",
           "slope NBE (lact)", "slope NBE recov (lact)",
           paste0("Upsilon (", 1:20, ", ", rep(1:20, each = 20), ")"), "log theta",
           "Var(hemo)", "Cov(hemo, hr)", "Cov(hemo, map)", "Cov(hemo, lact)", 
           "Cov(hr, hemo)", "Var(hr)", "Cov(hr, map)", "Cov(hr, lact)",
           "Cov(map, hemo)", "Cov(map, hr)", "Var(map)", "Cov(map, lact)",
           "Cov(lact, hemo)", "Cov(lact, hr)", "Cov(lact, map)", "Var(lact)",
           "intercept: S1 --> S2", "RBC_order: S1 --> S2", "intercept: S1 --> S4", "RBC_order: S1 --> S4",
           "intercept: S2 --> S3", "RBC_order: S2 --> S3", "intercept: S2 --> S4", "RBC_order: S2 --> S4",
           "intercept: S3 --> S1", "RBC_order: S3 --> S1", "intercept: S3 --> S2", "RBC_order: S3 --> S2",
           "intercept: S3 --> S4", "RBC_order: S3 --> S4", "intercept: S4 --> S2", "RBC_order: S4 --> S2",
           "intercept: S4 --> S5", "RBC_order: S4 --> S5", "intercept: S5 --> S1", "RBC_order: S5 --> S1",
           "intercept: S5 --> S2", "RBC_order: S5 --> S2", "intercept: S5 --> S4", "RBC_order: S5 --> S4",
           "logit Pr(init S2)", "logit Pr(init S3)", "logit Pr(init S4)", 
           "logit Pr(init S5)", 
           "log(lambda): intercept (hemo)", "log(lambda): slope bleeding (hemo)", "log(lambda): slope recovery (hemo)",
           "log(lambda): slope NBE (hemo)", "log(lambda): slope NBE recov (hemo)",
           "log(lambda): intercept (hr)", "log(lambda): slope bleeding (hr)", "log(lambda): slope recovery (hr)",
           "log(lambda): slope NBE (hr)", "log(lambda): slope NBE recov (hr)",
           "log(lambda): intercept (map)", "log(lambda): slope bleeding (map)", "log(lambda): slope recovery (map)",
           "log(lambda): slope NBE (map)", "log(lambda): slope NBE recov (map)",
           "log(lambda): intercept (lact)", "log(lambda): slope bleeding (lact)", "log(lambda): slope recovery (lact)",
           "log(lambda): slope NBE (lact)", "log(lambda): slope NBE recov (lact)") 


index_seeds = c(1:3)
trialNum = 1 # Change this everytime!!!! ****************
itNum = 4

# load('true_par1.rda')
# load('../_Developement/Model_out/mcmc_out_interm_4_11it10.rda')
# true_par = colMeans(mcmc_out_temp$chain)
# load('../Simulation/Data/Debug/true_pars_11.rda')
# load('../Simulation/Data/Debug/true_par_index_11.rda')
# true_par = true_pars
true_par = NULL

# Sigma = matrix(true_par[par_index$vec_sigma_upsilon], ncol = 12)
# Lambda = diag(exp(true_par[par_index$log_lambda]))
# Upsilon = Lambda %*% Sigma %*% Lambda
# true_par[par_index$vec_sigma_upsilon] = c(Upsilon)

par_index = NULL
accept_rat = rep(NA, length(index_seeds))

n_subjects = 177
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
    file_name = paste0(dir,'mcmc_out_interm_',toString(seed),'_', trialNum,'it', itNum, '.rda')
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

# Re-calculating the Upsilon matrix
for(i in 1:nrow(stacked_chains)) {
    Sigma = matrix(stacked_chains[i,par_index$vec_sigma_upsilon], ncol = 20)
    Lambda = diag(exp(stacked_chains[i,par_index$log_lambda]))

    Upsilon = Lambda %*% Sigma %*% Lambda
    stacked_chains[i, par_index$vec_sigma_upsilon] = c(Upsilon)
}

pdf(paste0('Plots/trace_plot_', trialNum, '_it', itNum, '.pdf'))
par(mfrow=c(3, 2))
lab_ind = 0
for(s in names(par_index)){

	temp_par = par_index[[s]]
    if (s == names(par_index)[3]) {
        temp_par = temp_par[c(1, 22, 43, 64, 85, 106, 127, 148, 169, 190, 211, 
                              232, 253, 274, 295, 316, 337, 358, 379, 400)]
    }

	for(r in temp_par){
        # lab_ind = lab_ind + 1
        lab_ind = r
        # stacked_chains[,r]
		plot( NULL, ylab=NA, main=labels[lab_ind], xlim=c(1,length(index_post)),
              ylim=range(stacked_chains[,r]), xlab = "")
              # xlab=paste0("true val: ", round(true_par[r], 4)))
            #   paste0('Accept Rat: ',toString(round( accept_rat, 4)))

        for(seed in 1:length(chain_list)) lines( chain_list[[seed]][,r], type='l', col=seed)

		parMean = round( mean(stacked_chains[,r]), 4)
		parMedian = round( median(stacked_chains[,r]), 4)
		upper = quantile( stacked_chains[,r], prob=.975)
		lower = quantile( stacked_chains[,r], prob=.025)

		hist( stacked_chains[,r], breaks=sqrt(nrow(stacked_chains)), ylab=NA, main=NA, freq=FALSE,
			  xlab=paste0('Mean =',toString(parMean),' Median =',toString(parMedian)))
		abline( v=upper, col='red', lwd=2, lty=2)
		# abline( v=true_par[r], col='green', lwd=2, lty=2)
		abline( v=lower, col='purple', lwd=2, lty=2)
	}
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

