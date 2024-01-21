# Model evaluation to determine bleeding threshold and optimal time window.

# -----------------------------------------------------------------------------
# Cumulative posterior probability plot # -------------------------------------
# -----------------------------------------------------------------------------
indices_i_new = which(indices_i == T)
cumulative_post_prob = matrix(nrow = 2, ncol = length(indices_i_new))
for(w in 1:length(indices_i_new)) {
	start_index = indices_i_new[1]
	end_index = indices_i_new[w] 
	if(w - 24 > 0) start_index = indices_i_new[w - 24]
	y_or_n_2 = apply(mcmc_out_temp$B_chain[, start_index:end_index, drop=F],
						1, function(x) (2 %in% x))
	prob_2 = mean(y_or_n_2)

	cumulative_post_prob[,w] = c(prob_2, 1-prob_2)
}
barplot( cumulative_post_prob, 
				col=c('firebrick1', 'black'), 
				main=paste0('cumulative prob.'), xlab='time', xaxt='n', space=0, 
				col.main='green', border=NA) 
grid( nx=NA, NULL, col='white')
legend( 'topright', inset=c(0,-.28), xpd=T, horiz=T, bty='n', x.intersp=.75,
		legend=c( 'bleeding', 'LIMBO'), pch=15, pt.cex=1.5, 
				col=c('firebrick1', 'black'))
axis( side=1, at=t_grid, col.axis='green', labels=t_grid / 4)
axis( side=2, at=0:1, col.axis='green')
abline(v = rbc_times, col = 'darkorchid2', lwd = 2)
if(simulation){
	abline( v=to_s1, col='cyan', lwd=3)
	abline( v=to_s2, col='brown4', lwd=3)
	abline( v=to_s3, col='darkgoldenrod3', lwd=3)
	col_choice = c('cyan', 'brown4', 'darkgoldenrod3')
	abline( v= 1, col = col_choice[b_i[1]], lwd = 2)
}



# -----------------------------------------------------------------------------
# ROC curve # -----------------------------------------------------------------
# -----------------------------------------------------------------------------
# Construct the ROC curve and determine what is a 'reasonable' threshold to say
# a patient is bleeding.

library(splines)

# Each run of the simulation has a different seed associated with it
index_seeds = c(1,2,3,6,7,8)
trialNum = 14
window_length = 1:12
optimal_c_win = vector(mode = "list", length = length(window_length))

for(www in 1:length(window_length)) {
    best_c = data.frame("c" = rep(NA, length(index_seeds)),
                        "true_pos" = rep(NA, length(index_seeds)),
                        "false_pos" = rep(NA, length(index_seeds)))
    best_c_ind = 1

    win_length = window_length[www]
    print(paste0(win_length, " (",www, ")"))

    for(seed_ind in index_seeds) {
        
        print(seed_ind)
        load(paste0('Data/use_data_small', seed_ind, '.rda'))
        load(paste0('Model_out/post_mcmc_out_dev', seed_ind, '_', trialNum, '.rda'))
        mcmc_out$B_chain = mcmc_out$B_chain[10000:15000, ]
        mcmc_out$hr_chain = mcmc_out$hr_chain[10000:15000, ]
        mcmc_out$bp_chain = mcmc_out$bp_chain[10000:15000, ]
        mcmc_out$hc_chain = mcmc_out$hc_chain[10000:15000, ]

        # Truth
        true_state = use_data[,'b_true']

        # Model output
        cumulative_post_prob = rep(NA, length(true_state))
        ind = 1
        for(i in unique(use_data[,'EID'])) {
            indices_i = (use_data[,'EID']==i)
            indices_i_new = which(indices_i == T)
            for(w in 1:length(indices_i_new)) {
                start_index = indices_i_new[1]
                end_index = indices_i_new[w] 
                if(w - win_length > 0) start_index = indices_i_new[w - win_length]

                y_or_n_2 = apply(mcmc_out$B_chain[, start_index:end_index, drop=F],
                                    1, function(x) (2 %in% x))
                prob_2 = mean(y_or_n_2)
                
                cumulative_post_prob[ind] = prob_2
                ind = ind + 1
            }
        }

        # Comparison and plotting the ROC curve
        state_2_binary = as.numeric(true_state == 2)

        c = seq(0,1,by=0.001)

        success_mat = matrix(nrow=length(c), ncol = length(state_2_binary))
        for(i in 1:nrow(success_mat)) {
            success_mat[i,] = as.numeric(cumulative_post_prob >= c[i])
        }

        true_positives = which(state_2_binary == 1)
        p = length(true_positives)                  # Number of true "positives"

        true_negatives = which(state_2_binary == 0)
        n = length(true_negatives)                  # Number of true "negatives"

        c_results = data.frame("c" = c, "true_pos" = rep(NA, length(c)),
                                "false_pos" = rep(NA, length(c)))

        for(i in 1:nrow(success_mat)) {
            test_positives = which(success_mat[i,] == 1)
            tpr = sum(test_positives %in% true_positives) / p
            test_negatives = which(success_mat[i,] == 0)
            fpr = 1 - (sum(test_negatives %in% true_negatives) / n)
            
            c_results[i,] = c(c[i], tpr, fpr)
        }

        # Reorder to be in increasing order for the false positive rate
        temp = c_results[order(c_results$false_pos), ]
        c_results = temp

        plot(c_results$false_pos, c_results$true_pos, xlim = c(0,1),
             xlab = "FPR", ylab = "TPR")

        # ----------------------------------------------------------------------
        sub_plot = c_results

        nonLinMod = lm(true_pos ~ ns(false_pos, 5), data = sub_plot)

        gridPoints = c_results$false_pos
        predictGrid = predict(nonLinMod, newdata=data.frame(false_pos=gridPoints))
        df = data.frame(gridPoints = gridPoints, predictGrid = predictGrid,
                        false_pos = sub_plot$false_pos, true_pos = sub_plot$true_pos)
        lines(gridPoints, predictGrid, col = "orange", lwd = 2)   
        # ----------------------------------------------------------------------

        # Determining the best 'c' by looking at the slope between adjacent lines
        delta_y = diff(c_results$true_pos)
        delta_x = diff(c_results$false_pos)
        slopes = delta_y / delta_x

        delta_y_smooth = diff(df$predictGrid)
        delta_x_smooth = diff(df$gridPoints)
        slopes_smooth = delta_y_smooth / delta_x_smooth

        print(c_results[which.min(abs(slopes_smooth - 1)), ])
        best_c[best_c_ind, ] = c_results[which.min(abs(slopes_smooth - 1)), ]
        points(c_results$false_pos[which.min(abs(slopes_smooth - 1))],
               c_results$true_pos[which.min(abs(slopes_smooth - 1))],
               col = 'red', cex = 2)
        points(df$gridPoints[which.min(abs(slopes_smooth - 1))],
               df$predictGrid[which.min(abs(slopes_smooth - 1))],
               col = 'green', cex = 2, pch = 2)
        
        title(main = paste0("Window: ", win_length, ", Seed: ", seed_ind, '\n',
                            "c: ", best_c[best_c_ind, 1],
                            ", TPR: ", round(best_c[best_c_ind, 2], 3),
                            ", FPR: ", round(best_c[best_c_ind, 3], 3)))
        best_c_ind = best_c_ind + 1
    }

    optimal_c_win[[www]] = best_c

}

save(optimal_c_win, file = "Plots/optimal_c_win_smooth.rda")