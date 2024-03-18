# Model evaluation to determine bleeding threshold and optimal time window.
library(matrixStats)
library(plotrix)
library(splines)

simulation = T

if(simulation) {
    trialNum = 1
    itNum = 5
    seed_list = 1
} else {
    trialNum = 8
    itNum = 5
}

# Load the data ---------------------------------------------------------------
if(simulation) {
    data_num = 1
    load(paste0('Data_sim/use_data1_', data_num, '.rda'))
} else {
    load('Data/data_format_new3.rda')
    use_data = data_format 
}

EIDs = unique(use_data[,"EID"])

# -----------------------------------------------------------------------------
# Testing different threshold values ------------------------------------------
# -----------------------------------------------------------------------------
window_length = 0:20
optimal_c_win = vector(mode = "list", length = length(window_length))

# Load the model output -------------------------------------------------------
B_chain = NULL

Dir = 'Model_out/'
for(seed_num in 1:length(seed_list)) {
    print(seed_list[seed_num])
    if(simulation) {
        load(paste0(Dir,'mcmc_out_interm_', seed_list[seed_num],'_', trialNum,'it',itNum,'_sim.rda'))   
    } else {
        load(paste0(Dir,'mcmc_out_interm_', seed_list[seed_num],'_', trialNum,'it',itNum,'.rda'))
    }
    
    if(seed_num == 1) {
        B_chain = mcmc_out_temp$B_chain
    } else {
        B_chain = rbind(B_chain, mcmc_out_temp$B_chain)
    }
}
#  ----------------------------------------------------------------------------

# # Revising the true state sequence to distinguish between bleeding and the FIRST
# # instance of bleeding.
# true_state = rep(NA, nrow(use_data))
# for(i in 1:length(EIDs)) {
#     indices_i = (use_data[,'EID']==EIDs[i])
#     sub_true_state = use_data[indices_i,'b_true']
# 
#     first_s_2 = T
#     for(j in 1:length(sub_true_state)) {
#         if(sub_true_state[j] == 2) {
#             if(first_s_2) {
#                 sub_true_state[j] = -2
#                 first_s_2 = F
#             }
#         } else {
#             first_s_2 = T
#         }
#     }
#     
#     true_state[indices_i] = sub_true_state
# }
true_state = use_data[,'b_true']

if(simulation) {
    plot_title = paste0("Plots/ROC_curve_", trialNum, "_it", itNum, "_sim.pdf")
} else {
    plot_title = paste0("Plots/ROC_curve_", trialNum, "_it", itNum, ".pdf")
}

pdf(plot_title)
par(mfrow=c(1,1))
for(www in 1:length(window_length)) {
    
    win_length = window_length[www]
    print(www)
    
    # Model output
    cumulative_post_prob = rep(NA, length(true_state))
    ind = 1
    for(i in unique(use_data[,'EID'])) {
        print(which(EIDs == i))
        indices_i = (use_data[,'EID']==i)
        indices_i_new = which(indices_i == T)
        for(w in 1:length(indices_i_new)) {
            start_index = indices_i_new[1]
            end_index = indices_i_new[w] 
            if(w - win_length > 0) start_index = indices_i_new[w - win_length]
            
            y_or_n_2 = apply(B_chain[, start_index:end_index, drop=F],
                             1, function(x) (2 %in% x))
            prob_2 = mean(y_or_n_2)
            
            cumulative_post_prob[ind] = prob_2
            ind = ind + 1
        }
    }
    
    # Comparison and plotting the ROC curve
    state_2_binary = as.numeric(true_state == 2)
    state_1_3_binary = as.numeric(true_state %in% c(1,3))
    
    c = seq(0,1,by=0.001)
    
    success_mat = matrix(nrow=length(c), ncol = length(state_2_binary))
    for(i in 1:nrow(success_mat)) {
        success_mat[i,] = as.numeric(cumulative_post_prob >= c[i])
    }
    
    true_positives = which(state_2_binary == 1)
    true_negatives = which(state_1_3_binary == 1)
    
    c_results = data.frame("c" = c, "true_pos"  = rep(NA, length(c)),
                                    "false_pos" = rep(NA, length(c)),
                                    "true_neg"  = rep(NA, length(c)),
                                    "false_neg" = rep(NA, length(c)))
    
    # Calculating the sensitivity and specificity information
    p = length(true_positives)  # Number of true "positives"
    n = length(true_negatives)  # Number of true "negatives"
    
    for(i in 1:nrow(success_mat)) {
        test_positives = which(success_mat[i,] == 1)
        test_negatives = which(success_mat[i,] == 0)
        
        tpr = sum(test_positives %in% true_positives) / p
        tnr = sum(test_negatives %in% true_negatives) / n
        
        fnr = 1 - tpr
        fpr = 1 - tnr
        
        c_results[i,] = c(c[i], tpr, fpr, tnr, fnr)
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
    points(c_results$false_pos[which.min(abs(slopes_smooth - 1))],
           c_results$true_pos[which.min(abs(slopes_smooth - 1))],
           col = 'red', cex = 2)
    points(df$gridPoints[which.min(abs(slopes_smooth - 1))],
           df$predictGrid[which.min(abs(slopes_smooth - 1))],
           col = 'green', cex = 2, pch = 2)
    
    title(main = paste0("Window: ", win_length, '\n',
                        "c: ", c_results[which.min(abs(slopes_smooth - 1)), 1],
                        ", TPR: ", round(c_results[which.min(abs(slopes_smooth - 1)), 2], 3),
                        ", FPR: ", round(c_results[which.min(abs(slopes_smooth - 1)), 3], 3),
                        ", TNR: ", round(c_results[which.min(abs(slopes_smooth - 1)), 4], 3),
                        ", FNR: ", round(c_results[which.min(abs(slopes_smooth - 1)), 5], 3)))
    optimal_c_win[[www]] = c_results
}
dev.off()

if(simulation) {
    save(optimal_c_win, file = paste0("Plots/optimal_c_win_", trialNum, 
                                      "_it", itNum, "_sim.rda"))
} else {
    save(optimal_c_win, file = paste0("Plots/optimal_c_win_", trialNum, 
                                      "_it", itNum, ".rda"))
}

# Calculating area under ROC curve to determine the best window and c
wind_area = matrix(nrow=length(optimal_c_win), ncol = 2)
wind_area[,1] = 1:length(optimal_c_win)
for(i in 1:length(optimal_c_win)) {
    c_results = optimal_c_win[[i]]
    
    temp = c_results[order(c_results$false_pos), ]
    c_results = temp
    area = 0
    
    # Estimate area using Trapezoid rule
    for(k in 2:nrow(c_results)) {
        area = area + 0.5 * (c_results$true_pos[k] + c_results$true_pos[k-1]) * 
                (c_results$false_pos[k] - c_results$false_pos[k-1])
    }
    
    wind_area[i,2] = area
}

print(wind_area)






