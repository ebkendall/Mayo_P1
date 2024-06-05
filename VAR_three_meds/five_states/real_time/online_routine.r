library(mvtnorm, quietly=T)
library(Matrix, quietly=T)
library(LaplacesDemon, quietly=T)
library(matrixStats)
library(plotrix)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
sourceCpp("online_func.cpp")
sourceCpp("../likelihood_fnc_arm.cpp")

# -----------------------------------------------------------------------------
# The mcmc algorithm
# -----------------------------------------------------------------------------
mcmc_online = function( par, par_index, A, W, B, Y, x, z, steps, burnin, ind, 
                        trialNum, Dn_omega, simulation, max_ind, df_num, n_cores, 
                        otype, pcov, pscale, s_time){

    EIDs = as.character(unique(Y[,'EID']))

    # Metropolis Parameter Index for MH within Gibbs updates
    # Ordering of the transition rate parameters:
    # 1->2, 1->4, 2->3, 2->4, 3->1, 3->2, 3->4, 4->2, 4->5, 5->1, 5->2, 5->4
    mpi = list(c(par_index$vec_init),
               c(par_index$vec_zeta[1:4]),
               c(par_index$vec_zeta[5:8]),
               c(par_index$vec_zeta[9:14]),
               c(par_index$vec_zeta[15:24]),
               c(par_index$vec_A[1:4]),
               c(par_index$vec_A[5:12]),
               c(par_index$vec_A[13:20]),
               c(par_index$vec_upsilon_omega[c(1:16, 36:57)]),
               c(par_index$vec_upsilon_omega[c(17:35, 58:88)]),
               c(par_index$vec_R))

    n_group = length(mpi)

    # # Reset pcov and pscale --------------------------------------------------
    # pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(mpi[[j]]))*.001
    # pscale = rep( 1, n_group)
  
    # Begin the MCMC algorithm -------------------------------------------------
    chain_length_MASTER = 5000
    chain = matrix( 0, chain_length_MASTER, length(par)) # steps
    B_chain = hc_chain = hr_chain = bp_chain = la_chain = matrix( 0, chain_length_MASTER, nrow(Y)) # steps-burnin
    accept = rep( 0, n_group)

    Dn_Xn = update_Dn_Xn_cpp( as.numeric(EIDs), B, Y, par, par_index, x, n_cores)
    Dn = Dn_Xn[[1]]; names(Dn) = EIDs
    Xn = Dn_Xn[[2]]
    
    # Keeping track of the sampled alpha_i
    A_chain = vector(mode = "list", length = length(EIDs))
    a_chain_id = 1:length(A_chain)
    
    for(a_ind in 1:length(A_chain)) {
        A_chain[[a_ind]] = matrix(nrow = length(par_index$vec_alpha_tilde), 
                                  ncol = chain_length_MASTER)
    }

    ttt = 0
    time_frame = T
    while(time_frame){
        ttt = ttt + 1
        
        # ----------------------------------------------------------------------
        # Updates constantly ---------------------------------------------------
        # ----------------------------------------------------------------------
        # Imputing the missing Y values
        Y = update_Y_i_cpp( as.numeric(EIDs), par, par_index, A, Y, Dn, Xn,
                            otype, Dn_omega, W, B, n_cores)
        colnames(Y) = c('EID','hemo', 'hr', 'map', 'lactate',
                        'RBC_rule', 'clinic_rule')
        # Gibbs updates of the alpha_i
        A = update_alpha_i_cpp( as.numeric(EIDs), par, par_index, Y, Dn, Xn,
                                Dn_omega, W, B, n_cores)
        names(A) = EIDs
        for(aaa in 1:length(a_chain_id)) {
            A_chain[[aaa]][,ttt] = A[[a_chain_id[aaa]]]
        }
        # Gibbs updates of the omega_i
        W = update_omega_i_cpp( as.numeric(EIDs), par, par_index, Y, Dn, Xn,
                                Dn_omega, A, B, n_cores)
        names(W) = EIDs

        # Metropolis-within-Gibbs update of the state space
        B_Dn = update_b_i_online(as.numeric(EIDs), par, par_index, A, B, Y, z, Dn,
                              Xn, Dn_omega, W, n_cores, x)
        B = B_Dn[[1]]; names(B) = EIDs
        Dn = B_Dn[[2]]; names(Dn) = EIDs
        # ----------------------------------------------------------------------
        # ----------------------------------------------------------------------
        

        # ----------------------------------------------------------------------
        # Updates every 100 steps ----------------------------------------------
        # ----------------------------------------------------------------------
        # if(ttt %% 100 == 0) {
        #     # Gibbs updates of the alpha~, omega~, beta, & Upsilon parameters
        #     par = update_beta_Upsilon_R_cpp( as.numeric(EIDs), par, par_index, A, Y,
        #                                      Dn, Xn, Dn_omega, W, B, n_cores)
        #     par = update_alpha_tilde_cpp( as.numeric(EIDs), par, par_index, A, Y)
        #     par = update_omega_tilde_cpp( as.numeric(EIDs), par, par_index, W, Y)
        #     
        #     # Save the parameter updates made in the Gibbs steps before Metropolis steps
        #     chain[ttt,] = par
        #     
        #     
        #     log_target_prev = log_post_cpp( as.numeric(EIDs), par, par_index, A, B,
        #                                     Y, z, Dn, Xn, Dn_omega, W, n_cores)
        #     
        #     if(!is.finite(log_target_prev)){
        #         print("Infinite log-posterior; Gibbs update went wrong")
        #         print(paste0("value: ", log_target_prev))
        #         break
        #     }
        #     
        #     # Metropolis-within-Gibbs updates -------------------------------------
        #     for(j in 1:n_group) {
        #         ind_j = mpi[[j]]
        #         
        #         proposal = par
        #         
        #         if(sum(ind_j %in% par_index$vec_R) == 0) {
        #             # logit_init, zeta, logit A1
        #             proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],
        #                                        sigma=pscale[[j]]*pcov[[j]])
        #             
        #             log_target = log_post_cpp( as.numeric(EIDs), proposal, par_index,
        #                                        A, B, Y, z, Dn, Xn, Dn_omega, W, n_cores)
        #             
        #             if(ttt < burnin){
        #                 while(!is.finite(log_target)){
        #                     print('bad proposal')
        #                     proposal = par
        #                     
        #                     proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],
        #                                                sigma=pcov[[j]]*pscale[j])
        #                     
        #                     log_target = log_post_cpp( as.numeric(EIDs), proposal,
        #                                                par_index, A, B, Y, z, Dn,
        #                                                Xn, Dn_omega, W, n_cores)
        #                 }
        #             }
        #             
        #             if(!is.finite(log_target) | is.nan(log_target)) {
        #                 # Ensuring that we do not have problems from C++
        #                 print(paste0("bad proposal post burnin: ", log_target))
        #                 log_target = -Inf
        #             }
        #             
        #             if( log_target - log_target_prev > log(runif(1,0,1)) ){
        #                 log_target_prev = log_target
        #                 par[ind_j] = proposal[ind_j]
        #                 accept[j] = accept[j] +1
        #             }
        #             
        #             chain[ttt,ind_j] = par[ind_j]
        #             
        #         } else {
        #             # Updating R ---------------------------------------------------
        #             
        #             # Prior for R ---------------------------
        #             nu_R = 10
        #             psi_R = diag(c(4.58, 98.2, 101.3, 7.6))
        #             psi_R = (nu_R - 4 - 1) * psi_R
        #             
        #             # q(R* | R(t)) -----------------------------
        #             curr_R = matrix(par[ind_j], nrow = 4)
        #             psi_nu_star_t = proposal_R_cpp_new(nu_R, psi_R, curr_R, Y, Dn, Xn, A, par, par_index, as.numeric(EIDs), B, Dn_omega, W)
        #             q_s_star_t = psi_nu_star_t[[1]] / pscale[j]
        #             q_nu_star_t = floor(psi_nu_star_t[[2]] / pscale[j])
        #             
        #             # Proposal R -----------------------------
        #             proposal[ind_j] = c(rinvwishart(nu = q_nu_star_t, S = q_s_star_t))
        #             
        #             # q(R(t) | R*) -----------------------------
        #             prop_R = matrix(proposal[ind_j], nrow = 4)
        #             psi_nu_t_star = proposal_R_cpp_new(nu_R, psi_R, prop_R, Y, Dn, Xn, A, par, par_index, as.numeric(EIDs), B, Dn_omega, W)
        #             q_s_t_star = psi_nu_t_star[[1]] / pscale[j]
        #             q_nu_t_star = floor(psi_nu_t_star[[2]] / pscale[j])
        #             
        #             log_prop      = dinvwishart(Sigma = prop_R, nu = q_nu_star_t, S = q_s_star_t, log = T)
        #             log_prop_prev = dinvwishart(Sigma = curr_R, nu = q_nu_t_star, S = q_s_t_star, log = T)
        #             
        #             log_target = log_post_cpp( as.numeric(EIDs), proposal, par_index, A, B, Y, z, Dn, Xn, Dn_omega, W, n_cores)
        #             
        #             if(!is.finite(log_target + log_prop_prev - log_target_prev - log_prop) | is.nan(log_target + log_prop_prev - log_target_prev - log_prop)) {
        #                 # Ensuring that we do not have problems from C++
        #                 print("Current R")
        #                 print(curr_R)
        #                 print("Prop R")
        #                 print(prop_R)
        #                 print(paste0("log_target: ", log_target))
        #                 print(paste0("log_prop_prev: ", log_prop_prev))
        #                 print(paste0("log_target_prev: ", log_target_prev))
        #                 print(paste0("log_prop: ", log_prop))
        #                 print(log_target + log_prop_prev - log_target_prev - log_prop)
        #             }
        #             
        #             if( log_target + log_prop_prev - log_target_prev - log_prop > log(runif(1,0,1)) ){
        #                 log_target_prev = log_target
        #                 par[ind_j] = proposal[ind_j]
        #                 accept[j] = accept[j] +1
        #             }
        #             
        #             if (ttt == 1 | ttt %% 100 == 0){
        #                 print(paste0('likelihood: ', log_target_prev))
        #             }
        #         }
        #     }
        # }
        # ----------------------------------------------------------------------
        # ----------------------------------------------------------------------
        
        
        # Restart the acceptance ratio at burnin
        if(ttt == burnin) accept = rep( 0, n_group)
        if(ttt > burnin){
            B_chain[ ttt,] = do.call( 'c', B)
            hc_chain[ ttt,] = Y[,'hemo']
            hr_chain[ ttt,] = Y[,'hr']
            bp_chain[ ttt,] = Y[,'map']
            la_chain[ ttt,] = Y[,'lactate']
        }
        # ----------------------------------------------------------------------

        if(ttt%%1==0)  cat('--->',ttt,'\n')
        
        e_time = proc.time();
        if(e_time[3] - s_time[3] > 15) time_frame = F
    }
    # ---------------------------------------------------------------------------

    return(list( chain=chain, B_chain=B_chain, hc_chain=hc_chain, hr_chain=hr_chain, 
                bp_chain=bp_chain, la_chain = la_chain, A_chain = A_chain, otype=otype,
                accept=accept/(steps-burnin), pscale=pscale, pcov = pcov,
                par_index=par_index, ttt = ttt))
}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

initialize_Y <- function(Y, EIDs) {
    
    mean_vital_val = c(9.991089, 84.823810, 78.106992, 3.980250)
    
    for(i in EIDs) {
        heading_names = c('hemo', 'hr', 'map', 'lactate')
        sub_dat = Y[Y[,"EID"] == i, ]
        
        for(k in 1:length(heading_names)) {
            if(sum(is.na(sub_dat[,heading_names[k]])) == nrow(sub_dat)) {
                sub_dat[,heading_names[k]] = mean_vital_val[k]
            } else {
                if(sum(!is.na(sub_dat[,heading_names[k]])) == 1) {
                    sub_dat[,heading_names[k]] = sub_dat[!is.na(sub_dat[,heading_names[k]]), heading_names[k]]
                } else {
                    obs_indices = which(!is.na(sub_dat[,heading_names[k]]))
                    miss_indices = which(is.na(sub_dat[,heading_names[k]]))
                    for(j in miss_indices) {
                        if(j < obs_indices[1]) {
                            sub_dat[j,heading_names[k]] = sub_dat[obs_indices[1], heading_names[k]]
                        } else if(j > tail(obs_indices,1)) {
                            sub_dat[j,heading_names[k]] = sub_dat[tail(obs_indices,1), heading_names[k]]
                        } else {
                            end_pts = c(max(obs_indices[obs_indices < j]),
                                        min(obs_indices[obs_indices > j]))
                            slope = (sub_dat[end_pts[2], heading_names[k]] - sub_dat[end_pts[1], heading_names[k]]) / diff(end_pts)
                            sub_dat[j,heading_names[k]] = slope * (j - end_pts[1]) + sub_dat[end_pts[1], heading_names[k]]
                        }
                    }
                }
            }
            Y[Y[,"EID"] == i, heading_names[k]] = sub_dat[,heading_names[k]]
        }
    }
    
    return(Y)
}

add_data = function(it_num, pred_id, data_format, b_curr, Y_curr_EIDs) {
    
    new_rows = 4 + it_num
    
    df_sub = data_format[data_format[,"EID"] %in% pred_id, ]
    
    df_sub_start = NULL
    no_more_data = T
    for(j in unique(df_sub[,"EID"])) {
        df_j = df_sub[df_sub[,'EID'] == j, ]
        if(new_rows <= nrow(df_j)) {
            no_more_data = F
            df_sub_start = rbind(df_sub_start, df_j[1:new_rows, ])   
        } else {
            df_sub_start = rbind(df_sub_start, df_j)
        }
    }
    
    if(no_more_data) {
        # No more new rows to learn
        Y = matrix(nrow=0,ncol=0)
        EIDs = x = z = otype = NULL
        
        return(list( Y=Y, EIDs=EIDs, x = x, z = z, otype = otype))
    }
    
    Y = df_sub_start[, c('EID','hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
    EIDs = as.character(unique(Y[,'EID']))
    
    x = df_sub_start[,c('n_RBC_admin'), drop=F]
    p = ncol(x)
    
    z = cbind(1, df_sub_start[,c('RBC_ordered'), drop=F])
    m = ncol(z)
    
    # 1 = observed, 0 = missing
    otype = !is.na(Y[, c('hemo','hr','map','lactate')])
    colnames(otype) = c('hemo','hr','map','lactate')
    
    A = list()
    W = list()
    B = list()
    
    for(i in EIDs){
        
        A[[i]] = matrix(par[par_index$vec_alpha_tilde], ncol =1)
        W[[i]] = matrix(par[par_index$omega_tilde], ncol =1)
        
        b_temp_init = b_curr[Y_curr_EIDs == i]
        b_temp = c(b_temp_init, b_temp_init[length(b_temp_init)])
        
        B[[i]] = matrix(b_temp, ncol = 1)
    }
    
    return(list( Y=Y, EIDs=EIDs, x = x, z = z, otype = otype, 
                 A = A, W = W, B = B, df_sub_start = df_sub_start))
}

plotting_fnc <- function(use_data, EIDs, B_chain, Hc_chain, Hr_chain, 
                         Map_chain, La_chain, it_num, prev_state_avg) {
    
    pdf_title = paste0("Plots/iteration_", it_num, ".pdf")
    pdf(pdf_title)
    panel_dim = c(4,1)
    inset_dim = c(0,-.18)
    par(mfrow=panel_dim, mar=c(2,4,2,4), bg='black', fg='green')
    for(i in EIDs){
    
        indices_i = (use_data[,'EID']==i)
        n_i = sum(indices_i)
        
        t_grid = round(use_data[indices_i, 'time'] / 60, digits = 3)
        t_grid_bar = 1:length(t_grid)
        rbc_times_bar = which(use_data[use_data[,'EID']==i, 'RBC_ordered'] != 0)
        
        rbc_admin_times_bar = which(use_data[use_data[,'EID']==i, 'RBC_admin'] != 0) 
        
        rbc_times = t_grid[rbc_times_bar]
        rbc_admin_times = t_grid[rbc_admin_times_bar]
        
        # Computing the posterior probability averages -------------------------
        curr_post_prob = rbind(colMeans(B_chain[, indices_i] == 1),
                               colMeans(B_chain[, indices_i] == 2),
                               colMeans(B_chain[, indices_i] == 3),
                               colMeans(B_chain[, indices_i] == 4),
                               colMeans(B_chain[, indices_i] == 5))
        prev_post_prob = prev_state_avg[2:6, prev_state_avg[1,] == i]
        curr_post_prob[,1:ncol(prev_post_prob)] = 0.5 * (curr_post_prob[,1:ncol(prev_post_prob)] + prev_post_prob)
        
        pb = barplot(curr_post_prob, 
                     col=c( 'dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray'), 
                     xlab='time', space=0, col.main='green', border=NA, axes = F, plot = F) 
        
        # Heart Rate and MAP double plot -----------------------------------------
        if(mean(use_data[indices_i, 'clinic_rule']) != 0) {
            title_name = paste0('Heart Rate & MAP: ', i, ', RBC Rule = ', 
                                mean(use_data[indices_i, 'RBC_rule']),
                                ', clinic = ', mean(use_data[indices_i, 'clinic_rule']))
        } else {
            title_name = paste0('Heart Rate & MAP: ', i, ', RBC Rule = ', 
                                mean(use_data[indices_i, 'RBC_rule']))
        }
        
        hr_upper = colQuantiles( Hr_chain[, indices_i, drop=F], probs=.975)
        hr_lower = colQuantiles( Hr_chain[, indices_i, drop=F], probs=.025)
        bp_upper = colQuantiles( Map_chain[, indices_i, drop=F], probs=.975)
        bp_lower = colQuantiles( Map_chain[, indices_i, drop=F], probs=.025)
        
        hr_map_ylim = c(min(hr_lower, bp_lower), max(hr_upper, bp_upper))
        
        # Make a new plot to add the background color
        plot(NULL, xlim=range(pb) + c(-0.5,0.5), ylim=hr_map_ylim, main=title_name,
             xlab='time', ylab=NA, xaxt='n', col.main='green',
             col.axis='green')
        
        plotCI( x = pb, y=colMeans(Hr_chain[, indices_i, drop=F]), ui=hr_upper, li=hr_lower,
                main=title_name,
                xlab='time', ylab=NA, xaxt='n', col.main='green',
                col.axis='green', pch=20, cex=1, sfrac=.0025, col = 'aquamarine',
                xlim = range(pb) + c(-0.5,0.5), ylim = hr_map_ylim, add =T) 
        plotCI( x = pb, y=colMeans(Map_chain[, indices_i, drop=F]), ui=bp_upper, li=bp_lower,
                main=title_name,
                xlab='time', ylab=NA, xaxt='n', pch=20, cex=1, sfrac=.0025,
                col = 'orange',
                xlim = range(pb) + c(-0.5,0.5), add = T) 
        legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
                legend=c( 'HR', 'MAP'), pch=15, pt.cex=1.5, 
                col=c( 'aquamarine', 'orange'))
        # grid( nx=20, NULL, col='white')
        grid(nx = NA, ny = NULL, col = "white")
        axis( side=1, at=pb, col.axis='green', labels=t_grid)
        
        abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
        abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
        
        
        # Hemoglobin and Lactate double plot -------------------------------------
        if(mean(use_data[indices_i, 'clinic_rule']) != 0) {
            title_name = paste0('Hemoglobin & Lactate: ', i, ', RBC Rule = ', 
                                mean(use_data[indices_i, 'RBC_rule']),
                                ', clinic = ', mean(use_data[indices_i, 'clinic_rule']))
        } else {
            title_name = paste0('Hemoglobin & Lactate: ', i, ', RBC Rule = ',
                                mean(use_data[indices_i, 'RBC_rule']))
        }
        
        hc_upper = colQuantiles( Hc_chain[, indices_i, drop=F], probs=.975)
        hc_lower = colQuantiles( Hc_chain[, indices_i, drop=F], probs=.025)
        la_upper = colQuantiles( La_chain[, indices_i, drop=F], probs=.975)
        la_lower = colQuantiles( La_chain[, indices_i, drop=F], probs=.025)
        
        hr_map_ylim = c(min(hc_lower, la_lower), max(hc_upper, la_upper))
        
        plot(NULL, xlim=range(pb) + c(-0.5,0.5), ylim=hr_map_ylim, main=title_name,
             xlab='time', ylab=NA, xaxt='n', col.main='green',
             col.axis='green')
        
        plotCI(x = pb, y = colMeans(Hc_chain[, indices_i, drop=F]), ui=hc_upper, li=hc_lower,
               main=title_name,
               xlab='time', ylab=NA, xaxt='n', col.main='green',
               col.axis='green', pch=20, cex=1, sfrac=.0025, col = 'aquamarine',
               xlim = range(pb) + c(-0.5,0.5), ylim = hr_map_ylim, add = T) 
        plotCI( x = pb, y=colMeans(La_chain[, indices_i, drop=F]), ui=la_upper, li=la_lower,
                main=title_name,
                xlab='time', ylab=NA, xaxt='n', pch=20, cex=1, sfrac=.0025,
                col = 'orange',
                xlim = range(pb) + c(-0.5,0.5), add = T) 
        legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
                legend=c( 'hemo', 'lactate'), pch=15, pt.cex=1.5, 
                col=c( 'aquamarine', 'orange'))
        # grid( nx=20, NULL, col='white')
        grid(nx = NA, ny = NULL, col = "white")
        axis( side=1, at=pb, col.axis='green', labels=t_grid)
        
        abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
        abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
        
        # BAR PLOTS --------------------------------------------------------------
        barplot(curr_post_prob, 
                col=c( 'dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray'), 
                xlab='time', space=0, col.main='green', border=NA,
                xlim=range(pb) + c(-0.5,0.5)) 
        # grid( nx=20, NULL, col='white')
        grid(nx = NA, ny = NULL, col = "white")
        legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
                legend=c( 'Baseline', 'State 2', 'State 3', 'State 4', 'State 5'), 
                pch=15, pt.cex=1.5, 
                col=c( 'dodgerblue', 'firebrick1', 'yellow2','green', 'darkgray'))
        legend( 'topleft', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
                legend=c( 'RBC order', 'RBC admin'), pch=15, pt.cex=1.5, 
                col=c( 'darkorchid1', 'aquamarine'))				
        axis( side=1, at=t_grid_bar-0.5, col.axis='green', labels = t_grid)
        axis( side=2, at=0:1, col.axis='green')
        
        
        abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
        abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
        
        # Cumulative PLOTS ---------------------------------------------------------
        cumulative_post_prob = matrix(nrow = 2, ncol = n_i)
        ind = 1
        win_length = 0
        c = 0.257
        
        indices_i_new = which(indices_i == T)
        for(w in 1:length(indices_i_new)) {
            start_index = indices_i_new[1]
            end_index = indices_i_new[w] 
            if(w - win_length > 0) start_index = indices_i_new[w - win_length]
            
            y_or_n_2 = apply(B_chain[, start_index:end_index, drop=F],
                             1, function(x) (2 %in% x))
            prob_2 = mean(y_or_n_2)
            
            cumulative_post_prob[, ind] = c(prob_2, 1-prob_2)
            ind = ind + 1
        }
        
        barplot( cumulative_post_prob,
                 col=c('darkred', 'black'),
                 main=paste0('cumulative prob.'), xlab='time', space=0, 
                 col.main='green', border=NA,
                 xlim=range(pb) + c(-0.5,0.5))
        # grid( nx=20, NULL, col='white')
        grid(nx = NA, ny = NULL, col = "white")
        legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
                legend=c('', 'State 2'), pch=15, pt.cex=1.5,
                col=c('black', 'darkred'))
        axis( side=1, at=t_grid_bar-0.5, col.axis='green', labels = t_grid)
        axis( side=2, at=0:1, col.axis='green')
        
        abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
        abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
        abline(h = c, col = 'yellow', lwd = 2)
    }
    dev.off()
}
