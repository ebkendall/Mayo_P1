library(mvtnorm, quietly=T)
library(Matrix, quietly=T)
library(LaplacesDemon, quietly=T)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
sourceCpp("likelihood_fnc_arm.cpp")

# Needed for OpenMP C++ parallel
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")


# -----------------------------------------------------------------------------
# The mcmc algorithm
# -----------------------------------------------------------------------------
mcmc_routine = function( par, par_index, A, W, B, Y, x, z, steps, burnin, ind, trialNum, Dn_omega, simulation){
  
    EIDs = as.character(unique(Y[,'EID']))
    
    # Index of observed versus missing data
    # 1 = observed, 0 = missing
    otype = !is.na(Y[, c('hemo','hr','map','lactate')])
    colnames(otype) = c('hemo','hr','map','lactate')
    
    # Metropolis Parameter Index for MH within Gibbs updates
    # Ordering of the transition rate parameters:
    # 1->2, 1->4, 2->3, 2->4, 3->1, 3->2, 3->4, 4->2, 4->5, 5->1, 5->2, 5->4
    mpi = list( c(par_index$vec_init),
                c(par_index$vec_zeta),
                # c(par_index$log_lambda)
                c(par_index$vec_A[1:4]),
                c(par_index$vec_A[5:8]),
                c(par_index$vec_A[9:12]),
                c(par_index$vec_R))

    n_group = length(mpi)
    
    # Loading an existing pcov and pscale ------------------------------
    pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(mpi[[j]]))*.001
    pscale = rep( 1, n_group)

    if(!simulation) {
        print('Real data analysis')
        load('Model_out/mcmc_out_interm_4_3it5.rda')

        pcov = mcmc_out_temp$pcov
        pscale = mcmc_out_temp$pscale

        # Setting initial values for Y
        # Y[, 'hemo'] = c(mcmc_out_temp$hc_chain[1000, ])
        # Y[, 'hr'] = c(mcmc_out_temp$hr_chain[1000, ])
        # Y[, 'map'] = c(mcmc_out_temp$bp_chain[1000, ])
        # Y[, 'lactate'] = c(mcmc_out_temp$la_chain[1000, ])
        rm(mcmc_out_temp)
        # Setting initial values for Y
        print("Initializing the missing Y's for imputation")
        for(i in EIDs) {
            heading_names = c('hemo','hr','map','lactate')
            sub_dat = Y[Y[,"EID"] == i, ]
            
            for(k in 1:length(heading_names)) {
                if(sum(is.na(sub_dat[,heading_names[k]])) == nrow(sub_dat)) {
                    sub_dat[,heading_names[k]] = mean(Y[,heading_names[k]], na.rm =T)
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
        print("Done initializing")
    }
  
  # Begin the MCMC algorithm -------------------------------------------------
  chain_length_MASTER = 1001
  chain = matrix( 0, chain_length_MASTER, length(par)) # steps
  B_chain = hc_chain = hr_chain = bp_chain = la_chain = matrix( 0, chain_length_MASTER, nrow(Y)) # steps-burnin
  A_chain = vector(mode = "list", length = 10)
  accept = rep( 0, n_group)
  
  Dn_Xn = update_Dn_Xn_cpp( as.numeric(EIDs), B, Y, par, par_index, x)
  Dn = Dn_Xn[[1]]; names(Dn) = EIDs
  Xn = Dn_Xn[[2]]
  
  # Dn_omega = list()
  # for(i in 1:length(EIDs)) Dn_omega[[i]] = diag(c(1,1,1,1)) %x% med_format[[i]]
  # names(Dn_omega) = EIDs
  
  for(ttt in 1:steps){
    
    chain_ind = ttt %% 10000
    if(chain_ind == 0) chain_ind = 10000

    # Thinning the saved chain index
    chain_ind = floor(chain_ind / 10) + 1
    A_check = 100
    
    Y = update_Y_i_cpp( as.numeric(EIDs), par, par_index, A, Y, Dn, Xn, otype, Dn_omega, W, B)
    colnames(Y) = c('EID','hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')

    # Gibbs updates of the alpha_i
    A = update_alpha_i_cpp( as.numeric(EIDs), par, par_index, Y, Dn, Xn, Dn_omega, W, B)
    names(A) = EIDs
    
    # # Gibbs updates of the omega_i
    # W = update_omega_i_cpp( as.numeric(EIDs), par, par_index, Y, Dn, Xn, Dn_omega, A)
    # names(W) = EIDs
    
    if(chain_ind %% A_check == 0) {
      A_chain[[chain_ind / A_check]] = A
    }
    
    # Debug information -------------------------------------------------------
    # debug_temp1 = matrix(nrow = 7, ncol = sum(Y[,'EID']==1) - 1)
    # debug_temp2 = matrix(nrow = 7, ncol = sum(Y[,'EID']==2) - 1)
    debug_temp1 = matrix(nrow = 7, ncol = 50)
    debug_temp2 = matrix(nrow = 7, ncol = 50)
    # -------------------------------------------------------
    
    # Metropolis-within-Gibbs update of the state space (*** VAR UPDATED ***)
    B_Dn = update_b_i_cpp(16, as.numeric(EIDs), par, par_index, A, B, Y, z, Dn, Xn, Dn_omega, W,
                          debug_temp1, debug_temp2)
    B = B_Dn[[1]]; names(B) = EIDs
    Dn = B_Dn[[2]]; names(Dn) = EIDs

    # Gibbs updates of the alpha_tilde, beta, Upsilon, & R parameters (*** VAR UPDATED ***)
    par = update_beta_Upsilon_R_cpp( as.numeric(EIDs), par, par_index, A, Y, Dn, Xn, Dn_omega, W, B)
    par = update_alpha_tilde_cpp( as.numeric(EIDs), par, par_index, A, Y)
    # par = update_omega_tilde_cpp( as.numeric(EIDs), par, par_index, W, Y)
    
    # Save the parameter updates made in the Gibbs steps before Metropolis steps
    chain[chain_ind,] = par

    # Printing updates --------------------------------------------------------
    if (ttt %% 100 == 0){
      print("alpha_tilde")
      print(round(chain[chain_ind, par_index$vec_alpha_tilde], 3))
      
      print("diag of Sigma_Upsilon")
      Sigma_t = matrix(chain[chain_ind,par_index$vec_sigma_upsilon], ncol = 12)
      # Lambda_t = diag(exp(chain[chain_ind,par_index$log_lambda]))
      # Upsilon_t = Lambda_t %*% Sigma_t %*% Lambda_t
      print(round(diag(Sigma_t), 3))
      
      print("A1")
      vec_A_t_logit = chain[chain_ind, par_index$vec_A]
      vec_A_t = exp(vec_A_t_logit) / (1 + exp(vec_A_t_logit))
      mat_A_t = matrix(vec_A_t, nrow = 4)
      print(mat_A_t)
      
      print("R")
      R_t = matrix(chain[chain_ind, par_index$vec_R], ncol = 4)
      print(R_t)
      
      print("zeta")
      zed = matrix(chain[chain_ind, par_index$vec_zeta], nrow = 2)
      print(zed)
      
      print("acceptance")
      print(accept / (ttt %% 480))
    }
    # -------------------------------------------------------------------------
    
    log_target_prev = log_post_cpp( as.numeric(EIDs), par, par_index, A, B, Y, z, Dn, Xn, Dn_omega, W)

    if(!is.finite(log_target_prev)){
        print("Infinite log-posterior; Gibbs update went wrong")
        print(paste0("value: ", log_target_prev))
        break
    }
    
    # Metropolis-within-Gibbs update of the theta and zeta parameters
    for(j in 1:n_group) {
        ind_j = mpi[[j]]
        proposal = par
      
        if(j <= 5) {
            # logit_init, zeta, log_lambda, logit A1
            proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j], sigma=pscale[[j]]*pcov[[j]])
            
            log_target = log_post_cpp( as.numeric(EIDs), proposal, par_index, A, B, Y, z, Dn, Xn, Dn_omega, W)

            if(ttt < burnin){
                while(!is.finite(log_target)){
                    print('bad proposal')
                    proposal = par

                    proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],sigma=pcov[[j]]*pscale[j])
                    
                    log_target = log_post_cpp( as.numeric(EIDs), proposal, par_index, A, B, Y, z, Dn, Xn, Dn_omega, W)
                }
            }
      
            if(!is.finite(log_target) | is.nan(log_target)) {
                # Ensuring that we do not have problems from C++
                print(paste0("bad proposal post burnin: ", log_target))
                log_target = -Inf
            }
          
            if( log_target - log_target_prev > log(runif(1,0,1)) ){
                log_target_prev = log_target
                par[ind_j] = proposal[ind_j]
                accept[j] = accept[j] +1
            }
          
            chain[chain_ind,ind_j] = par[ind_j]
          
            # Proposal tuning scheme ------------------------------------------------
            if(ttt < burnin){
                # During the burnin period, update the proposal covariance in each step
                # to capture the relationships within the parameters vectors for each
                # transition.  This helps with mixing.
                if(ttt == 100)  pscale[j] = 1
                
                if (length(ind_j) > 1) {
                    if(100 <= ttt & ttt <= 2000){
                        temp_chain = chain[1:(floor(ttt/10) + 1),ind_j]
                        pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
                        
                    } else if(2000 < ttt){
                        temp_chain = chain[(floor((ttt-2000) / 10) + 1):(floor(ttt/10) + 1),ind_j]
                        pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
                    }
                } else {
                    if(100 <= ttt & ttt <= 2000){
                        temp_chain = chain[1:(floor(ttt/10) + 1),ind_j]
                        pcov[[j]] = matrix(var(temp_chain[ !duplicated(temp_chain)]))
                        
                    } else if(2000 < ttt){
                        temp_chain = chain[(floor((ttt-2000) / 10) + 1):(floor(ttt/10) + 1),ind_j]
                        pcov[[j]] = matrix(var(temp_chain[ !duplicated(temp_chain)]))
                    }
                }
                
                if( sum( is.na(pcov[[j]]) ) > 0)  pcov[[j]] = diag( length(ind_j) )
                
                # Tune the proposal covariance for each transition to achieve
                # reasonable acceptance ratios.
                if(ttt %% 30 == 0){
                    if(ttt %% 480 == 0){
                        accept[j] = 0
                        
                    } else if( accept[j] / (ttt %% 480) < .4 ){
                        pscale[j] = (.75^2)*pscale[j]
                        
                    } else if( accept[j] / (ttt %% 480) > .5 ){
                        pscale[j] = (1.25^2)*pscale[j]
                    }
                }
            }
            # -----------------------------------------------------------------------
          
        } else {
            # Updating R ----------------------------------------------------
            # Changing the proposal distribution and therefore the Metrop Ratio
            nu_R = 6
            psi_R = diag(4)
            
            curr_psi_nu = proposal_R_cpp(nu_R, psi_R, Y, Dn, Xn, A, par, par_index, as.numeric(EIDs), B)
            
            proposal[ind_j] = c(rinvwishart(nu = curr_psi_nu[[2]],
                                            S = curr_psi_nu[[1]]))
            
            curr_R = matrix(par[ind_j], nrow = 4)
            prop_R = matrix(proposal[ind_j], nrow = 4)
            
            log_prop = dinvwishart(Sigma = prop_R, nu = curr_psi_nu[[2]],
                                    S = curr_psi_nu[[1]], log = T)
            log_prop_prev = dinvwishart(Sigma = curr_R, nu = curr_psi_nu[[2]],
                                        S = curr_psi_nu[[1]], log = T)
            
            log_target = log_post_cpp( as.numeric(EIDs), proposal, par_index, A, B, Y, z, Dn, Xn, Dn_omega, W)

            if(!is.finite(log_target + log_prop_prev - log_target_prev - log_prop) | is.nan(log_target + log_prop_prev - log_target_prev - log_prop)) {
                # Ensuring that we do not have problems from C++
                print("Current R")
                print(curr_R)
                print("Prop R")
                print(prop_R)
                print(paste0("log_target: ", log_target))
                print(paste0("log_prop_prev: ", log_prop_prev))
                print(paste0("log_target_prev: ", log_target_prev))
                print(paste0("log_prop: ", log_prop))
                print(log_target + log_prop_prev - log_target_prev - log_prop)
            }
            
            if( log_target + log_prop_prev - log_target_prev - log_prop > log(runif(1,0,1)) ){
                log_target_prev = log_target
                par[ind_j] = proposal[ind_j]
                accept[j] = accept[j] +1
            }
        }

    }
    
    # Restart the acceptance ratio at burnin
    if(ttt == burnin) accept = rep( 0, n_group)
    if(ttt > burnin){
      B_chain[ chain_ind,] = do.call( 'c', B)
      hc_chain[ chain_ind,] = Y[,'hemo']
      hr_chain[ chain_ind,] = Y[,'hr']
      bp_chain[ chain_ind,] = Y[,'map']
      la_chain[ chain_ind,] = Y[,'lactate']
    }
    # -------------------------------------------------------------------------
    
    if(ttt%%1==0)  cat('--->',ttt,'\n')
    if(ttt > burnin & ttt%%10000==0) {
      mcmc_out_temp = list( chain=chain, B_chain=B_chain, hc_chain=hc_chain, 
                            hr_chain=hr_chain, bp_chain=bp_chain, 
                            la_chain = la_chain, A_chain = A_chain,
                            otype=otype, accept=accept/length(burnin:ttt), 
                            pscale=pscale, pcov = pcov, par_index=par_index)
      if(ttt/10000 > 1) {
        if(simulation) {
          save(mcmc_out_temp, file = paste0('Model_out/mcmc_out_interm_',ind,'_', 
                                        trialNum, 'it', ttt/10000, '_sim.rda'))
        } else {
          save(mcmc_out_temp, file = paste0('Model_out/mcmc_out_interm_',ind,'_', 
                                        trialNum, 'it', ttt/10000, '.rda'))
        }
      }
      # Reset the chains
      chain = matrix( NA, chain_length_MASTER, length(par)) 
      B_chain = hc_chain = hr_chain = bp_chain = la_chain = matrix( NA, chain_length_MASTER, nrow(Y))
    }
  }
  # ---------------------------------------------------------------------------
  
  return(list( chain=chain, B_chain=B_chain, hc_chain=hc_chain, hr_chain=hr_chain, 
               bp_chain=bp_chain, la_chain = la_chain, A_chain = A_chain, otype=otype,
               accept=accept/(steps-burnin), pscale=pscale, pcov = pcov,
               par_index=par_index))
}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------