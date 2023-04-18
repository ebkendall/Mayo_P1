library(mvtnorm, quietly=T)
library(Matrix, quietly=T)
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
mcmc_routine = function( par, par_index, A, W, B, Y, x, z, steps, burnin, ind, trialNum, med_format){
  
  EIDs = as.character(unique(Y[,'EID']))
  
  # Index of observed versus missing data
  # 1 = observed, 0 = missing
  otype = !is.na(Y[, c('hemo','hr','map','lactate')])
  colnames(otype) = c('hemo','hr','map','lactate')
  
  Xn = list()
  for(i in EIDs)  Xn[[i]] = diag(4) %x% x[ Y[,'EID']==as.numeric(i),, drop=F]
  
  # Metropolis Parameter Index for MH within Gibbs updates
  # Ordering of the transition rate parameters:
  # 1->2, 1->4, 2->3, 2->4, 3->1, 3->2, 3->4, 4->2, 4->5, 5->1, 5->2, 5->4
  mpi = list( c(par_index$log_theta, par_index$vec_init), 
              c(par_index$vec_zeta[1:2]), c(par_index$vec_zeta[3:4]),
              c(par_index$vec_zeta[5:6]), c(par_index$vec_zeta[7:8]),
              c(par_index$log_lambda[c(1,4,7,10)]),
              c(par_index$log_lambda[c(2,5,8,11)]),
              c(par_index$log_lambda[c(3,6,9,12)]))
  n_group = length(mpi)
  
  # Loading an existing pcov and pscale ------------------------------
  # load('Model_out/mcmc_out_interm_3_12it5.rda')
  # pcov = mcmc_out_temp$pcov
  # pscale = mcmc_out_temp$pscale
  # rm(mcmc_out_temp)
  pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(mpi[[j]]))*.001
  pscale = rep( 1, n_group)
  
  # Begin the MCMC algorithm -------------------------------------------------
  chain_length_MASTER = 1001
  chain = matrix( 0, chain_length_MASTER, length(par)) # steps
  B_chain = hc_chain = hr_chain = bp_chain = la_chain = matrix( 0, chain_length_MASTER, nrow(Y)) # steps-burnin
  A_chain = vector(mode = "list", length = 10)
  accept = rep( 0, n_group)
  
  # # Debugging strategy -------------------------------------------------------
  # # id's of interest: 100950, 747025
  # debug_info_100950 = vector(mode = "list", length = 2)
  # debug_info_747025 = vector(mode = "list", length = 2)
  # # The first index will be a matrix will all state updates
  # debug_info_100950[[1]] = matrix(nrow = 5000, ncol = sum(Y[,'EID']==100950))
  # debug_info_747025[[1]] = matrix(nrow = 5000, ncol = sum(Y[,'EID']==747025))
  # # The second index will track the likelihood before, after, the proposed state, and whether it accepted
  # debug_info_100950[[2]] = vector(mode = "list", length = 5000)
  # debug_info_747025[[2]] = vector(mode = "list", length = 5000)
  # # --------------------------------------------------------------------------
  
  Dn = update_Dn_cpp( as.numeric(EIDs), B, Y)
  names(Dn) = EIDs
  
  Dn_omega = list()
  for(i in 1:length(EIDs)) Dn_omega[[i]] = diag(4) %x% med_format[[i]]
  names(Dn_omega) = EIDs
  
  for(ttt in 1:steps){
    
    chain_ind = ttt %% 10000
    if(chain_ind == 0) chain_ind = 10000

    # Thinning the saved chain index
    chain_ind = floor(chain_ind / 10) + 1
    A_check = 100
    
    # Update the inverse OU covariance matrix (function of theta)
    invKn = update_invKn_cpp( as.numeric(EIDs), par, par_index, Y)
    names(invKn) = EIDs
    
    Y = update_Y_i_cpp( as.numeric(EIDs), par, par_index, A, Y, Dn, Xn, invKn, otype, Dn_omega, W)
    colnames(Y) = c('EID','hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')

    # Gibbs updates of the alpha_i
    A = update_alpha_i_cpp( as.numeric(EIDs), par, par_index, Y, Dn, Xn, invKn, Dn_omega, W)
    names(A) = EIDs
    
    # Gibbs updates of the omega_i
    W = update_omega_i_cpp( as.numeric(EIDs), par, par_index, Y, Dn, Xn, invKn, Dn_omega, A)
    names(W) = EIDs
    
    if(chain_ind %% A_check == 0) {
      A_chain[[chain_ind / A_check]] = A
    }
    
    # Debug information -------------------------------------------------------
    debug_temp1 = matrix(nrow = 5, ncol = sum(Y[,'EID']==100950) - 1)
    debug_temp2 = matrix(nrow = 5, ncol = sum(Y[,'EID']==747025) - 1)
    # -------------------------------------------------------
    
    # Metropolis-within-Gibbs update of the state space
    B_Dn = update_b_i_cpp(16, as.numeric(EIDs), par, par_index, A, B, Y, z, Dn, Xn, invKn, Dn_omega, W,
                          debug_temp1, debug_temp2)
    B = B_Dn[[1]]; names(B) = EIDs
    Dn = B_Dn[[2]]; names(Dn) = EIDs
    
    # # Debug information -------------------------------------------------------
    # if(ttt %% 5000 == 0) {
    #     final_debug = list("l_100950" = debug_info_100950,
    #                        "l_747025" = debug_info_747025)
    #     save(final_debug, file = paste0('Model_out/final_debug',ind,'_it', ttt/5000, '_11.rda'))
        
    #     debug_info_100950[[1]] = matrix(nrow = 5000, ncol = sum(Y[,'EID']==100950))
    #     debug_info_747025[[1]] = matrix(nrow = 5000, ncol = sum(Y[,'EID']==747025))
    #     debug_info_100950[[2]] = vector(mode = "list", length = 5000)
    #     debug_info_747025[[2]] = vector(mode = "list", length = 5000)
    # } else {
    #     debug_info_100950[[1]][ttt %% 5000, ] = c(B[['100950']])
    #     debug_info_747025[[1]][ttt %% 5000, ] = c(B[['747025']])
        
    #     debug_info_100950[[2]][[ttt %% 5000]] = B_Dn[[3]]
    #     y_sub = t(Y[Y[,'EID'] == 100950, c('hemo', 'hr', 'map', 'lactate')])
    #     debug_info_100950[[2]][[ttt %% 5000]] = rbind(debug_info_100950[[2]][[ttt %% 5000]],
    #                                                   y_sub[,-1])
        
    #     debug_info_747025[[2]][[ttt %% 5000]] = B_Dn[[4]]
    #     y_sub = t(Y[Y[,'EID'] == 747025, c('hemo', 'hr', 'map', 'lactate')])
    #     debug_info_747025[[2]][[ttt %% 5000]] = rbind(debug_info_747025[[2]][[ttt %% 5000]],
    #                                                   y_sub[,-1])
    # }
    # # -------------------------------------------------------
    
    # Gibbs updates of the alpha_tilde, beta, Upsilon, & R parameters
    par = update_beta_Upsilon_R_cpp( as.numeric(EIDs), par, par_index, A, Y, Dn, Xn, invKn, Dn_omega, W) 
    par = update_alpha_tilde_cpp( as.numeric(EIDs), par, par_index, A, Y)
    par = update_omega_tilde_cpp( as.numeric(EIDs), par, par_index, W, Y)
    
    # Save the parameter updates made in the Gibbs steps before Metropolis steps
    chain[chain_ind,] = par
    if (ttt %% 100 == 0) print(round(chain[chain_ind,], 3))

    # Testing Data for development of functions ------------------------------
    # test_post = list("EIDs" = EIDs, "par" = par, "par_index" = par_index,
    #                  "A" = A, "B" = B, "Y" = Y, "z" = z, "Dn" = Dn,
    #                  "Xn" = Xn, "invKn" = invKn, "otype" = otype)
    # save(test_post, file = "test_post_24.rda")
    # return(0);
    
    # Metropolis-within-Gibbs update of the theta and zeta parameters
    for(j in 1:n_group) {
      
      ind_j = mpi[[j]]
      proposal = par
      proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j], sigma=pscale[[j]]*pcov[[j]])
      
      log_target_prev = log_post_cpp( as.numeric(EIDs), par, par_index, A, B, Y, z, Dn, Xn, invKn, Dn_omega, W)
      
      log_target = log_post_cpp( as.numeric(EIDs), proposal, par_index, A, B, Y, z, Dn, Xn, invKn, Dn_omega, W)
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
        
        if(100 <= ttt & ttt <= 2000){
          temp_chain = chain[1:(floor(ttt/10) + 1),ind_j]
          pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
          
        } else if(2000 < ttt){
          temp_chain = chain[(floor((ttt-2000) / 10) + 1):(floor(ttt/10) + 1),ind_j]
          pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
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
        save(mcmc_out_temp, file = paste0('Model_out/mcmc_out_interm_',ind,'_', 
                                        trialNum, 'it', ttt/10000, '.rda'))
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