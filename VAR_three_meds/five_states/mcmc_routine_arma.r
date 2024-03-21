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
mcmc_routine = function( par, par_index, A, W, B, Y, x, z, steps, burnin, ind, 
                         trialNum, Dn_omega, simulation, bleed_indicator, max_ind){
  
    n_cores = strtoi(Sys.getenv(c("LSB_DJOB_NUMPROC")))

    print(paste0("Number of cores: ", n_cores))

    EIDs = as.character(unique(Y[,'EID']))

    # Index of observed versus missing data
    # 1 = observed, 0 = missing
    otype = !is.na(Y[, c('hemo','hr','map','lactate')])
    colnames(otype) = c('hemo','hr','map','lactate')

    # Metropolis Parameter Index for MH within Gibbs updates
    # Ordering of the transition rate parameters:
    # 1->2, 1->4, 2->3, 2->4, 3->1, 3->2, 3->4, 4->2, 4->5, 5->1, 5->2, 5->4
    mpi = list(c(par_index$vec_init),
               c(par_index$vec_zeta),
               c(par_index$vec_A[1:4]),
               c(par_index$vec_A[5:12]),
               c(par_index$vec_A[13:20]),
               # c(par_index$vec_upsilon_omega[c(1,2,4,7,8,11,12,13,15,23,24,27,
               #                                 30,32)]),
               # c(par_index$vec_upsilon_omega[c(3,5,6,9,10,14,16,17,18,19,20,21,
               #                                 22,25,26,28,29,31,33,34,35,36)]),
               # c(par_index$vec_upsilon_omega[c(38,40,42,43,44,45,47,48,50,52,55,
               #                                 56,57,58,60,61,62,64,66,67,70,71,
               #                                 73,75,76,77,78,79,80,83,84,85,86,
               #                                 87,88)]),                                        
               # c(par_index$vec_upsilon_omega[c(37,39,41,46,49,51,53,54,59,63,65,
               #                                 68,69,72,74,81,82)]),
               c(par_index$vec_upsilon_omega[c(1:16)]),
               c(par_index$vec_upsilon_omega[c(17:35)]),
               c(par_index$vec_upsilon_omega[c(36:57)]),
               c(par_index$vec_upsilon_omega[c(58:88)]),
               c(par_index$vec_R))

    n_group = length(mpi)

    # Loading an existing pcov and pscale ------------------------------
    pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(mpi[[j]]))*.001
    pscale = rep( 1, n_group)

    if(!simulation) {
        print('Real data analysis')
        load(paste0('Model_out/mcmc_out_interm_', 3, '_1it', 2, '.rda'))
        # pcov = mcmc_out_temp$pcov
        # pscale = mcmc_out_temp$pscale
        
        # Setting initial values for Y
        load(paste0('../Data/data_format_new', real_dat_num, '.rda'))
        old_EIDs = unique(data_format[,"EID"])
        old_time = data_format[,"time"]
        old_ids = data_format[,"EID"]
        
        load('Data_updates/data_format.rda')
        new_EIDs = unique(data_format[,'EID'])
        
        data_format = data_format[data_format[,"EID"] %in% old_EIDs, ]
    
        for(i in EIDs){
            old_t = old_time[old_ids == as.numeric(i)]
            curr_t = data_format[data_format[,"EID"] == as.numeric(i), "time"]
            if(sum(floor(old_t) %in% floor(curr_t)) == length(old_t)) {
                Y[Y[,'EID'] == as.numeric(i), 'hemo'] = c(mcmc_out_temp$hc_chain[nrow(mcmc_out_temp$hc_chain), old_ids == as.numeric(i)])
                Y[Y[,'EID'] == as.numeric(i), 'hr'] = c(mcmc_out_temp$hr_chain[nrow(mcmc_out_temp$hr_chain), old_ids == as.numeric(i)])
                Y[Y[,'EID'] == as.numeric(i), 'map'] = c(mcmc_out_temp$bp_chain[nrow(mcmc_out_temp$bp_chain), old_ids == as.numeric(i)])
                Y[Y[,'EID'] == as.numeric(i), 'lactate'] = c(mcmc_out_temp$la_chain[nrow(mcmc_out_temp$la_chain), old_ids == as.numeric(i)])
            } else {
                b_index = which(floor(old_t) %in% floor(curr_t))
                b_temp_init = which(old_ids == as.numeric(i))
                b_temp = b_temp_init[b_index]
                
                Y[Y[,'EID'] == as.numeric(i), 'hemo'] = c(mcmc_out_temp$hc_chain[nrow(mcmc_out_temp$hc_chain), b_temp])
                Y[Y[,'EID'] == as.numeric(i), 'hr'] = c(mcmc_out_temp$hr_chain[nrow(mcmc_out_temp$hr_chain), b_temp])
                Y[Y[,'EID'] == as.numeric(i), 'map'] = c(mcmc_out_temp$bp_chain[nrow(mcmc_out_temp$bp_chain), b_temp])
                Y[Y[,'EID'] == as.numeric(i), 'lactate'] = c(mcmc_out_temp$la_chain[nrow(mcmc_out_temp$la_chain), b_temp])
            }
        }
        
        rm(mcmc_out_temp)
        rm(data_format)
        # Y[, 'hemo'] = c(mcmc_out_temp$hc_chain[nrow(mcmc_out_temp$hc_chain), ])
        # Y[, 'hr'] = c(mcmc_out_temp$hr_chain[nrow(mcmc_out_temp$hr_chain), ])
        # Y[, 'map'] = c(mcmc_out_temp$bp_chain[nrow(mcmc_out_temp$bp_chain), ])
        # Y[, 'lactate'] = c(mcmc_out_temp$la_chain[nrow(mcmc_out_temp$la_chain), ])
    } 
  
    # Begin the MCMC algorithm -------------------------------------------------
    chain_length_MASTER = 1001
    chain = matrix( 0, chain_length_MASTER, length(par)) # steps
    B_chain = hc_chain = hr_chain = bp_chain = la_chain = matrix( 0, chain_length_MASTER, nrow(Y)) # steps-burnin
    accept = rep( 0, n_group)

    Dn_Xn = update_Dn_Xn_cpp( as.numeric(EIDs), B, Y, par, par_index, x, n_cores)
    Dn = Dn_Xn[[1]]; names(Dn) = EIDs
    Xn = Dn_Xn[[2]]
    
    # Keeping track of the sampled alpha_i
    A_chain = vector(mode = "list", length = 10)
    a_chain_id = c(3, 86, 163, 237, 427, 521, 632, 646, 692, 713)
    
    for(a_ind in 1:10) {
        A_chain[[a_ind]] = matrix(nrow = 20, ncol = chain_length_MASTER)
    }

    for(ttt in 1:steps){

        chain_ind = ttt %% 10000
        if(chain_ind == 0) chain_ind = 10000

        # Thinning the saved chain index
        chain_ind = floor(chain_ind / 10) + 1

        if(!simulation) {
            # Imputing the missing Y values ----------------------------------------
            # print("Update Y"); s_time = Sys.time()
            Y = update_Y_i_cpp( as.numeric(EIDs), par, par_index, A, Y, Dn, Xn, 
                                otype, Dn_omega, W, B, n_cores)
            colnames(Y) = c('EID','hemo', 'hr', 'map', 'lactate', 
                            'RBC_rule', 'clinic_rule')
            # e_time = Sys.time() - s_time; print(e_time)
        }

        # Gibbs updates of the alpha_i -----------------------------------------
        # print("Update alpha_i"); s_time = Sys.time()
        A = update_alpha_i_cpp( as.numeric(EIDs), par, par_index, Y, Dn, Xn, 
                                Dn_omega, W, B, n_cores)
        names(A) = EIDs
        # e_time = Sys.time() - s_time; print(e_time)

        for(aaa in 1:length(a_chain_id)) {
            A_chain[[aaa]][,chain_ind] = A[[a_chain_id[aaa]]]
        }
        
        # Gibbs updates of the omega_i -----------------------------------------
        # print("Update omega_i"); s_time = Sys.time()
        W = update_omega_i_cpp( as.numeric(EIDs), par, par_index, Y, Dn, Xn, 
                                Dn_omega, A, B, n_cores)
        names(W) = EIDs
        # e_time = Sys.time() - s_time; print(e_time)

        # Metropolis-within-Gibbs update of the state space --------------------
        # print("Update b_i"); s_time = Sys.time()
        B_Dn = update_b_i_cpp(as.numeric(EIDs), par, par_index, A, B, Y, z, Dn, 
                              Xn, Dn_omega, W, bleed_indicator, n_cores)
        B = B_Dn[[1]]; names(B) = EIDs
        Dn = B_Dn[[2]]; names(Dn) = EIDs
        # e_time = Sys.time() - s_time; print(e_time)

        # Gibbs updates of the alpha~, omega~, beta, & Upsilon parameters ------
        # print("Update beta_upsilon"); s_time = Sys.time()
        par = update_beta_Upsilon_R_cpp( as.numeric(EIDs), par, par_index, A, Y, 
                                         Dn, Xn, Dn_omega, W, B, n_cores)
        # e_time = Sys.time() - s_time; print(e_time)
        # print("Update alpha tilde"); s_time = Sys.time()
        par = update_alpha_tilde_cpp( as.numeric(EIDs), par, par_index, A, Y)
        # e_time = Sys.time() - s_time; print(e_time)
        # print("Update omega tilde"); s_time = Sys.time()
        par = update_omega_tilde_cpp( as.numeric(EIDs), par, par_index, W, Y)
        # e_time = Sys.time() - s_time; print(e_time)

        # Save the parameter updates made in the Gibbs steps before Metropolis steps
        chain[chain_ind,] = par

        # Printing updates ----------------------------------------------------
        if (ttt %% 100 == 0){
            print("alpha_tilde")
            print(round(chain[chain_ind, par_index$vec_alpha_tilde], 3))
            
            print("diag of Sigma_Upsilon")
            Sigma_t = matrix(chain[chain_ind,par_index$vec_sigma_upsilon], ncol = 20)
            print(round(diag(Sigma_t), 3))
            
            print("A1")
            vec_A_t_logit = chain[chain_ind, par_index$vec_A]
            vec_A_t = exp(vec_A_t_logit) / (1 + exp(vec_A_t_logit))
            mat_A_t = matrix(vec_A_t, nrow = 4)
            print(mat_A_t)
            print(matrix(vec_A_t_logit, nrow = 4))
            
            print("R")
            R_t = matrix(chain[chain_ind, par_index$vec_R], ncol = 4)
            print(R_t)
            
            print("zeta")
            zed = matrix(chain[chain_ind, par_index$vec_zeta], nrow = 2)
            print(zed)
            
            print("omega_tilde")
            print(chain[chain_ind, par_index$omega_tilde])
            
            print("log upsilon omega")
            print(chain[chain_ind, par_index$vec_upsilon_omega])
            
            print("acceptance")
            print(accept / (ttt %% 480))
        }
        # ---------------------------------------------------------------------
        
        log_target_prev = log_post_cpp( as.numeric(EIDs), par, par_index, A, B, 
                                        Y, z, Dn, Xn, Dn_omega, W, n_cores)

        if(!is.finite(log_target_prev)){
            print("Infinite log-posterior; Gibbs update went wrong")
            print(paste0("value: ", log_target_prev))
            break
        }

        # print("Update metropolis step"); s_time = Sys.time()
        # Metropolis-within-Gibbs updates -------------------------------------
        for(j in 1:n_group) {
            ind_j = mpi[[j]]
            proposal = par
            
            if(sum(ind_j %in% par_index$vec_R) == 0) {
                # logit_init, zeta, logit A1
                proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j], 
                                           sigma=pscale[[j]]*pcov[[j]])
                
                log_target = log_post_cpp( as.numeric(EIDs), proposal, par_index, 
                                           A, B, Y, z, Dn, Xn, Dn_omega, W, n_cores)

                if(ttt < burnin){
                    while(!is.finite(log_target)){
                        print('bad proposal')
                        proposal = par

                        proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],
                                                   sigma=pcov[[j]]*pscale[j])
                        
                        log_target = log_post_cpp( as.numeric(EIDs), proposal, 
                                                   par_index, A, B, Y, z, Dn, 
                                                   Xn, Dn_omega, W, n_cores)
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
                
                # Proposal tuning scheme ---------------------------------------
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
                # --------------------------------------------------------------
                
            } else {
                # Updating R ---------------------------------------------------
                # Changing the proposal distribution and therefore the Metrop Ratio
                # emp_R_est = diag(4.58, 98.2, 101.3, 7.6)
                nu_R = 100
                psi_R = diag(c(4.58, 98.2, 101.3, 7.6))
                psi_R = (nu_R - 4 - 1) * psi_R
                
                curr_psi_nu = proposal_R_cpp(nu_R, psi_R, Y, Dn, Xn, A, par, par_index, as.numeric(EIDs), B, Dn_omega, W)
                
                q_s  = curr_psi_nu[[1]]
                q_nu = curr_psi_nu[[2]]
                
                proposal[ind_j] = c(rinvwishart(nu = q_nu, S = q_s))
                
                curr_R = matrix(par[ind_j], nrow = 4)
                prop_R = matrix(proposal[ind_j], nrow = 4)
                
                log_prop = dinvwishart(Sigma = prop_R, nu = q_nu, S = q_s, log = T)
                log_prop_prev = dinvwishart(Sigma = curr_R, nu = q_nu, S = q_s, log = T)
                
                log_target = log_post_cpp( as.numeric(EIDs), proposal, par_index, A, B, Y, z, Dn, Xn, Dn_omega, W, n_cores)

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
        # e_time = Sys.time() - s_time; print(e_time)

        # Restart the acceptance ratio at burnin
        if(ttt == burnin) accept = rep( 0, n_group)
        if(ttt > burnin){
            B_chain[ chain_ind,] = do.call( 'c', B)
            hc_chain[ chain_ind,] = Y[,'hemo']
            hr_chain[ chain_ind,] = Y[,'hr']
            bp_chain[ chain_ind,] = Y[,'map']
            la_chain[ chain_ind,] = Y[,'lactate']
        }
        # ----------------------------------------------------------------------

        if(ttt%%1==0)  cat('--->',ttt,'\n')
        if(ttt > burnin & ttt%%10000==0) {
            mcmc_out_temp = list( chain=chain, B_chain=B_chain, hc_chain=hc_chain, 
                                hr_chain=hr_chain, bp_chain=bp_chain, 
                                la_chain = la_chain, A_chain = A_chain,
                                otype=otype, accept=accept/length(burnin:ttt), 
                                pscale=pscale, pcov = pcov, par_index=par_index)
            if(simulation) {
                save(mcmc_out_temp, file = paste0('Model_out/mcmc_out_interm_',ind,'_', 
                                                  trialNum, 'it', ttt/10000 + (max_ind - 5), '_sim.rda'))
            } else {
                save(mcmc_out_temp, file = paste0('Model_out/mcmc_out_interm_',ind,'_', 
                                                  trialNum, 'it', ttt/10000 + (max_ind - 5), '.rda'))
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

b_ind_fnc <- function(data_format) {
    bleed_pat = unique(data_format[data_format[,"RBC_rule"] != 0, "EID"])
    bleed_indicator = rep(0, nrow(data_format))
    for(i in 1:length(bleed_pat)) {
        sub_dat = data_format[data_format[,"EID"] == bleed_pat[i], ]
        
        # Check in any 12 hour period
        max_time = tail(sub_dat[,"time"], 1)
        when_rbc = c(1, which(diff(sub_dat[,"n_RBC_admin"]) != 0))
        
        for(j in 1:length(when_rbc)) {
            s_time = sub_dat[when_rbc[j], "time"]
            e_time_12 = s_time + 720
            e_time_24 = s_time + 1440
            RBC_diff_12 = RBC_diff_24 = 0
            
            if (e_time_12 <= max_time) {
                s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
                ind_12 = order(abs(sub_dat[,"time"] - e_time_12))[1]
                RBC_diff_12 = sub_dat[ind_12, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
            } else {
                s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
                e_ind = order(abs(sub_dat[,"time"] - max_time))[1]
                RBC_diff_12 = sub_dat[e_ind, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
            }
            if (e_time_24 <= max_time) {
                s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
                ind_24 = order(abs(sub_dat[,"time"] - e_time_24))[1]
                RBC_diff_24 = sub_dat[ind_24, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
            } else {
                s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
                e_ind = order(abs(sub_dat[,"time"] - max_time))[1]
                RBC_diff_24 = sub_dat[e_ind, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
            }
            
            if(RBC_diff_12 >=3 | RBC_diff_24 >= 6) {
                admin_times = sub_dat[sub_dat[,"RBC_admin"] != 0, "time"]
                if(RBC_diff_12 >=3) {
                    a_t = which(admin_times >= s_time & admin_times < e_time_12)
                    first_time = admin_times[a_t[1]]
                    order_times = sub_dat[sub_dat[,"RBC_ordered"] != 0, "time"]
                    if(sum(order_times <= first_time) == 0) {
                        print(paste0(i, ", ", sub_dat[1,"EID"]))
                        first_order_time = first_time
                    } else {
                        first_order_time = max(order_times[order_times <= first_time])   
                    }
                } else if (RBC_diff_24 >= 6) {
                    a_t = which(admin_times >= s_time & admin_times < e_time_24)
                    first_time = admin_times[a_t[1]]  
                    order_times = sub_dat[sub_dat[,"RBC_ordered"] != 0, "time"]
                    if(sum(order_times <= first_time) == 0) {
                        print(paste0(i, ", ", sub_dat[1,"EID"]))
                        first_order_time = first_time
                    } else {
                        first_order_time = max(order_times[order_times <= first_time])   
                    }
                }
                
                bleed_indicator[data_format[,"EID"] == bleed_pat[i] & 
                                    data_format[,"time"] == first_order_time] = 1
                break
            }
            
        }
    }
    return(bleed_indicator)
}


var_R_calc <- function(psi, nu, p) {
    var_mat = matrix(0, p, p)
    for(i in 1:p) {
        for(j in 1:p) {
            num_ij = (nu - p + 1) * (psi[i,j] * psi[i,j]) + (nu - p - 1) * psi[i,i] * psi[j,j]
            den_ij = (nu - p) * (nu - p - 1) * (nu - p - 1) * (nu - p - 3)
            var_mat[i,j] = num_ij / den_ij
        }
    }
    return(var_mat)
}
