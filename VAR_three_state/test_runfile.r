library(MASS, quietly=T)
library(mvtnorm, quietly=T)
library(LaplacesDemon, quietly=T)
library(Matrix, quietly=T)
library(RcppArmadillo)
library(RcppDist)
library(Matrix, quietly=T)
library(rbenchmark)
library(foreach, quietly=T)
library(doParallel, quietly=T)
library(bayesSurv)

library(Rcpp)
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")

# load('test_post.rda')
# load('test_list.rda')
Rcpp::sourceCpp("likelihood_fnc_arm.cpp")
set.seed(5)


# -----------------------------------------------------------------------------
# Older functions
# -----------------------------------------------------------------------------
log_f_i = function( i, t_pts, par, par_index, A, B, Y, z, Dn, Xn, invKn){
  
  in_value = 0
  
  vec_beta = matrix( par[par_index$vec_beta], ncol=1)
  R = matrix(par[par_index$vec_R], ncol = 4)
  zeta = matrix( par[par_index$vec_zeta], ncol=4)
  
  # The time-homogeneous probability transition matrix
  z_i = z[ Y[,'EID']==as.numeric(i),, drop=F]
  n_i = sum(Y[,'EID']==as.numeric(i))
  
  # The state transition likelihood component for current iterate of b_i
  b_i = B[[i]]
  init_logit = c(1, exp(par[par_index$vec_init][1]), exp(par[par_index$vec_init][2]))
  P_i = init_logit / sum(init_logit) # Initial state vector will need to be estimated
  # Full likelihood evaluation is not needed for updating pairs of b_i compnents
  if(any(t_pts=='full'))  t_pts = 1:n_i
  for(k in t_pts){
    if(k==1){
      in_value = in_value + log(P_i[b_i[1,]])
    } else{
      # Transition rate from nominal to bleeding
      q12 = exp(z_i[k,, drop=F] %*% zeta[,1, drop=F])
      # Transition rate from bleeding to recovering
      q23 = exp(z_i[k,, drop=F] %*% zeta[,2, drop=F])
      # Transition rate from recovering to nominal
      q31 = exp(z_i[k,, drop=F] %*% zeta[,3, drop=F])
      # Transition rate from recovering to bleeding
      q32 = exp(z_i[k,, drop=F] %*% zeta[,4, drop=F])
      
      Q = matrix(c(   1, q12,   0,
                      0,   1, q23,
                      q31, q32,   1), ncol=3, byrow=T)
      P_i = Q / rowSums(Q)
      
      in_value = in_value + log(P_i[ b_i[k-1,], b_i[k,]])
    }
  }
  
  
  # The likelihood component for hemo, hr, and map
  vecY_i = matrix( Y[ Y[,'EID']==as.numeric(i), c('hemo', 'hr', 'map', 'lactate')], ncol=1)
  vec_alpha_i = matrix( A[[i]], ncol=1)
  
  dev = vecY_i - Dn[[i]]%*%vec_alpha_i - Xn[[i]]%*%vec_beta
  precision = solve(R) %x% invKn[[i]]
  
  log_det1 = -n_i * determinant(R, logarithm=T)$modulus
  log_det2 = 4 * determinant(invKn[[i]], logarithm=T)$modulus
  log_det_precision = log_det1 + log_det2
  
  in_value = in_value + .5*( log_det_precision - t(dev)%*%precision%*%dev )
  
  return(as.numeric(in_value))
}

log_post = function( EIDs, par, par_index, A, B, Y, z, Dn, Xn, invKn){
  
  # Compute the likelihood ----------------------------------------------------
  value = foreach( i=EIDs, .combine='+', .export='log_f_i') %dopar% {
    return(log_f_i(i, 'full', par, par_index, A, B, Y, z, Dn, Xn, invKn))
  }
  print(value)
  # ---------------------------------------------------------------------------
  
  # Compute prior densities of all necessary model parameters -----------------
  theta = exp(par[par_index$log_theta])
  prior_theta = dlnorm( x=theta, meanlog=0, sdlog=0.5, log=T)
  print("Prior Theta")
  print(prior_theta)
  
  zeta = matrix( par[par_index$vec_zeta], ncol=4)
  m = nrow(zeta)
  # zeta_mean = c(matrix( c(-8.1708, -3.8595, -3.8595, -6.042217, rep(0, 4*(m-1))), ncol=4, byrow=T))
  zeta_mean = c(matrix( c(-8.1708, -2.652, -10.77766, -8.1708, rep(0, 4*(m-1))), ncol=4, byrow=T)) 
  zeta_sd = 1*diag(4*m)
  
  zeta_sd[1,1] = 0.5; zeta_sd[4,4] = 0.237
  zeta_sd[7,7] = 0.158; zeta_sd[10,10] = 0.229
  prior_zeta = dmvnorm( x=c(zeta), mean= zeta_mean, sigma=zeta_sd,log=T)
  print("Prior Zeta")
  print(prior_zeta)
  # ---------------------------------------------------------------------------
  
  value = value + prior_theta + prior_zeta
  return(value)
}

update_Dn = function( EIDs, B, Y){
  
  Dn = foreach(i=EIDs) %dopar% {
    
    n_i = sum(Y[,'EID']==as.numeric(i))
    b_i = B[[i]]
    D_i = diag(4) %x% cbind( 1, cumsum(b_i==2), cumsum(b_i==3))
    
    return(D_i)
  }
  names(Dn) = EIDs
  
  return(Dn)
}

update_alpha_i = function( EIDs, par, par_index, Y, Dn, Xn, invKn){
  
  invR = solve( matrix(par[par_index$vec_R], ncol = 4) )
  vec_alpha_tilde = matrix( par[par_index$vec_alpha_tilde], ncol=1)
  vec_beta = matrix( par[par_index$vec_beta], ncol=1)
  inv_Upsilon = solve(matrix( par[par_index$vec_Upsilon], ncol=12))
  
  A = foreach( i=EIDs, .packages='bayesSurv') %dopar% {		
    
    vecY_i = matrix( Y[ Y[,'EID']==as.numeric(i), c('hemo', 'hr', 'map','lactate')], ncol=1)
    
    hold = t(Dn[[i]]) %*% ( invR %x% invKn[[i]] )
    V_i = as.matrix( hold%*%(vecY_i - Xn[[i]]%*%vec_beta) 
                     + inv_Upsilon%*%vec_alpha_tilde)
    
    inv_W_i = as.matrix(hold %*% Dn[[i]] + inv_Upsilon)
    W_i = solve(inv_W_i)
    
    vec_alpha_i = rMVNorm( n=1, mean= W_i %*% V_i, Q=inv_W_i, param='standard') 
    
    return(c(vec_alpha_i))
  }
  names(A) = EIDs
  
  return(A)
}

update_alpha_tilde = function( EIDs, par, par_index, A, Y){
  
  # The prior mean for vec_alpha_tilde
  vec_alpha_tilde_0 = c(matrix( c(  10, 80, 80,   1,
                                    -1,  1, -1,   1,
                                    1, -1,  1,-0.5), ncol=4, byrow=T))
  
  # The prior precision matrix for vec_alpha_tilde
  inv_Sigma_alpha = diag(1/c(matrix( c( 1, 20, 20, 1,
                                        1, 10, 10, 1,
                                        1, 10, 10, 1)^2, ncol=4, byrow=T)))
  
  inv_Upsilon = solve(matrix( par[par_index$vec_Upsilon], ncol=12))
  
  inv_U = inv_Upsilon*length(EIDs) + inv_Sigma_alpha
  U = solve(inv_U)
  
  sum_alpha = matrix( rowSums(do.call( 'cbind', A)), ncol=1)
  
  hold = inv_Sigma_alpha %*% vec_alpha_tilde_0 + inv_Upsilon %*% sum_alpha
  
  par[par_index$vec_alpha_tilde] = rMVNorm( n=1, mean= U %*% hold, Q=inv_U, 
                                            param='standard')
  
  return(par)
}

update_beta_Upsilon_R = function( EIDs, par, par_index, A, Y, Dn, Xn, invKn){
    
    # The prior mean for vec_beta
    vec_beta_0 = matrix( 0, length(par_index$vec_beta), 1)
    
    # The prior precision matrix for vec_beta
    inv_Sigma_beta = diag(length(par_index$vec_beta)) / 10^2
    
    # The prior degrees of freedom for Upsilon
    nu_Upsilon = length(EIDs)
    
    # The prior scale matrix for Upsilon
    Psi_Upsilon = .0001 * diag(20)
    
    # The prior degrees of freedom for R
    nu_R = 100
    
    # The prior scale matrix for R
    Psi_R = 1 * diag(4)
    
    R = matrix(par[par_index$vec_R], ncol = 4)
    
    vec_alpha_tilde = matrix( par[par_index$vec_alpha_tilde], ncol=1)
    vec_beta = matrix( par[par_index$vec_beta], ncol=1)
    
    parallel_out = foreach( i=EIDs) %dopar% {			
        
        vecY_i = matrix( Y[ Y[,'EID']==as.numeric(i), c('hemo', 'hr', 'map','lactate')], ncol=1)
        vec_alpha_i = matrix( A[[i]], ncol=1)
        
        hold1 = t(Xn[[i]]) %*% ( solve(R) %x% invKn[[i]] )
        in_V = as.matrix(hold1 %*% (vecY_i - Dn[[i]]%*%vec_alpha_i))
        in_inv_W = as.matrix(hold1 %*% Xn[[i]])
        
        hold2 = vec_alpha_i - vec_alpha_tilde
        in_Upsilon_cov = hold2 %*% t(hold2)
        
        Y_i = Y[ Y[,'EID']==as.numeric(i),c('hemo', 'hr', 'map','lactate'), drop=F]
        vec_E_Y_i = Dn[[i]]%*%vec_alpha_i + Xn[[i]]%*%vec_beta
        E_Y_i = matrix(vec_E_Y_i, ncol = 4)
        hold3 = Y_i - E_Y_i
        in_R_cov = t(hold3) %*% invKn[[i]] %*% hold3
        
        return(list( in_V, in_inv_W, in_Upsilon_cov, in_R_cov))
    }
    # 'foreach' returns a list in which each component, i, is a list of 3 
    # components.  For each of these 3 components, the sum over all i=EIDs
    # is needed.  This is achieved as follows.
    
    V = inv_Sigma_beta %*% vec_beta_0 + Reduce( '+', sapply(parallel_out,'[',1))
    
    inv_W = inv_Sigma_beta + Reduce( '+', sapply(parallel_out,'[',2))
    W = solve(inv_W)
    
    Upsilon_cov = Psi_Upsilon + Reduce( '+', sapply(parallel_out,'[',3))
    
    R_cov = Psi_R + Reduce( '+', sapply(parallel_out,'[',4))
    N = nrow(Y)
    
    par[par_index$vec_beta] = rMVNorm( n=1,mean= W %*% V,Q=inv_W,param='standard')
    par[par_index$vec_Upsilon]= c( rinvwishart(nu=nu_Upsilon+1, S=Upsilon_cov) )
    par[par_index$vec_R]= c( rinvwishart(nu=nu_R + N, S=R_cov) )
    
    return(par)
}

update_Y_i = function( EIDs, par, par_index, A, Y, Dn, Xn, invKn, otype){
    
    newY = foreach( i=EIDs, .combine='rbind', .packages='bayesSurv') %dopar% {
        
        vec_beta = matrix( par[par_index$vec_beta], ncol=1)
        R = matrix(par[par_index$vec_R], ncol = 4)
        
        # Index of observed versus missing data
        # 1 = observed, 0 = missing
        otype_i = c( otype[ Y[,'EID']==as.numeric(i),] )
        vecY_i = matrix( Y[ Y[,'EID']==as.numeric(i), c('hemo','hr','map','lactate')], ncol=1)
        vec_alpha_i = matrix( A[[i]], ncol=1)
        
        loc = Dn[[i]]%*%vec_alpha_i + Xn[[i]]%*%vec_beta 
        loc_0 = loc[ otype_i==0,, drop=F]
        loc_1 = loc[ otype_i==1,, drop=F]
        
        precision = solve(R) %x% invKn[[i]]
        cond_prec = as.matrix(precision[ otype_i==0, otype_i==0])
        block_12 = precision[ otype_i==0, otype_i==1]
        
        dev = vecY_i[ otype_i==1,, drop=F] - loc_1
        cond_loc = as.matrix(loc_0 - chol2inv(chol(cond_prec)) %*% block_12 %*% dev)
        
        vecY_i[otype_i==0] = rMVNorm(n=1,mean=cond_loc,Q=cond_prec,param='standard')
        
        return(matrix( vecY_i, ncol=4))
    }
    
    Y[, c('hemo', 'hr', 'map', 'lactate')] = newY
    return(Y)
}

update_invKn = function( EIDs, par, par_index, Y){
    
    theta = exp(par[par_index$log_theta])
    
    invKn = foreach( i=EIDs, .packages='Matrix') %dopar% {
        
        n_i = sum(Y[,'EID']==as.numeric(i))
        
        diagonals = list( rep( (1 + exp(-2*theta)), n_i), rep( -exp(-theta), n_i-1))
        invK_i = bandSparse( n=n_i, k=0:1, diagonals=diagonals, symmetric=T)
        invK_i[1,1] = invK_i[n_i,n_i] = 1
        invK_i = invK_i * 2*theta /(1 - exp(-2*theta))
        
        return(invK_i)
        # Recall the form of Kn[[i]]
        #Kn[[i]] = matrix( 0, n_i, n_i)
        #for(k in 1:n_i){
        #	for(l in 1:n_i)  Kn[[i]][k,l] = .5 * exp(-theta *abs(k-l)) /theta
        #}
    }
    
    print(theta)
    print(as.matrix(invKn[[1]]))
    
    names(invKn) = EIDs
    
    return(invKn)
}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# load('Model_out/mcmc_out_interm_1_1it6.rda')
test_fnc(5);

# EIDs = names(test_post$Dn)

# n_cores = 4
# cl <- makeCluster(n_cores, outfile="")
# registerDoParallel(cl)

# s_time = Sys.time()
# value = log_post(EIDs, test_post$par, test_post$par_index, test_post$A, test_post$B,
#                          test_post$Y, test_post$z, test_post$Dn, test_post$Xn, test_post$invKn)
# e_time = Sys.time() - s_time; print(e_time)
# 
# s_time = Sys.time()
# Dn_R = update_Dn(EIDs, test_post$B, test_post$Y)
# e_time = Sys.time() - s_time; print(e_time)
# 
# s_time = Sys.time()
# alpha_R = update_alpha_i(EIDs, test_post$par, test_post$par_index, test_post$Y,
#                          test_post$Dn, test_post$Xn, test_post$invKn)
# e_time = Sys.time() - s_time; print(e_time)
# 
# s_time = Sys.time()
# alpha_tilde_R = update_alpha_tilde(EIDs, test_post$par, test_post$par_index, 
#                                    test_post$A, test_post$Y)
# e_time = Sys.time() - s_time; print(e_time)

# s_time = Sys.time()
# beta_up_r = update_beta_Upsilon_R(EIDs, test_post$par, test_post$par_index, 
#                                   test_post$A, test_post$Y, test_post$Dn, test_post$Xn,
#                                   test_post$invKn)
# e_time = Sys.time() - s_time; print(e_time)

# s_time = Sys.time()
# Y_R = update_Y_i(EIDs, test_post$par, test_post$par_index, test_post$A, test_post$Y,
#                  test_post$Dn, test_post$Xn, test_post$invKn, test_post$otype)
# e_time = Sys.time() - s_time; print(e_time)

# s_time = Sys.time()
# invKn_R = update_invKn(EIDs, test_post$par, test_post$par_index, test_post$Y)
# e_time = Sys.time() - s_time; print(e_time)

# stopCluster(cl)


# EIDs_c = as.numeric(names(test_post$Dn))

# s_time = Sys.time()
# value_cpp = log_post_cpp(EIDs_c, test_post$par, test_post$par_index, test_post$A, test_post$B,
#                          test_post$Y, test_post$z, test_post$Dn, test_post$Xn, test_post$invKn)
# e_time = Sys.time() - s_time; print(e_time)
# 
# s_time = Sys.time()
# Dn_C = update_Dn_cpp(EIDs_c, test_post$B, test_post$Y)
# e_time = Sys.time() - s_time; print(e_time)
# 
# s_time = Sys.time()
# alpha_C = update_alpha_i_cpp(EIDs_c, test_post$par, test_post$par_index, test_post$Y,
#                          test_post$Dn, test_post$Xn, test_post$invKn)
# e_time = Sys.time() - s_time; print(e_time)
# 
# s_time = Sys.time()
# alpha_tilde_C = update_alpha_tilde_cpp(EIDs_c, test_post$par, test_post$par_index, 
#                                        test_post$A, test_post$Y)
# e_time = Sys.time() - s_time; print(e_time)

# s_time = Sys.time()
# beta_up_r_cpp = update_beta_Upsilon_R_cpp(EIDs_c, test_post$par, test_post$par_index, 
#                                          test_post$A, test_post$Y, test_post$Dn, test_post$Xn,
#                                          test_post$invKn)
# e_time = Sys.time() - s_time; print(e_time)

# s_time = Sys.time()
# Y_C =  update_Y_i_cpp(EIDs_c, test_post$par, test_post$par_index, test_post$A, test_post$Y,
#                       test_post$Dn, test_post$Xn, test_post$invKn, test_post$otype)
# e_time = Sys.time() - s_time; print(e_time)

# s_time = Sys.time()
# invKn_cpp = update_invKn_cpp(EIDs_c, test_post$par, test_post$par_index, test_post$Y)
# e_time = Sys.time() - s_time; print(e_time)

# print(sum(is.na(Y_R)))
# print(sum(is.na(Y_C)))

# print(value)
# print(value_cpp)
# 
# print(dim(Dn_R[[1]]))
# print(dim(Dn_C[[1]]))
# print(sum(Dn_R[[1]] == Dn_C[[1]]))
# 
# print(alpha_R[[1]])
# print(alpha_C[[1]])
# 
# print(alpha_tilde_R[test_post$par_index$vec_alpha_tilde])
# print(alpha_tilde_C[test_post$par_index$vec_alpha_tilde])

# print(head(test_post$Y))
# print(head(Y_R))
