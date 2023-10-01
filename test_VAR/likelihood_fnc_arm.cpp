#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// Defining the Omega_List as a global variable when pre-compiling ----------
const arma::mat adj_mat = { {1, 1, 0},
                            {0, 1, 1},
                            {1, 1, 1} };

const arma::mat adj_mat_sub = { {1, 0, 0},
                                {0, 1, 0},
                                {1, 0, 1} };

arma::field<arma::field<arma::mat>> Omega_set(const arma::mat &G) {
  int N = G.n_cols; // dimension of adj matrix
  
  arma::field<arma::mat> c(N);
  arma::field<arma::mat> b(N, N);
  arma::field<arma::mat> a(N);
  
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      b(i, j) = arma::mat(1, 2, arma::fill::zeros);
    }
  }
  
  // a -------------------------------------------------------
  for(int i = 0; i < N; i++) {
    arma::mat a_i(1, 2, arma::fill::zeros);
    arma::uvec sec_elem = arma::find(G.row(i) == 1);
    
    for(int j = 0; j < sec_elem.n_elem; j++) {
      int sec_ind = sec_elem(j);
      arma::uvec third_elem = arma::find(G.row(sec_ind) == 1);
      
      for(int k = 0; k < third_elem.n_elem; k++) {
        int third_ind = third_elem(k);
        arma::mat temp(1,2);
        temp(0,0) = sec_ind; temp(0,1) = third_ind;
        a_i = arma::join_vert(a_i, temp+1);
      }
    }
    
    a_i = a_i.rows(1, a_i.n_rows-1);
    a(i) = a_i;
  }
  
  // b -------------------------------------------------------
  for(int i = 0; i < N; i++) {
    arma::uvec sec_elem = arma::find(G.row(i) == 1);
    
    for(int j = 0; j < sec_elem.n_elem; j++) {
      int sec_ind = sec_elem(j);
      arma::uvec third_elem = arma::find(G.row(sec_ind) == 1);
      
      for(int k = 0; k < third_elem.n_elem; k++) {
        int third_ind = third_elem(k);
        arma::uvec fourth_elem = arma::find(G.row(third_ind) == 1);
        
        for(int l = 0; l < fourth_elem.n_elem; l++) {
          int fourth_ind = fourth_elem(l);
          arma::mat temp(1,2);
          temp(0,0) = sec_ind; temp(0,1) = third_ind;
          b(i, fourth_ind) = arma::join_vert(b(i, fourth_ind), temp+1);
        }
      }
    }
  }
  
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      if (b(i,j).n_rows > 1) {
        b(i, j) = b(i, j).rows(1, b(i,j).n_rows - 1);
      } else{
        arma::mat temp(1,2); temp(0,0) = -1; temp(0,1) = -1;
        b(i,j) = temp;
      }
    }
  }
  
  // c -------------------------------------------------------
  for(int i = 0; i < N; i++) {
    arma::mat c_i(1, 2, arma::fill::zeros);
    arma::uvec sec_elem = arma::find(G.col(i) == 1);
    
    for(int j = 0; j < sec_elem.n_elem; j++) {
      int sec_ind = sec_elem(j);
      arma::uvec third_elem = arma::find(G.col(sec_ind) == 1);
      
      for(int k = 0; k < third_elem.n_elem; k++) {
        int third_ind = third_elem(k);
        arma::mat temp(1,2);
        temp(0,0) = third_ind; temp(0,1) = sec_ind;
        c_i = arma::join_vert(c_i, temp+1);
      }
    }
    
    c_i = c_i.rows(1, c_i.n_rows-1);
    c(i) = c_i;
  }
  
  arma::field<arma::field<arma::mat>> Omega_List(3);
  Omega_List(0) = c; Omega_List(1) = b; Omega_List(2) = a;
  
  return Omega_List;
}

const arma::field<arma::field<arma::mat>> Omega_List_GLOBAL = Omega_set(adj_mat);

const arma::field<arma::field<arma::mat>> Omega_List_GLOBAL_sub = Omega_set(adj_mat_sub);

// All other functions -------------------------------------------------------

arma::mat Omega_fun_cpp_new(const int k, const int n_i, const arma::vec &b_i,
                            const bool sub) {
  
  arma::mat Omega_set;
  // We are in a reduced dimension case
  if(sub) { 
    // b(k) is either 1, 2, 3, 4, 5 therefore subtract 1 for the index
    if (k == 1){
      // () -> () -> 1-5
      Omega_set = Omega_List_GLOBAL_sub(0)(b_i(2) - 1);
    } else if (k <= n_i - 2) {
      // 1-5 -> () -> () -> 1-5
      Omega_set = Omega_List_GLOBAL_sub(1)(b_i(k - 2) - 1, b_i(k + 1) - 1);
    } else if (k == n_i - 1) {
      // 1-5 -> () -> ()
      Omega_set = Omega_List_GLOBAL_sub(2)(b_i(n_i - 3) - 1);
    }
  } else {
    // b(k) is either 1, 2, 3, 4, 5 therefore subtract 1 for the index
    if (k == 1) {
      // () -> () -> 1-5
      Omega_set = Omega_List_GLOBAL(0)(b_i(2) - 1);
    } else if (k <= n_i - 2) {
      // 1-5 -> () -> () -> 1-5
      Omega_set = Omega_List_GLOBAL(1)(b_i(k - 2) - 1, b_i(k + 1) - 1);
    } else if (k == n_i - 1) {
      // 1-5 -> () -> ()
      Omega_set = Omega_List_GLOBAL(2)(b_i(n_i - 3) - 1);
    }
  }
  
  return Omega_set;
  
}

double log_f_i_cpp(const int i, const int ii, arma::vec t_pts, const arma::vec &par, 
                   const arma::field<arma::uvec> &par_index, const arma::vec &A, const arma::vec &B, 
                   const arma::mat &Y, const arma::mat &z, const arma::field<arma::mat> &Dn, 
                   const arma::field<arma::mat> &Xn, const arma::mat &Dn_omega, const arma::vec &W) {
  
  // par_index KEY: (0) beta, (1) alpha_tilde, (2) upsilon, (3) vec_A, (4) R, (5) zeta, 
  //                (6) init, (7) R_bleed, (8) omega_tilde, (9) vec_upsilon_omega
  
  // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
  // "i" is the numeric EID number
  // "ii" is the index of the EID
  double in_value = 0;
  
  arma::vec eids = Y.col(0);

  arma::vec vec_A_total = par.elem(par_index(3) - 1);
  arma::vec vec_A_scale = {exp(vec_A_total(0)) / (1 + exp(vec_A_total(0))),
                           exp(vec_A_total(1)) / (1 + exp(vec_A_total(1))),
                           exp(vec_A_total(2)) / (1 + exp(vec_A_total(2))),
                           exp(vec_A_total(3)) / (1 + exp(vec_A_total(3)))};
  arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 1);

  arma::vec vec_R_base = par.elem(par_index(4) - 1);
  arma::vec vec_R_bleed = par.elem(par_index(7) - 1);
  arma::mat R_base = arma::reshape(vec_R_base, 4, 4);
  arma::mat R_bleed = arma::reshape(vec_R_bleed, 4, 4);
  
  arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
  arma::mat zeta = arma::reshape(vec_zeta_content, 2, 4); // THREE STATE
  
  // The time-homogeneous probability transition matrix
  arma::uvec sub_ind = arma::find(eids == i);
  arma::mat z_i = z.rows(sub_ind.min(), sub_ind.max());
  int n_i = z_i.n_rows;
  
  // The state transition likelihood component for current iterate of b_i
  arma::vec b_i = B;
  arma::vec vec_init_content = par.elem(par_index(6) - 1);
  arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1))}; // THREE STATE
  arma::vec P_init = init_logit / arma::accu(init_logit); 
  
  // Subsetting the data
  arma::mat Y_temp = Y.rows(sub_ind);
  arma::mat Y_i = Y_temp.cols(1, 4);
  Y_i = Y_i.t();
  
  arma::field<arma::mat> Dn_alpha_full = Dn;
  arma::field<arma::mat> Xn_full = Xn;
  arma::vec vec_alpha_ii = A;
  arma::mat vec_beta = par.elem(par_index(0) - 1);
  
  // Full likelihood evaluation is not needed for updating pairs of b_i components
  if (any(t_pts == -1)) { t_pts = arma::linspace(1, n_i, n_i);}
  
  for(int w=0; w < t_pts.n_elem;++w){
    int k = t_pts(w);
    if(k==1){
        // State space component
        int b_0 = b_i(0);
        
        arma::vec vec_A = A_all_state.col(0);
        arma::mat R;
        if(b_0 == 1) {
            R = R_base;
        } else {
            R = R_bleed;
        }
        
        arma::mat Gamma     = {{R(0,0) / (1 - vec_A(0) * vec_A(0)), 
                                R(0,1) / (1 - vec_A(0) * vec_A(1)), 
                                R(0,2) / (1 - vec_A(0) * vec_A(2)), 
                                R(0,3) / (1 - vec_A(0) * vec_A(3))},
                                {R(1,0) / (1 - vec_A(1) * vec_A(0)), 
                                 R(1,1) / (1 - vec_A(1) * vec_A(1)), 
                                 R(1,2) / (1 - vec_A(1) * vec_A(2)), 
                                 R(1,3) / (1 - vec_A(0) * vec_A(3))},
                                 {R(2,0) / (1 - vec_A(2) * vec_A(0)), 
                                  R(2,1) / (1 - vec_A(2) * vec_A(1)), 
                                  R(2,2) / (1 - vec_A(2) * vec_A(2)), 
                                  R(2,3) / (1 - vec_A(0) * vec_A(3))},
                                  {R(3,0) / (1 - vec_A(3) * vec_A(0)), 
                                   R(3,1) / (1 - vec_A(3) * vec_A(1)), 
                                   R(3,2) / (1 - vec_A(3) * vec_A(2)), 
                                   R(3,3) / (1 - vec_A(0) * vec_A(3))}};
      
        arma::vec y_1 = Y_i.col(0);
        arma::vec nu_1 = Dn_alpha_full(0) * vec_alpha_ii + Xn_full(0) * vec_beta;
        arma::vec log_y_pdf = dmvnorm(y_1.t(), nu_1, Gamma, true);
        
        in_value = in_value + log(P_init(b_0 - 1)) + arma::as_scalar(log_y_pdf);
    }
    else{
        // State space component
        double q1_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(0));
        double q1 = exp(q1_sub);
        
        double q2_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(1));
        double q2 = exp(q2_sub);
        
        double q3_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(2));
        double q3 = exp(q3_sub);
        
        double q4_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(3));
        double q4 = exp(q4_sub);
        
        arma::mat Q = {{  1,  q1,  0},
                     {  0,   1, q2},
                     { q3,  q4,  1}}; // THREE STATE
        
        arma::vec q_row_sums = arma::sum(Q, 1);
        arma::mat P_i = Q.each_col() / q_row_sums;
        int b_k_1 = b_i(k-2);
        int b_k = b_i(k-1);
        
        arma::mat R;
        if(b_k == 1) {
            R = R_base;
        } else {
            R = R_bleed;
        }
        
        arma::vec vec_A = A_all_state.col(0);
        arma::mat A_1 = arma::diagmat(vec_A);
        
        arma::vec y_k_1 = Y_i.col(k-2);
        arma::vec y_k = Y_i.col(k-1);
        arma::vec nu_k_1 = Dn_alpha_full(k-2) * vec_alpha_ii + Xn_full(k-2) * vec_beta;
        arma::vec nu_k = Dn_alpha_full(k-1) * vec_alpha_ii + Xn_full(k-1) * vec_beta;
        
        arma::vec mean_k = nu_k + A_1 * (y_k_1 - nu_k_1);
        
        arma::vec log_y_k_pdf = dmvnorm(y_k.t(), mean_k, R, true);
        
        in_value = in_value + log(P_i( b_k_1 - 1, b_k - 1)) + arma::as_scalar(log_y_k_pdf);
    }
  }

  return in_value;
}

double log_f_i_cpp_total(const arma::vec &EIDs, arma::vec t_pts, const arma::vec &par, const arma::field<arma::uvec> &par_index, 
                         const arma::field <arma::vec> &A, const arma::field <arma::vec> &B, 
                         const arma::mat &Y, const arma::mat &z, const arma::field<arma::field<arma::mat>> &Dn, 
                         const arma::field <arma::field<arma::mat>> &Xn, const arma::field <arma::mat> &Dn_omega, 
                         const arma::field <arma::vec> &W) {

    // par_index KEY: (0) beta, (1) alpha_tilde, (2) upsilon, (3) vec_A, (4) R, (5) zeta, 
    //                (6) init, (7) R_bleed, (8) omega_tilde, (9) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);
      
    omp_set_num_threads(25);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        in_vals(ii) = log_f_i_cpp(i, ii, t_pts, par, par_index, A(ii), B(ii), Y, z, Dn(ii), Xn(ii), Dn_omega(ii), W(ii));
    }
    
    double in_value = arma::accu(in_vals);
    
    return in_value;
}

// [[Rcpp::export]]
double log_post_cpp(const arma::vec &EIDs, const arma::vec &par, const arma::field<arma::uvec> &par_index,
                    const arma::field<arma::vec> &A, const arma::field<arma::vec> &B,
                    const arma::mat &Y, const arma::mat &z, const arma::field<arma::field<arma::mat>> &Dn,
                    const arma::field<arma::field<arma::mat>> &Xn, const arma::field <arma::mat> &Dn_omega, 
                    const arma::field<arma::vec> &W) {

  // par_index KEY: (0) beta, (1) alpha_tilde, (2) upsilon, (3) vec_A, (4) R, (5) zeta, 
  //                (6) init, (7) R_bleed, (8) omega_tilde, (9) vec_upsilon_omega
  // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
  // "i" is the numeric EID number
  // "ii" is the index of the EID

  // Compute the likelihood ----------------------------------------------------
  double value;
  arma::vec t_pts = {-1};
  value = log_f_i_cpp_total(EIDs, t_pts, par, par_index, A, B, Y, z, Dn, Xn, Dn_omega, W);
  // ---------------------------------------------------------------------------

  // Compute prior densities of all necessary model parameters -----------------
  // Zeta priors ---------------------------------------------------------------
  arma::uvec vec_zeta_ind = par_index(5);
  arma::vec vec_zeta_content = par.elem(vec_zeta_ind - 1);
  arma::mat zeta = arma::reshape(vec_zeta_content, 2, 4); // THREE STATE
  int m = zeta.n_rows;

  arma::vec vec_zeta_mean(8, arma::fill::zeros);
  arma::vec scalar_1 = {10, 10, 10, 10, 10, 10, 10, 10}; // UNINFORMATIVE
  arma::mat zeta_sd = arma::diagmat(scalar_1);

  arma::vec prior_zeta = dmvnorm(vec_zeta_content.t(), vec_zeta_mean, zeta_sd, true);
  double prior_zeta_val = arma::as_scalar(prior_zeta);

  // Initial Probabilities priors ----------------------------------------------
  arma::uvec vec_init_ind = par_index(6);
  arma::vec vec_init_content = par.elem(vec_init_ind - 1);
  arma::vec vec_init_mean = {0, 0}; // THREE STATE
  arma::vec scalar_2 = {10,10}; // THREE STATE
  arma::mat init_sd = arma::diagmat(scalar_2);

  arma::vec prior_init = dmvnorm(vec_init_content.t(), vec_init_mean, init_sd, true);
  double prior_init_val = arma::as_scalar(prior_init);
  
  // A_1 priors ----------------------------------------------------------------
  arma::vec vec_A1_content = par.elem(par_index(3) - 1);
  
  arma::vec vec_A1_mean = {0, 0, 0, 0};
  arma::vec A1_scalar   = {2, 2, 2, 2}; 
  arma::mat A1_sd = arma::diagmat(A1_scalar);
  
  arma::vec prior_A1 = dmvnorm(vec_A1_content.t(), vec_A1_mean, A1_sd, true);
  double prior_A1_val = arma::as_scalar(prior_A1);
  
  // R priors ------------------------------------------------------------------
  arma::vec vec_R_base = par.elem(par_index(4) - 1);
  arma::vec vec_R_bleed = par.elem(par_index(7) - 1);
  arma::mat R_base = arma::reshape(vec_R_base, 4, 4);
  arma::mat R_bleed = arma::reshape(vec_R_bleed, 4, 4);
  
  arma::mat Psi_R(4,4,arma::fill::eye);
  int nu_R = 6;
  
  double prior_R_base_val = diwish(R_base, nu_R, Psi_R, true);
  double prior_R_bleed_val = diwish(R_bleed, nu_R, Psi_R, true);
  
  // Omega priors -------------------------------------------------------------
  // arma::vec vec_omega_content = par.elem(par_index(8) - 1);
  // arma::vec omega_mean(8, arma::fill::zeros);
  // arma::vec diag_omega_sd(8, arma::fill::ones);
  // diag_omega_sd = 20 * diag_omega_sd;
  // arma::mat omega_sd = arma::diagmat(diag_omega_sd);
  // 
  // arma::vec prior_omega = dmvnorm(vec_omega_content.t(), omega_mean, omega_sd, true);
  // double prior_omega_val = arma::as_scalar(prior_omega);
  // ---------------------------------------------------------------------------
  
  
  value = value + prior_A1_val + prior_R_base_val + prior_R_bleed_val + prior_zeta_val + prior_init_val;
  return value;
}

// [[Rcpp::export]]
Rcpp::List update_b_i_cpp(const int t, const arma::vec EIDs, const arma::vec par, const arma::field<arma::uvec> par_index, 
                          const arma::field <arma::vec> A, arma::field <arma::vec> B, 
                          const arma::mat Y, const arma::mat z, arma::field<arma::field<arma::mat>> Dn, 
                          const arma::field <arma::field<arma::mat>> Xn, const arma::field <arma::mat> Dn_omega, 
                          const arma::field <arma::vec> W, arma::mat &l1, arma::mat &l2) {

  // par_index KEY: (0) beta, (1) alpha_tilde, (2) upsilon, (3) vec_A, (4) R, (5) zeta,
  //                (6) init, (7) R_bleed, (8) omega_tilde, (9) vec_upsilon_omega
  // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
  // "i" is the numeric EID number
  // "ii" is the index of the EID

  arma::vec eids = Y.col(0); 
  arma::vec rbc_rule_vec = Y.col(5);
  arma::vec clinic_rule_vec = Y.col(6); 
  
  arma::field<arma::vec> B_return(EIDs.n_elem);
  arma::field<arma::field<arma::mat>> Dn_return(EIDs.n_elem);
  
  omp_set_num_threads(t) ;
  # pragma omp parallel for
  for (int ii = 0; ii < EIDs.n_elem; ii++) {
    int i = EIDs(ii);
    arma::uvec sub_ind = arma::find(eids == i);
    
    int n_i = sub_ind.n_elem;
    
    int rbc_rule = rbc_rule_vec(sub_ind.min());
    int clinic_rule = clinic_rule_vec(sub_ind.min());
    
    // Subsetting fields
    arma::vec B_temp = B(ii);
    arma::mat Dn_omega_temp = Dn_omega(ii);
    arma::vec A_temp = A(ii);
    arma::vec W_temp = W(ii);
    
    arma::field<arma::mat> Dn_temp = Dn(ii);
    arma::field<arma::mat> Xn_temp = Xn(ii);
    
    // Subsetting the remaining data
    arma::mat Y_temp = Y.rows(sub_ind);
    arma::mat z_temp = z.rows(sub_ind);
    
    // keep the last state at state 1 so we only iterate to n_i-2 (CHANGED*****)
    for (int k = 0; k < n_i - 1; k++) {
      
      arma::vec t_pts;
      if (k == n_i - 2) {
          t_pts = arma::linspace(k+1, k+2, 2);
      } else {
          t_pts = arma::linspace(k+1, k+3, 3);
      }
      
      arma::vec pr_B = B_temp;
      arma::field<arma::mat> pr_Dn = Dn_temp;
      
      // Sample and update the two neighboring states
      arma::mat Omega_set;
      if (clinic_rule >= 0) {
        Omega_set = Omega_fun_cpp_new(k + 1, n_i, B_temp, false);
      } else {
        Omega_set = Omega_fun_cpp_new(k + 1, n_i, B_temp, true);
      }

      int sampled_index = arma::randi(arma::distr_param(1, Omega_set.n_rows));
      
      pr_B.rows(k, k+1) = Omega_set.row(sampled_index-1).t();
      
      // DEBUG ----------------------------------------------------------------
      // Rows: likelihood b4, likelihood after, p1, p2, accept
      if(i == -1) {
          l1(arma::span(2,3), k) = Omega_set.row(sampled_index-1).t();
      }
      if(i == -2) {
          l2(arma::span(2,3), k) = Omega_set.row(sampled_index-1).t();
      }
      // ----------------------------------------------------------------------
      
      // Adding clinical review
      bool valid_prop = false;
      
      if(clinic_rule >= 0) {
        bool b_i_rule = arma::any(arma::vectorise(pr_B)==2);
        if (clinic_rule == 1) {
          if(b_i_rule) {valid_prop = true;}
        } else {
          if (rbc_rule == 0 || (rbc_rule == 1 && b_i_rule)) {valid_prop = true;}
        }
      } else {
        valid_prop = true; // evaluate likelihood anyways because S1&S3
      }
      
      if(valid_prop) {
        double log_target_prev = log_f_i_cpp(i, ii, t_pts, par, par_index,A_temp,
                                             B_temp,Y_temp,z_temp,Dn_temp,Xn_temp,
                                             Dn_omega_temp, W_temp);
        
        arma::vec twos(pr_B.n_elem, arma::fill::zeros);
        arma::vec threes = twos; // THREE STATE
        
        twos.elem(arma::find(pr_B == 2)) += 1;
        threes.elem(arma::find(pr_B == 3)) += 1; // THREE STATE
        
        arma::vec ones(pr_B.n_elem, arma::fill::ones);
        
        arma::mat bigB = arma::join_rows(ones, arma::cumsum(twos));
        bigB = arma::join_rows(bigB, arma::cumsum(threes)); // THREE STATE
        
        arma::mat I = arma::eye(4,4);
        for(int jj = 0; jj < n_i; jj++) {
            pr_Dn(jj) = arma::kron(I, bigB.row(jj));
        }
        
        double log_target = log_f_i_cpp(i,ii,t_pts,par,par_index,A_temp,
                                        pr_B,Y_temp,z_temp,pr_Dn,Xn_temp,
                                        Dn_omega_temp, W_temp);
        
        // DEBUG ----------------------------------------------------------------
        // Rows: likelihood b4, likelihood after, p1, p2, accept
        if(i == -1) {
            l1(0, k) = log_target_prev;
            l1(1, k) = log_target;
            l1(4, k) = 0;
            l1(5, k) = B_temp(k);
            l1(6, k) = B_temp(k+1);
        }
        if(i == -2) {
            l2(0, k) = log_target_prev;
            l2(1, k) = log_target;
            l2(4, k) = 0;
            l2(5, k) = B_temp(k);
            l2(6, k) = B_temp(k + 1);
        }
        // ----------------------------------------------------------------------
        
        
        // Note that the proposal probs cancel in the MH ratio
        double diff_check = log_target - log_target_prev;
        double min_log = log(arma::randu(arma::distr_param(0,1)));
        if(diff_check > min_log){
          B_temp = pr_B;
          Dn_temp = pr_Dn;
          
          // DEBUG ----------------------------------------------------------------
          // Rows: likelihood b4, likelihood after, p1, p2, accept
          if(i == -1)  { l1(4, k) = 1;}
          if(i == -2) { l2(4, k) = 1;}
          // ----------------------------------------------------------------------
        }
      }
    }
    B_return(ii) = B_temp;
    Dn_return(ii) = Dn_temp;
  }
  List B_Dn = List::create(B_return, Dn_return, l1, l2);
  
  return B_Dn;
}

// [[Rcpp::export]]
Rcpp::List update_Dn_Xn_cpp( const arma::vec EIDs, arma::field <arma::vec> B, 
                             const arma::mat Y, const arma::vec par, 
                             const arma::field<arma::uvec> par_index,
                             const arma::vec x) {

  // par_index KEY: (0) beta, (1) alpha_tilde, (2) upsilon, (3) vec_A, (4) R, (5) zeta,
  //                (6) init, (7) R_bleed, (8) omega_tilde, (9) vec_upsilon_omega
  // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
  // "i" is the numeric EID number
  // "ii" is the index of the EID

  arma::vec eids = Y.col(0);
  arma::field<arma::field<arma::mat>> Dn(EIDs.n_elem);
  arma::field<arma::field<arma::mat>> Xn(EIDs.n_elem);
  
  omp_set_num_threads(8) ;
  # pragma omp parallel for
  for (int ii = 0; ii < EIDs.n_elem; ii++) {
    int i = EIDs(ii);
    arma::uvec sub_ind = arma::find(eids == i);
    int n_i = sub_ind.n_elem;
    
    arma::field<arma::mat> Dn_i(n_i);
    arma::field<arma::mat> Xn_i(n_i);
    
    arma::vec b_i = B(ii);
    arma::vec x_i = x.elem(sub_ind);
    
    arma::vec twos(b_i.n_elem, arma::fill::zeros);
    arma::vec threes = twos; // THREE STATE
    
    twos.elem(arma::find(b_i == 2)) += 1;
    threes.elem(arma::find(b_i == 3)) += 1; // THREE STATE
    
    arma::vec ones(b_i.n_elem, arma::fill::ones);
    
    arma::mat bigB = arma::join_rows(ones, arma::cumsum(twos));
    bigB = arma::join_rows(bigB, arma::cumsum(threes)); // THREE STATE
    
    arma::mat I = arma::eye(4,4);

    for(int jj = 0; jj < n_i; jj++) {
        Dn_i(jj) = arma::kron(I, bigB.row(jj));
        Xn_i(jj) = x_i(jj) * I;
    }
    
    Dn(ii) = Dn_i;
    Xn(ii) = Xn_i;
  }
  
  List Dn_Xn = List::create(Dn, Xn);
  
  return Dn_Xn;
}

// [[Rcpp::export]]
arma::field <arma::vec> update_alpha_i_cpp( const arma::vec EIDs, const arma::vec par, 
                                            const arma::field<arma::uvec> par_index,
                                            const arma::mat Y, arma::field<arma::field<arma::mat>> Dn, 
                                            const arma::field <arma::field<arma::mat>> Xn,
                                            const arma::field <arma::mat> Dn_omega, 
                                            const arma::field <arma::vec> W, 
                                            arma::field <arma::vec> B){

    // par_index KEY: (0) beta, (1) alpha_tilde, (2) upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) R_bleed, (8) omega_tilde, (9) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID

    arma::vec eids = Y.col(0);

    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = {exp(vec_A_total(0)) / (1 + exp(vec_A_total(0))),
                             exp(vec_A_total(1)) / (1 + exp(vec_A_total(1))),
                             exp(vec_A_total(2)) / (1 + exp(vec_A_total(2))),
                             exp(vec_A_total(3)) / (1 + exp(vec_A_total(3)))};
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 1);
    
    arma::vec vec_R_base = par.elem(par_index(4) - 1);
    arma::vec vec_R_bleed = par.elem(par_index(7) - 1);
    arma::mat R_base = arma::reshape(vec_R_base, 4, 4);
    arma::mat R_bleed = arma::reshape(vec_R_bleed, 4, 4);
    arma::mat invR_base = arma::inv_sympd(R_base);
    arma::mat invR_bleed = arma::inv_sympd(R_bleed);

    arma::vec vec_alpha_tilde = par.elem(par_index(1) - 1);

    arma::vec vec_beta = par.elem(par_index(0) - 1);

    arma::vec Upsilon_vec = par.elem(par_index(2) - 1);
    arma::mat Upsilon = arma::reshape(Upsilon_vec, 12, 12); // THREE STATE
    arma::mat inv_Upsilon = arma::inv_sympd(Upsilon);
    
    arma::field<arma::vec> A(EIDs.n_elem);
  
    omp_set_num_threads(16) ;
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
      int i = EIDs(ii);
    
      arma::uvec sub_ind = arma::find(eids == i);
      int n_i = sub_ind.n_elem;
      
      arma::vec b_i = B(ii);
    
      arma::mat Y_temp = Y.rows(sub_ind);
      arma::mat Y_i = Y_temp.cols(1, 4);
      Y_i = Y_i.t();
      
      arma::field<arma::mat> Dn_alpha_full = Dn(ii);
      arma::field<arma::mat> Xn_full = Xn(ii);
      
      arma::mat interm_W(12,12,arma::fill::zeros);
      arma::vec interm_V(12,arma::fill::zeros);
      for(int k=1; k < n_i; k++) {
          arma::mat D_k = Dn_alpha_full(k);
          arma::mat D_k_1 = Dn_alpha_full(k-1);
          
          arma::mat invR;
          if(b_i(k) == 1) {
              invR = invR_base;
          } else {
              invR = invR_bleed;
          }
          
          arma::vec vec_A_k = A_all_state.col(0);
          arma::mat A_1_k = arma::diagmat(vec_A_k);
          
          arma::mat diff_hold_d = D_k - A_1_k * D_k_1;
          interm_W = interm_W + diff_hold_d.t() * invR * diff_hold_d;
          
          arma::mat X_k = Xn_full(k);
          arma::mat X_k_1 = Xn_full(k-1);
          arma::mat diff_hold_x = X_k - A_1_k * X_k_1;
          
          arma::vec m_k = Y_i.col(k) - A_1_k * Y_i.col(k-1) - diff_hold_x * vec_beta;
          interm_V = interm_V + diff_hold_d.t() * invR * m_k;
      }
      
      arma::mat R;
      if(b_i(0) == 1) {
          R = R_base;
      } else {
          R = R_bleed;
      }
      
      arma::vec vec_A = A_all_state.col(0);
      arma::mat Gamma     = {{R(0,0) / (1 - vec_A(0) * vec_A(0)), 
                              R(0,1) / (1 - vec_A(0) * vec_A(1)), 
                              R(0,2) / (1 - vec_A(0) * vec_A(2)), 
                              R(0,3) / (1 - vec_A(0) * vec_A(3))},
                             {R(1,0) / (1 - vec_A(1) * vec_A(0)), 
                              R(1,1) / (1 - vec_A(1) * vec_A(1)), 
                              R(1,2) / (1 - vec_A(1) * vec_A(2)), 
                              R(1,3) / (1 - vec_A(0) * vec_A(3))},
                             {R(2,0) / (1 - vec_A(2) * vec_A(0)), 
                              R(2,1) / (1 - vec_A(2) * vec_A(1)), 
                              R(2,2) / (1 - vec_A(2) * vec_A(2)), 
                              R(2,3) / (1 - vec_A(0) * vec_A(3))},
                             {R(3,0) / (1 - vec_A(3) * vec_A(0)), 
                              R(3,1) / (1 - vec_A(3) * vec_A(1)), 
                              R(3,2) / (1 - vec_A(3) * vec_A(2)), 
                              R(3,3) / (1 - vec_A(0) * vec_A(3))}};
    
      arma::mat inv_Gamma = arma::inv_sympd(Gamma);
      
      arma::mat W_i_inv = inv_Upsilon + Dn_alpha_full(0).t() * inv_Gamma * Dn_alpha_full(0) + interm_W;
    
      arma::mat W_i = inv(W_i_inv);
      
      arma::mat V_i = inv_Upsilon * vec_alpha_tilde + 
          Dn_alpha_full(0).t() * inv_Gamma * (Y_i.col(0) - Xn_full(0) * vec_beta) + interm_V;
    
      arma::vec mu = W_i * V_i;
    
      arma::mat alpha_i = rmvnorm(1, mu, W_i);
    
      arma::vec vec_alpha_i = alpha_i.t();
    
      A(ii) = vec_alpha_i;
    }
    
    return A;
  
}

// [[Rcpp::export]]
arma::field <arma::vec> update_omega_i_cpp( const arma::vec EIDs, const arma::vec par, 
                                            const arma::field<arma::uvec> par_index,
                                            const arma::mat Y, arma::field <arma::mat> Dn, 
                                            const arma::field <arma::mat> Xn, const arma::field <arma::sp_mat> invKn,
                                            const arma::field <arma::mat> Dn_omega, 
                                            const arma::field <arma::vec> A){

  // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
  //                (6) init, (7) log_lambda, (8) omega_tilde, (9) vec_upsilon_omega
  // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
  // "i" is the numeric EID number
  // "ii" is the index of the EID

  arma::vec eids = Y.col(0);

  arma::uvec vec_R_ind = par_index(4);
  arma::vec vec_R_content = par.elem(vec_R_ind - 1);
  arma::mat R = arma::reshape(vec_R_content, 4, 4);

  arma::mat invR = inv(R);

  arma::vec vec_omega_tilde = par.elem(par_index(8) - 1);

  arma::uvec vec_beta_ind = par_index(0);
  arma::vec vec_beta = par.elem(vec_beta_ind - 1);

  arma::vec vec_upsilon_omega = par.elem(par_index(9) - 1);
  arma::mat Upsilon = arma::reshape(vec_upsilon_omega, 8, 8);
  arma::mat inv_Upsilon = arma::inv_sympd(Upsilon);

  arma::field<arma::vec> W(EIDs.n_elem);

  // omp_set_num_threads(16) ;
  // # pragma omp parallel for
  for (int ii = 0; ii < EIDs.n_elem; ii++)
  {
    int i = EIDs(ii);

    arma::uvec sub_ind = arma::find(eids == i);

    arma::mat Y_temp = Y.rows(sub_ind);
    arma::mat Y_i = Y_temp.cols(1, 4);
    arma::vec vecY_i = arma::vectorise(Y_i);

    arma::mat Dn_ii = Dn(ii);
    arma::mat Dn_omega_ii = Dn_omega(ii);
    arma::mat Xn_ii = Xn(ii);
    arma::sp_mat invKn_ii = invKn(ii);

    arma::sp_mat inv_R_fill = arma::sp_mat(invR);
    arma::sp_mat precision = arma::kron(inv_R_fill, invKn_ii);

    arma::mat hold = Dn_omega_ii.t() * precision;

    // Rcpp::Rcout << size(vecY_i) << std::endl;
    // Rcpp::Rcout << size(Dn_ii) << std::endl;
    // Rcpp::Rcout << size(A(ii)) << std::endl;
    // Rcpp::Rcout << size(Xn_ii) << std::endl;
    // Rcpp::Rcout << size(vec_beta) << std::endl;
    // Rcpp::Rcout << size(inv_Upsilon) << std::endl;
    // Rcpp::Rcout << size(vec_omega_tilde) << std::endl;
    // arma::mat temp_hold = hold*(vecY_i - Dn_ii*A(ii) - Xn_ii*vec_beta);

    arma::mat V_i = hold * (vecY_i - Dn_ii * A(ii) - Xn_ii * vec_beta) + inv_Upsilon * vec_omega_tilde;

    // Rcpp::Rcout << "issue2" << std::endl;
    arma::mat inv_W_i = hold * Dn_omega_ii + inv_Upsilon;
    // Rcpp::Rcout << "issue3" << std::endl;
    arma::mat W_i = inv(inv_W_i);

    arma::vec mu = W_i * V_i;

    arma::mat omega_i = rmvnorm(1, mu, W_i);

    arma::vec vec_omega_i = omega_i.t();

    W(ii) = vec_omega_i;
    }
    
    return W;
    
}

// [[Rcpp::export]]
arma::vec update_alpha_tilde_cpp( const arma::vec EIDs, arma::vec par, 
                                  const arma::field<arma::uvec> par_index,
                                  const arma::field <arma::vec> A, const arma::mat Y){

    // par_index KEY: (0) beta, (1) alpha_tilde, (2) upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) R_bleed, (8) omega_tilde, (9) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID

    // The prior mean for vec_alpha_tilde
    arma::vec vec_alpha_tilde_0 = {9.57729783, -1, 0.1,
                                   88.69780576, 5.04150472, -4,
                                   79.74903940, -5.04150472, 4,
                                   5.2113319, 0.5360813, -0.6866748};
    
    // The prior PRECISION matrix for vec_alpha_tilde
    // arma::vec inv_Sigma_alpha_diag = {1, 1, 1, 0.0025, 0.01, 0.01, 0.0025, 0.01, 0.01, 1, 1, 1};
    arma::vec inv_Sigma_alpha_diag = {0.1, 0.3, 0.5, 0.05, 0.1, 0.1,
                                      0.05, 0.1, 0.1, 0.1, 0.3, 0.5};
    
    arma::mat inv_Sigma_alpha = arma::diagmat(inv_Sigma_alpha_diag);
    
    arma::vec Upsilon_vec = par.elem(par_index(2) - 1);
    arma::mat Upsilon = arma::reshape(Upsilon_vec, 12, 12); // THREE STATE
    arma::mat inv_Upsilon = arma::inv_sympd(Upsilon);
    
    int length_EIDs = EIDs.n_elem;
    
    arma::mat inv_U = inv_Upsilon * length_EIDs + inv_Sigma_alpha;
    arma::mat U = inv(inv_U);
    
    arma::mat total_alpha = A(0);
    for (int i = 1; i < A.n_elem; i++){
        total_alpha = arma::join_horiz(total_alpha, A(i));
    }
  
    arma::mat sum_alpha = arma::sum(total_alpha, 1);
    
    arma::mat hold = inv_Sigma_alpha * vec_alpha_tilde_0 + inv_Upsilon * sum_alpha;
    
    arma::vec mu = U * hold;
    
    arma::uvec vec_alpha_tilde_ind = par_index(1);
    par.elem(vec_alpha_tilde_ind - 1) = rmvnorm(1, mu, U);
    
    return par;
}

// [[Rcpp::export]]
arma::vec update_omega_tilde_cpp( const arma::vec EIDs, arma::vec par, 
                                  const arma::field<arma::uvec> par_index,
                                  const arma::field <arma::vec> W, const arma::mat Y){

  // par_index KEY: (0) beta, (1) alpha_tilde, (2) upsilon, (3) vec_A, (4) R, (5) zeta, 
  //                (6) init, (7) R_bleed, (8) omega_tilde, (9) vec_upsilon_omega
  // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
  // "i" is the numeric EID number
  // "ii" is the index of the EID

  // The prior mean for vec_omega_tilde
  arma::vec vec_omega_tilde_0 = {4, -4, -4, 4, 4, -4, -4, 4};

  // The prior PRECISION matrix for vec_omega_tilde
  arma::vec inv_Sigma_omega_diag = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};

  arma::mat inv_Sigma_omega = arma::diagmat(inv_Sigma_omega_diag);

  arma::vec vec_upsilon_omega = par.elem(par_index(9) - 1);
  arma::mat Upsilon = arma::reshape(vec_upsilon_omega, 8, 8);
  arma::mat inv_Upsilon = arma::inv_sympd(Upsilon);

  int length_EIDs = EIDs.n_elem;

  arma::mat inv_U = inv_Upsilon * length_EIDs + inv_Sigma_omega;
  arma::mat U = inv(inv_U);

  arma::mat total_omega = W(0);
  for (int i = 1; i < W.n_elem; i++)
  {
    total_omega = arma::join_horiz(total_omega, W(i));
    }
    
    arma::mat sum_omega = arma::sum(total_omega, 1);
    
    arma::mat hold = inv_Sigma_omega * vec_omega_tilde_0 + inv_Upsilon * sum_omega;
    
    arma::vec mu = U * hold;
    
    arma::uvec vec_omega_tilde_ind = par_index(8);
    par.elem(vec_omega_tilde_ind - 1) = rmvnorm(1, mu, U);
    
    return par;
}

// [[Rcpp::export]]
arma::vec update_beta_Upsilon_R_cpp( const arma::vec EIDs, arma::vec par, 
                                     const arma::field<arma::uvec> par_index,
                                     const arma::field <arma::vec> A, const arma::mat Y,
                                     arma::field<arma::field<arma::mat>> Dn, 
                                     const arma::field <arma::field<arma::mat>> Xn, 
                                     const arma::field <arma::mat> Dn_omega, 
                                     const arma::field <arma::vec> W, arma::field <arma::vec> B) {
    // Conjugate updates for beta and sigma_upsilon
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) R_bleed, (8) omega_tilde, (9) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID

    // The prior info for vec_beta --------------------------------------------
    arma::uvec vec_beta_ind = par_index(0);
    arma::mat vec_beta_0(vec_beta_ind.n_elem, 1, arma::fill::zeros);

    arma::vec scalar_mult(vec_beta_ind.n_elem, arma::fill::ones);
    scalar_mult.fill(0.01);
    arma::mat inv_Sigma_beta = arma::diagmat(scalar_mult);
    
    arma::vec vec_R_base = par.elem(par_index(4) - 1);
    arma::vec vec_R_bleed = par.elem(par_index(7) - 1);
    arma::mat R_base = arma::reshape(vec_R_base, 4, 4);
    arma::mat R_bleed = arma::reshape(vec_R_bleed, 4, 4);
    arma::mat invR_base = arma::inv_sympd(R_base);
    arma::mat invR_bleed = arma::inv_sympd(R_bleed);

    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = {exp(vec_A_total(0)) / (1 + exp(vec_A_total(0))),
                             exp(vec_A_total(1)) / (1 + exp(vec_A_total(1))),
                             exp(vec_A_total(2)) / (1 + exp(vec_A_total(2))),
                             exp(vec_A_total(3)) / (1 + exp(vec_A_total(3)))};
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 1);
    
    arma::vec vec_beta = par.elem(vec_beta_ind - 1);

    // The prior info for sigma_upsilon ----------------------------------------
    int nu_Upsilon = 14;
    arma::vec scalar_mult2(12, arma::fill::ones);
    arma::mat Psi_Upsilon = arma::diagmat(scalar_mult2);

    arma::uvec vec_sigma_upsilon_ind = par_index(2);

    // The prior info for Upsilon_omega ----------------------------------------
    // int nu_upsilon_omega = 10;
    // arma::vec upsilon_omega_scalar(8, arma::fill::ones);
    // arma::mat Psi_omega = arma::diagmat(upsilon_omega_scalar);
    // arma::uvec vec_upsilon_omega_ind = par_index(9);
    // arma::vec vec_omega_tilde = par.elem(par_index(8) - 1);


    arma::uvec vec_alpha_tilde_ind = par_index(1);
    arma::vec vec_alpha_tilde = par.elem(vec_alpha_tilde_ind - 1);

    arma::vec eids = Y.col(0);

    arma::field<arma::mat> in_V_main(EIDs.n_elem);
    arma::field<arma::mat> in_inv_W_main(EIDs.n_elem);
    arma::field<arma::mat> in_Upsilon_cov_main(EIDs.n_elem);
    // arma::field<arma::mat> in_Upsilon_omega(EIDs.n_elem);

    omp_set_num_threads(16) ;
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;

        arma::vec b_i = B(ii);
        arma::vec vec_alpha_i = A(ii);
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1, 4);
        Y_i = Y_i.t();
        
        arma::field<arma::mat> Dn_alpha_full = Dn(ii);
        arma::field<arma::mat> Xn_full = Xn(ii);
        
        arma::mat interm_W(4,4,arma::fill::zeros);
        arma::vec interm_V(4,arma::fill::zeros);
        for(int k=1; k < n_i; k++) {
            arma::mat X_k = Xn_full(k);
            arma::mat X_k_1 = Xn_full(k-1);
            arma::vec vec_A_k = A_all_state.col(0);
            arma::mat A_1_k = arma::diagmat(vec_A_k);
            arma::mat diff_hold_x = X_k - A_1_k * X_k_1;
            
            arma::mat invR;
            if(b_i(k) == 1) {
                invR = invR_base;
            } else {
                invR = invR_bleed;
            }
            
            interm_W = interm_W + diff_hold_x.t() * invR * diff_hold_x;
            
            arma::mat D_k = Dn_alpha_full(k);
            arma::mat D_k_1 = Dn_alpha_full(k-1);
            arma::mat diff_hold_d = D_k - A_1_k * D_k_1;
            
            arma::vec m_k = Y_i.col(k) - A_1_k * Y_i.col(k-1) - diff_hold_d * vec_alpha_i;
            interm_V = interm_V + diff_hold_x.t() * invR * m_k;
        }
        
        arma::mat R;
        if(b_i(0) == 1) {
            R = R_base;
        } else {
            R = R_bleed;
        }
        arma::vec vec_A = A_all_state.col(0);
        arma::mat Gamma     = {{R(0,0) / (1 - vec_A(0) * vec_A(0)), 
                                R(0,1) / (1 - vec_A(0) * vec_A(1)), 
                                R(0,2) / (1 - vec_A(0) * vec_A(2)), 
                                R(0,3) / (1 - vec_A(0) * vec_A(3))},
                                {R(1,0) / (1 - vec_A(1) * vec_A(0)), 
                                 R(1,1) / (1 - vec_A(1) * vec_A(1)), 
                                 R(1,2) / (1 - vec_A(1) * vec_A(2)), 
                                 R(1,3) / (1 - vec_A(0) * vec_A(3))},
                                 {R(2,0) / (1 - vec_A(2) * vec_A(0)), 
                                  R(2,1) / (1 - vec_A(2) * vec_A(1)), 
                                  R(2,2) / (1 - vec_A(2) * vec_A(2)), 
                                  R(2,3) / (1 - vec_A(0) * vec_A(3))},
                                  {R(3,0) / (1 - vec_A(3) * vec_A(0)), 
                                   R(3,1) / (1 - vec_A(3) * vec_A(1)), 
                                   R(3,2) / (1 - vec_A(3) * vec_A(2)), 
                                   R(3,3) / (1 - vec_A(0) * vec_A(3))}};
        
        arma::mat inv_Gamma = arma::inv_sympd(Gamma);
        
        arma::mat W_i_inv = Xn_full(0).t() * inv_Gamma * Xn_full(0) + interm_W;
        
        arma::mat V_i = Xn_full(0).t() * inv_Gamma * (Y_i.col(0) - Dn_alpha_full(0)*vec_alpha_i) + interm_V;
        
        arma::mat hold2 = vec_alpha_i - vec_alpha_tilde;
        arma::mat in_Upsilon_cov = hold2 * hold2.t();
        
        in_V_main(ii) = V_i;
        in_inv_W_main(ii) = W_i_inv;
        in_Upsilon_cov_main(ii) = in_Upsilon_cov;
    }
    
    arma::mat sum_in_V = in_V_main(0);
    arma::mat sum_in_inv_W = in_inv_W_main(0);
    arma::mat sum_in_Upsilon_cov = in_Upsilon_cov_main(0);
    // arma::mat sum_in_Upsilon_omega = in_Upsilon_omega(0);
    for(int ii = 1; ii < EIDs.n_elem; ii++) {
        sum_in_V = sum_in_V + in_V_main(ii);
        sum_in_inv_W = sum_in_inv_W + in_inv_W_main(ii);
        sum_in_Upsilon_cov = sum_in_Upsilon_cov + in_Upsilon_cov_main(ii);
        // sum_in_Upsilon_omega = sum_in_Upsilon_omega + in_Upsilon_omega(ii);
    }
    
    arma::mat V = inv_Sigma_beta * vec_beta_0 + sum_in_V;
    
    arma::mat inv_W = inv_Sigma_beta + sum_in_inv_W;
    arma::mat W_b = arma::inv_sympd(inv_W);

    arma::mat Upsilon_cov = Psi_Upsilon + sum_in_Upsilon_cov;

    int n_sub = EIDs.n_elem;
    
    // arma::mat upsilon_omega_cov = Psi_omega + sum_in_Upsilon_omega;
    
    par.elem(vec_beta_ind - 1) = arma::mvnrnd(W_b * V, W_b);
    par.elem(vec_sigma_upsilon_ind - 1) = arma::vectorise(riwish(nu_Upsilon + n_sub, Upsilon_cov));
    // par.elem(vec_upsilon_omega_ind - 1) = arma::vectorise(riwish(nu_upsilon_omega + n_sub, upsilon_omega_cov))
    
    return par;
}

// arma::mat update_Y_i_cpp( const arma::vec EIDs, const arma::vec par, 
//                           const arma::field<arma::uvec> par_index, 
//                           const arma::field <arma::vec> A, arma::mat Y,
//                           arma::field <arma::field<arma::mat>> Dn, 
//                           const arma::field <arma::field<arma::mat>> Xn, const arma::mat otype,
//                           const arma::field <arma::mat> Dn_omega,
//                           const arma::field <arma::vec> W, arma::field <arma::vec> B) {
// 
//     // par_index KEY: (0) beta, (1) alpha_tilde, (2) upsilon, (3) vec_A, (4) R, (5) zeta,
//     //                (6) init, (7) R_bleed, (8) omega_tilde, (9) vec_upsilon_omega
//     // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
//     // "i" is the numeric EID number
//     // "ii" is the index of the EID
// 
//     arma::mat newY(Y.n_rows, 4); 
//     arma::vec eids = Y.col(0);
//     
//     arma::vec vec_beta = par.elem(par_index(0) - 1);
//     
//     arma::vec vec_R = par.elem(par_index(4) - 1);
//     arma::mat R = arma::reshape(vec_R, 4, 4);
//     arma::mat invR = arma::inv_sympd(R);
// 
//     arma::vec vec_A_total = par.elem(par_index(3) - 1);
//     arma::vec vec_A_scale = { exp(vec_A_total(0)) / (1+exp(vec_A_total(0))),
//                               exp(vec_A_total(1)) / (1+exp(vec_A_total(1))),
//                               exp(vec_A_total(2)) / (1+exp(vec_A_total(2))),
//                               exp(vec_A_total(3)) / (1+exp(vec_A_total(3))),
//                               exp(vec_A_total(4)) / (1+exp(vec_A_total(4))),
//                               exp(vec_A_total(5)) / (1+exp(vec_A_total(5))),
//                               exp(vec_A_total(6)) / (1+exp(vec_A_total(6))),
//                               exp(vec_A_total(7)) / (1+exp(vec_A_total(7))),
//                               exp(vec_A_total(8)) / (1+exp(vec_A_total(8))),
//                               exp(vec_A_total(9)) / (1+exp(vec_A_total(9))),
//                               exp(vec_A_total(10)) / (1+exp(vec_A_total(10))),
//                               exp(vec_A_total(11)) / (1+exp(vec_A_total(11)))};
//     arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 3);
//     
//     omp_set_num_threads(14) ;
//     # pragma omp parallel for
//     for (int ii = 0; ii < EIDs.n_elem; ii++) {			
//         int i = EIDs(ii);
//         
//         arma::uvec sub_ind = arma::find(eids == i);
// 
//         arma::vec b_i = B(ii);
// 
//         arma::field<arma::mat> Dn_ii = Dn(ii);
//         arma::field<arma::mat> Xn_ii = Xn(ii);
//         arma::vec vec_alpha_ii = A(ii);
// 
//         // Index of observed versus missing data
//         // 1 = observed, 0 = missing
//         arma::mat otype_i = otype.rows(sub_ind);
//         arma::mat Y_temp = Y.rows(sub_ind);
//         arma::mat Y_i = Y_temp.cols(1,4);
//         Y_i = Y_i.t();
//         otype_i = otype_i.t();
//         arma::mat Y_i_new = Y_i;
// 
//         for(int k = 0; k < Y_i.n_cols; k++) {
//             if(all(otype_i.col(k) == 1)) {
//                 Y_i_new.col(k) = Y_i.col(k);
//             } else {
//                 if(k == 0) {
//                     arma::vec vec_A = A_all_state.col(b_i(k) - 1);
//                     arma::mat Gamma     = {{R(0,0) / (1 - vec_A(0) * vec_A(0)), 
//                                             R(0,1) / (1 - vec_A(0) * vec_A(1)), 
//                                             R(0,2) / (1 - vec_A(0) * vec_A(2)), 
//                                             R(0,3) / (1 - vec_A(0) * vec_A(3))},
//                                             {R(1,0) / (1 - vec_A(1) * vec_A(0)), 
//                                              R(1,1) / (1 - vec_A(1) * vec_A(1)), 
//                                              R(1,2) / (1 - vec_A(1) * vec_A(2)), 
//                                              R(1,3) / (1 - vec_A(0) * vec_A(3))},
//                                              {R(2,0) / (1 - vec_A(2) * vec_A(0)), 
//                                               R(2,1) / (1 - vec_A(2) * vec_A(1)), 
//                                               R(2,2) / (1 - vec_A(2) * vec_A(2)), 
//                                               R(2,3) / (1 - vec_A(0) * vec_A(3))},
//                                               {R(3,0) / (1 - vec_A(3) * vec_A(0)), 
//                                                R(3,1) / (1 - vec_A(3) * vec_A(1)), 
//                                                R(3,2) / (1 - vec_A(3) * vec_A(2)), 
//                                                R(3,3) / (1 - vec_A(0) * vec_A(3))}};
//                     arma::mat inv_Gamma = arma::inv_sympd(Gamma);
//                     
//                     arma::vec vec_A_p1 = A_all_state.col(b_i(k+1) - 1);
//                     arma::mat A_p1 = arma::diagmat(vec_A_p1);
//                     
//                     arma::vec nu_k = Dn_ii(k) * vec_alpha_ii + Xn_ii(k) * vec_beta;
//                     arma::vec nu_k_p1 = Dn_ii(k+1) * vec_alpha_ii + Xn_ii(k+1) * vec_beta;
//                     
//                     arma::vec y_val_kp1 = Y_i_new.col(k+1);
// 
//                     arma::mat inv_W_i = inv_Gamma + A_p1.t() * invR * A_p1;
//                     arma::mat W_i = inv(inv_W_i);
//                     arma::vec V_i = inv_Gamma*nu_k + A_p1.t()*invR*(y_val_kp1 - nu_k_p1 + A_p1*nu_k);
// 
//                     arma::vec y_i_mean = W_i * V_i;
// 
//                     arma::vec new_value = arma::mvnrnd(y_i_mean, W_i, 1);
//                     arma::vec update_value = Y_i_new.col(k);
//                     arma::uvec ind_replace = arma::find(otype_i.col(k) == 0);
//                     update_value.elem(ind_replace) = new_value.elem(ind_replace);
// 
//                     Y_i_new.col(k) = update_value;
//                 } else if(k == Y_i.n_cols - 1) {
//                     
//                     arma::vec vec_A = A_all_state.col(b_i(k) - 1);
//                     arma::mat A_k = arma::diagmat(vec_A);
//                     
//                     arma::vec nu_k = Dn_ii(k) * vec_alpha_ii + Xn_ii(k) * vec_beta;
//                     arma::vec nu_k_m1 = Dn_ii(k-1) * vec_alpha_ii + Xn_ii(k-1) * vec_beta;
//                     
//                     arma::vec y_val_km1 = Y_i_new.col(k-1);
// 
//                     arma::vec y_i_mean = nu_k + A_k * (y_val_km1 - nu_k_m1);
// 
//                     arma::vec new_value = arma::mvnrnd(y_i_mean, R, 1);
//                     arma::vec update_value = Y_i_new.col(k);
//                     arma::uvec ind_replace = arma::find(otype_i.col(k) == 0);
//                     update_value.elem(ind_replace) = new_value.elem(ind_replace);
// 
//                     Y_i_new.col(k) = update_value;
//                 } else {
// 
//                     arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
//                     arma::vec vec_A_p1 = A_all_state.col(b_i(k+1) - 1);
//                     arma::mat A_k = arma::diagmat(vec_A_k);                    
//                     arma::mat A_p1 = arma::diagmat(vec_A_p1);
// 
//                     arma::vec nu_k    = Dn_ii(k) * vec_alpha_ii + Xn_ii(k) * vec_beta;
//                     arma::vec nu_k_m1 = Dn_ii(k-1) * vec_alpha_ii + Xn_ii(k-1) * vec_beta;
//                     arma::vec nu_k_p1 = Dn_ii(k+1) * vec_alpha_ii + Xn_ii(k+1) * vec_beta;
//                     
//                     arma::vec y_val_km1 = Y_i_new.col(k-1);
//                     arma::vec y_val_kp1 = Y_i_new.col(k+1);
// 
//                     arma::mat inv_W_i = invR + A_p1.t() * invR * A_p1;
//                     arma::mat W_i = inv(inv_W_i);
// 
//                     arma::vec V_i = invR * (nu_k + A_k * (y_val_km1 - nu_k_m1)) + 
//                                         A_p1.t() * invR * (y_val_kp1 - nu_k_p1 + A_p1 * nu_k);
// 
//                     arma::vec y_i_mean = W_i * V_i;
// 
//                     arma::vec new_value = arma::mvnrnd(y_i_mean, W_i, 1);
//                     arma::vec update_value = Y_i_new.col(k);
//                     arma::uvec ind_replace = arma::find(otype_i.col(k) == 0);
//                     update_value.elem(ind_replace) = new_value.elem(ind_replace);
// 
//                     Y_i_new.col(k) = update_value;
//                 }
//             }
//         }
// 
//         Y_i_new = Y_i_new.t();
//         newY.rows(sub_ind) = Y_i_new;
//     }
//     
//     Y.cols(1,4) = newY;
//     return Y;
// }

// [[Rcpp::export]]
Rcpp::List proposal_R_cpp(const int nu_R, const arma::mat psi_R, 
                          const arma::mat &Y, arma::field<arma::field<arma::mat>> &Dn, 
                          const arma::field<arma::field<arma::mat>> &Xn, const arma::field <arma::vec> A, 
                          const arma::vec par, const arma::field<arma::uvec> par_index, 
                          const arma::vec EIDs, arma::field <arma::vec> B){

    // par_index KEY: (0) beta, (1) alpha_tilde, (2) upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) R_bleed, (8) omega_tilde, (9) vec_upsilon_omega

    arma::vec eids = Y.col(0);
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = { exp(vec_A_total(0)) / (1+exp(vec_A_total(0))),
                              exp(vec_A_total(1)) / (1+exp(vec_A_total(1))),
                              exp(vec_A_total(2)) / (1+exp(vec_A_total(2))),
                              exp(vec_A_total(3)) / (1+exp(vec_A_total(3)))};
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 1);

    arma::vec vec_R_base = par.elem(par_index(4) - 1);
    arma::vec vec_R_bleed = par.elem(par_index(7) - 1);
    arma::mat R_base = arma::reshape(vec_R_base, 4, 4);
    arma::mat R_bleed = arma::reshape(vec_R_bleed, 4, 4);

    arma::mat psi_prop_R_interm(4, 4, arma::fill::zeros);
    int df = 0;
    arma::mat psi_prop_R_interm2(4, 4, arma::fill::zeros);
    int df2 = 0;

    // omp_set_num_threads(8) ;
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {

        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1,4);
        Y_i = Y_i.t();
        
        arma::vec b_i = B(ii);
        arma::vec vec_alpha_i = A(ii);
        arma::vec vec_beta = par.elem(par_index(0) - 1);
        arma::field<arma::mat> Xn_i = Xn(ii);
        arma::field<arma::mat> Dn_i = Dn(ii);
        
        // Removing the first index
        arma::mat Y_i_minus_1 = Y_i.cols(1,Y_i.n_cols-1);
        arma::vec b_i_minus_1 = b_i.subvec(1, b_i.n_elem-1);
        
        arma::uvec b_base = arma::find(b_i_minus_1 == 1);
        df += b_base.n_elem;
        arma::uvec b_bleed = arma::find(b_i_minus_1 != 1);
        df2 += b_bleed.n_elem;
        
        // Rcpp::Rcout << ii << " " << b_base.n_elem << " " << b_bleed.n_elem << std::endl;
        if(b_base.n_elem > 0) {
            arma::mat Y_i_base = Y_i_minus_1.cols(b_base);
            arma::mat M(4, Y_i_base.n_cols, arma::fill::zeros);
            for(int k = 0; k < Y_i_base.n_cols; k++) {
                int kk = b_base(k)+1;
                arma::vec nu_k = Dn_i(kk) * vec_alpha_i + Xn_i(kk) * vec_beta;
                arma::vec nu_k_1 = Dn_i(kk-1) * vec_alpha_i + Xn_i(kk-1) * vec_beta;
                arma::vec vec_A_k = A_all_state.col(0);
                arma::mat A_1_k = arma::diagmat(vec_A_k);
                
                arma::vec m_temp = nu_k + A_1_k * (Y_i.col(kk-1) - nu_k_1);
                M.col(k) = m_temp;
            }
            
            arma::mat hold = Y_i_base - M;
            arma::mat hold2 = hold * hold.t();
            
            psi_prop_R_interm += hold2;
        }
        
        if(b_bleed.n_elem > 0) {
            arma::mat Y_i_bleed = Y_i_minus_1.cols(b_bleed);
            arma::mat M(4, Y_i_bleed.n_cols, arma::fill::zeros);
            for(int k = 0; k < Y_i_bleed.n_cols; k++) {
                int kk = b_bleed(k)+1;
                arma::vec nu_k = Dn_i(kk) * vec_alpha_i + Xn_i(kk) * vec_beta;
                arma::vec nu_k_1 = Dn_i(kk-1) * vec_alpha_i + Xn_i(kk-1) * vec_beta;
                arma::vec vec_A_k = A_all_state.col(0);
                arma::mat A_1_k = arma::diagmat(vec_A_k);
                
                arma::vec m_temp = nu_k + A_1_k * (Y_i.col(kk-1) - nu_k_1);
                M.col(k) = m_temp;
            }
            
            arma::mat hold = Y_i_bleed - M;
            arma::mat hold2 = hold * hold.t();
            
            psi_prop_R_interm2 += hold2;
        }
    }

    arma::mat psi_prop_R = psi_prop_R_interm + psi_R;
    arma::mat psi_prop_R2 = psi_prop_R_interm2 + psi_R;
    
    int nu_prop_R = nu_R + df;
    int nu_prop_R2 = nu_R + df2;
    
    List nu_psi_R = List::create(psi_prop_R, nu_prop_R, psi_prop_R2, nu_prop_R2);
    return nu_psi_R;
}

// [[Rcpp::export]]
void test_fnc() {
    
    arma::vec test = {1,2,3,4,5,2};
    
    arma::uvec sub_ind = arma::find(test == 1);
    Rcpp::Rcout << sub_ind << std::endl;
    Rcpp::Rcout << sub_ind.n_elem << std::endl;
    arma::uvec sub_ind2 = arma::find(test == 2);
    Rcpp::Rcout << sub_ind2 << std::endl;
    Rcpp::Rcout << sub_ind2.n_elem << std::endl;
    arma::uvec sub_ind3 = arma::find(test == -1);
    Rcpp::Rcout << sub_ind3 << std::endl;
    Rcpp::Rcout << sub_ind3.n_elem << std::endl;
    
    
    arma::mat mat1 = {{.2,   0,   0,   0},
                        {0,   2,   1,   0},
                        {0,   1,   4,   0},
                        {0,   0,   0,  .2}};
    arma::mat mat2 = {{.2, -.1,  .1, -.1},     
                      { -.1,   2, -1,  .1},     
                      {.1, -1,   2, -.1},     
                      {-.1,  .1, -.1,  .2}};
    Rcpp::Rcout << mat1 << std::endl;
    Rcpp::Rcout << mat1.is_sympd() << std::endl;
    Rcpp::Rcout << mat1.is_symmetric() << std::endl;
    
    Rcpp::Rcout << mat2 << std::endl;
    Rcpp::Rcout << mat2.is_sympd() << std::endl;
    Rcpp::Rcout << mat2.is_symmetric() << std::endl;
    
    
}
