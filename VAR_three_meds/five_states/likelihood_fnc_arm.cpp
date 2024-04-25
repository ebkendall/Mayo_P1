#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// #include <omp.h>
// // [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// Defining the Omega_List as a global variable when pre-compiling ----------
const arma::mat adj_mat = { {1, 1, 0, 1, 0},
                            {0, 1, 1, 1, 0},
                            {1, 1, 1, 1, 0},
                            {0, 1, 0, 1, 1},
                            {1, 1, 0, 1, 1} };

const arma::mat adj_mat_sub = { {1, 0, 0, 1, 0},
                                {0, 1, 0, 0, 0},
                                {1, 0, 1, 1, 0},
                                {0, 0, 0, 1, 1},
                                {1, 0, 0, 1, 1} };

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

arma::vec log_f_i_cpp(const int i, const int ii, arma::vec t_pts, const arma::vec &par, 
                   const arma::field<arma::uvec> &par_index, const arma::vec &A, const arma::vec &B, 
                   const arma::mat &Y, const arma::mat &z, const arma::field<arma::mat> &Dn, 
                   const arma::field<arma::mat> &Xn, const arma::field<arma::mat> &Dn_omega, const arma::vec &W) {
  
  // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta, 
  //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
  
  // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
  // "i" is the numeric EID number
  // "ii" is the index of the EID
  double in_value = 0;
  
  double in_value_init = 0;
  
  arma::vec eids = Y.col(0);
  
  arma::vec vec_R = par.elem(par_index(4) - 1);
  arma::mat R = arma::reshape(vec_R, 4, 4);
  arma::mat invR = arma::inv_sympd(R);
  
  arma::vec vec_A_total = par.elem(par_index(3) - 1);
  arma::vec vec_A_scale = { exp(vec_A_total(0)) / (1+exp(vec_A_total(0))),
                            exp(vec_A_total(1)) / (1+exp(vec_A_total(1))),
                            exp(vec_A_total(2)) / (1+exp(vec_A_total(2))),
                            exp(vec_A_total(3)) / (1+exp(vec_A_total(3))),
                            exp(vec_A_total(4)) / (1+exp(vec_A_total(4))),
                            exp(vec_A_total(5)) / (1+exp(vec_A_total(5))),
                            exp(vec_A_total(6)) / (1+exp(vec_A_total(6))),
                            exp(vec_A_total(7)) / (1+exp(vec_A_total(7))),
                            exp(vec_A_total(8)) / (1+exp(vec_A_total(8))),
                            exp(vec_A_total(9)) / (1+exp(vec_A_total(9))),
                           exp(vec_A_total(10)) / (1+exp(vec_A_total(10))),
                           exp(vec_A_total(11)) / (1+exp(vec_A_total(11))),
                           exp(vec_A_total(12)) / (1+exp(vec_A_total(12))),
                           exp(vec_A_total(13)) / (1+exp(vec_A_total(13))),
                           exp(vec_A_total(14)) / (1+exp(vec_A_total(14))),
                           exp(vec_A_total(15)) / (1+exp(vec_A_total(15))),
                           exp(vec_A_total(16)) / (1+exp(vec_A_total(16))),
                           exp(vec_A_total(17)) / (1+exp(vec_A_total(17))),
                           exp(vec_A_total(18)) / (1+exp(vec_A_total(18))),
                           exp(vec_A_total(19)) / (1+exp(vec_A_total(19)))};
  arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); // THREE STATE
  
  arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
  arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12); // THREE STATE
  
  // The time-homogeneous probability transition matrix
  arma::uvec sub_ind = arma::find(eids == i);
  arma::mat z_i = z.rows(sub_ind.min(), sub_ind.max());
  int n_i = z_i.n_rows;
  
  // The state transition likelihood component for current iterate of b_i
  arma::vec b_i = B;
  arma::vec vec_init_content = par.elem(par_index(6) - 1);
  arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                             exp(vec_init_content(2)), exp(vec_init_content(3))}; // THREE STATE
  arma::vec P_init = init_logit / arma::accu(init_logit); 
  
  // Subsetting the data
  arma::mat Y_temp = Y.rows(sub_ind);
  arma::mat Y_i = Y_temp.cols(1, 4);
  Y_i = Y_i.t();
  
  arma::field<arma::mat> Dn_alpha_full = Dn;
  arma::field<arma::mat> Dn_omega_full = Dn_omega;
  arma::field<arma::mat> Xn_full = Xn;
  arma::vec vec_alpha_ii = A;
  arma::vec vec_omega_ii = W;
  
  arma::mat vec_beta = par.elem(par_index(0) - 1);
  
  // Full likelihood evaluation is not needed for updating pairs of b_i components
  if (any(t_pts == -1)) { t_pts = arma::linspace(1, n_i, n_i);}
  
  for(int w=0; w < t_pts.n_elem;++w){
    int k = t_pts(w);
    if(k==1){
        // State space component
        int b_0 = b_i(0);
        
        arma::vec vec_A = A_all_state.col(b_0 - 1);
        
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
        arma::vec nu_1 = Dn_alpha_full(0) * vec_alpha_ii + 
                            Dn_omega_full(0) * vec_omega_ii +
                              Xn_full(0) * vec_beta;
        arma::vec log_y_pdf = dmvnorm(y_1.t(), nu_1, Gamma, true);
        
        in_value_init = in_value_init + log(P_init(b_0 - 1)) + arma::as_scalar(log_y_pdf);
        
        in_value = in_value + log(P_init(b_0 - 1)) + arma::as_scalar(log_y_pdf);
    } else{
        // State space component
        double q1_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(0));
        double q1 = exp(q1_sub);
        double q2_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(1));
        double q2 = exp(q2_sub);
        double q3_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(2));
        double q3 = exp(q3_sub);
        double q4_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(3));
        double q4 = exp(q4_sub);
        
        double q5_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(4));
        double q5 = exp(q5_sub);
        double q6_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(5));
        double q6 = exp(q6_sub);
        double q7_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(6));
        double q7 = exp(q7_sub);
        double q8_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(7));
        double q8 = exp(q8_sub);
        
        double q9_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(8));
        double q9 = exp(q9_sub);
        double q10_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(9));
        double q10 = exp(q10_sub);
        double q11_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(10));
        double q11 = exp(q11_sub);
        double q12_sub = arma::as_scalar(z_i.row(k-1) * zeta.col(11));
        double q12 = exp(q12_sub);
        
        arma::mat Q = {{   1,   q1,  0,  q2,  0},
                       {   0,    1, q3,  q4,  0},
                       {  q5,   q6,  1,  q7,  0},
                       {   0,   q8,  0,   1, q9},
                       { q10,  q11,  0, q12,  1}}; // THREE STATE
        
        arma::vec q_row_sums = arma::sum(Q, 1);
        arma::mat P_i = Q.each_col() / q_row_sums;
        int b_k_1 = b_i(k-2);
        int b_k = b_i(k-1);
        
        arma::vec vec_A = A_all_state.col(b_k - 1);
        arma::mat A_1 = arma::diagmat(vec_A);
        
        arma::vec y_k_1 = Y_i.col(k-2);
        arma::vec y_k = Y_i.col(k-1);
        arma::vec nu_k_1 = Dn_alpha_full(k-2) * vec_alpha_ii + 
                                Dn_omega_full(k-2) * vec_omega_ii +
                                    Xn_full(k-2) * vec_beta;
        arma::vec nu_k = Dn_alpha_full(k-1) * vec_alpha_ii + 
                            Dn_omega_full(k-1) * vec_omega_ii +
                                Xn_full(k-1) * vec_beta;
        
        arma::vec mean_k = nu_k + A_1 * (y_k_1 - nu_k_1);
        
        arma::vec log_y_k_pdf = dmvnorm(y_k.t(), mean_k, R, true);
        
        in_value = in_value + log(P_i( b_k_1 - 1, b_k - 1)) + arma::as_scalar(log_y_k_pdf);
    }
  }
  
  arma::vec in_value_vec = {in_value, in_value_init}; 
  
  return in_value_vec;
}

double log_f_i_cpp_total(const arma::vec &EIDs, arma::vec t_pts, const arma::vec &par, const arma::field<arma::uvec> &par_index, 
                         const arma::field <arma::vec> &A, const arma::field <arma::vec> &B, 
                         const arma::mat &Y, const arma::mat &z, const arma::field<arma::field<arma::mat>> &Dn, 
                         const arma::field <arma::field<arma::mat>> &Xn, const arma::field<arma::field<arma::mat>> &Dn_omega, 
                         const arma::field <arma::vec> &W, int n_cores) {

  // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
  //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
  // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
  // "i" is the numeric EID number
  // "ii" is the index of the EID
  arma::vec in_vals(EIDs.n_elem, arma::fill::zeros);
  arma::vec in_vals_init(EIDs.n_elem, arma::fill::zeros);

    // omp_set_num_threads(n_cores);
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        arma::vec in_val_vec = log_f_i_cpp(i, ii, t_pts, par, par_index, A(ii), B(ii), Y, z, Dn(ii), Xn(ii), Dn_omega(ii), W(ii));
        in_vals(ii) = in_val_vec(0);
        in_vals_init(ii) = in_val_vec(1);
    }
    
    double in_value = arma::accu(in_vals);
    double in_value_init = arma::accu(in_vals_init);
    
    // Rcpp::Rcout << "likelihood from initial time point = " << in_value_init << std::endl;
    // Rcpp::Rcout << "likelihood from all time points    = " << in_value << std::endl;
    
    return in_value;
}

// [[Rcpp::export]]
double log_post_cpp(const arma::vec &EIDs, const arma::vec &par, const arma::field<arma::uvec> &par_index,
                    const arma::field<arma::vec> &A, const arma::field<arma::vec> &B,
                    const arma::mat &Y, const arma::mat &z, const arma::field<arma::field<arma::mat>> &Dn,
                    const arma::field<arma::field<arma::mat>> &Xn, const arma::field<arma::field<arma::mat>> &Dn_omega, 
                    const arma::field<arma::vec> &W, int n_cores) {

  // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta, 
  //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
  // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
  // "i" is the numeric EID number
  // "ii" is the index of the EID

  // Compute the likelihood ----------------------------------------------------
  double value;
  arma::vec t_pts = {-1};
  value = log_f_i_cpp_total(EIDs, t_pts, par, par_index, A, B, Y, z, Dn, Xn, Dn_omega, W, n_cores);
  // ---------------------------------------------------------------------------

  // Compute prior densities of all necessary model parameters -----------------
  // Zeta priors ---------------------------------------------------------------
  arma::uvec vec_zeta_ind = par_index(5);
  arma::vec vec_zeta_content = par.elem(vec_zeta_ind - 1);
  arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12); // THREE STATE

  
  // transitions:                    1->2,         1->4,         2->3,         2->4, 
  //                                 3->1,         3->2,         3->4,         4->2, 
  //                                 4->5,         5->1,         5->2,         5->4
  arma::vec vec_zeta_mean = {-7.2405, 1.5, -5.2152,   1, -2.6473,  -1, -5.1475,  -1, 
                             -9.4459,  -1, -7.2404,   2, -5.2151,   1, -7.1778, 1.5, 
                             -2.6523,   0, -9.4459,  -1, -7.2404, 1.5, -5.2151,   1}; 
  arma::vec scalar_1 = {0.25, 1, 0.25, 1, 0.50, 1, 0.50, 1, 
                        0.70, 1, 0.50, 1, 0.70, 1, 0.25, 1, 
                        0.50, 1, 0.70, 1, 0.50, 1, 0.70, 1};
  arma::mat zeta_sd = arma::diagmat(scalar_1);

  arma::vec prior_zeta = dmvnorm(vec_zeta_content.t(), vec_zeta_mean, zeta_sd, true);
  double prior_zeta_val = arma::as_scalar(prior_zeta);

  // Initial Probabilities priors ----------------------------------------------
  arma::uvec vec_init_ind = par_index(6);
  arma::vec vec_init_content = par.elem(vec_init_ind - 1);
  arma::vec vec_init_mean = {0, 0, 0, 0}; // THREE STATE
  arma::vec scalar_2(vec_init_content.n_elem, arma::fill::ones); // THREE STATE
  scalar_2 = 10 * scalar_2;
  arma::mat init_sd = arma::diagmat(scalar_2);

  arma::vec prior_init = dmvnorm(vec_init_content.t(), vec_init_mean, init_sd, true);
  double prior_init_val = arma::as_scalar(prior_init);
  
  // A_1 priors ----------------------------------------------------------------
  arma::vec vec_A1_content = par.elem(par_index(3) - 1);
  
  arma::vec vec_A1_mean(vec_A1_content.n_elem, arma::fill::zeros);
  arma::vec A1_scalar(vec_A1_content.n_elem, arma::fill::ones);
  A1_scalar = 2 * A1_scalar;
  arma::mat A1_sd = arma::diagmat(A1_scalar);
  
  arma::vec prior_A1 = dmvnorm(vec_A1_content.t(), vec_A1_mean, A1_sd, true);
  double prior_A1_val = arma::as_scalar(prior_A1);
  
  // R priors ------------------------------------------------------------------
  arma::vec vec_R_content = par.elem(par_index(4) - 1);
  arma::mat R = arma::reshape(vec_R_content, 4, 4);
  
  int nu_R = 10;
  //   arma::mat Psi_R(4,4,arma::fill::eye);
  arma::vec scalar_vec_R = {4.58, 98.2, 101.3, 7.6};
  scalar_vec_R = (nu_R - 4 - 1) * scalar_vec_R;
  arma::mat Psi_R = arma::diagmat(scalar_vec_R);
  
  double prior_R_val = diwish(R, nu_R, Psi_R, true);
  
  // Upsilon omega priors -------------------------------------------------------------
  arma::vec vec_up_omega_content = par.elem(par_index(8) - 1);
  arma::vec omega_mean(vec_up_omega_content.n_elem, arma::fill::zeros);
  
  arma::vec diag_omega_sd(vec_up_omega_content.n_elem, arma::fill::ones);
  diag_omega_sd = 2 * diag_omega_sd;
  arma::mat omega_sd = arma::diagmat(diag_omega_sd);

  arma::vec prior_omega = dmvnorm(vec_up_omega_content.t(), omega_mean, omega_sd, true);
  double prior_omega_val = arma::as_scalar(prior_omega);
  // ---------------------------------------------------------------------------
  
  
  value = value + prior_A1_val + prior_R_val + prior_zeta_val + prior_init_val + prior_omega_val;
  
  // Rcpp::Rcout << "total log posterior = " << value << std::endl;
  return value;
}

// [[Rcpp::export]]
Rcpp::List update_b_i_cpp(const arma::vec EIDs, const arma::vec &par, 
                          const arma::field<arma::uvec> &par_index, 
                          const arma::field <arma::vec> &A, 
                          arma::field <arma::vec> &B, 
                          const arma::mat &Y, const arma::mat &z, 
                          arma::field<arma::field<arma::mat>> &Dn, 
                          const arma::field <arma::field<arma::mat>> &Xn, 
                          const arma::field<arma::field<arma::mat>> &Dn_omega, 
                          const arma::field <arma::vec> &W,
                          const arma::vec &bleed_indicator, int n_cores) {

  // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
  //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
  // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
  // "i" is the numeric EID number
  // "ii" is the index of the EID

  arma::vec eids = Y.col(0); 
  arma::vec rbc_rule_vec = Y.col(5);
  arma::vec clinic_rule_vec = Y.col(6); 
  
  arma::field<arma::vec> B_return(EIDs.n_elem);
  arma::field<arma::field<arma::mat>> Dn_return(EIDs.n_elem);
  
  // omp_set_num_threads(n_cores);
  // # pragma omp parallel for
  for (int ii = 0; ii < EIDs.n_elem; ii++) {
    int i = EIDs(ii);
    arma::uvec sub_ind = arma::find(eids == i);
    
    int n_i = sub_ind.n_elem;
    
    int rbc_rule = rbc_rule_vec(sub_ind.min());
    int clinic_rule = clinic_rule_vec(sub_ind.min());
    
    // Subsetting fields
    arma::vec B_temp = B(ii);
    arma::vec A_temp = A(ii);
    arma::vec W_temp = W(ii);
    arma::vec bleed_ind_i = bleed_indicator.elem(sub_ind);
    
    arma::field<arma::mat> Dn_temp = Dn(ii);
    arma::field<arma::mat> Dn_omega_temp = Dn_omega(ii);
    arma::field<arma::mat> Xn_temp = Xn(ii);
    
    // Subsetting the remaining data
    arma::mat Y_temp = Y.rows(sub_ind);
    arma::mat z_temp = z.rows(sub_ind);
    
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
      
      // State sampling using RBC and Clinic information -----------------------
      bool valid_prop = false;
      bool b_i_rule = arma::any(arma::vectorise(pr_B)==2);
      
      if(clinic_rule == 1) {
          if(rbc_rule == 1) {
              // clinic=1, rbc=1 -> NEED S2, *yes* restriction on time of S2
              int pos_bleed = arma::as_scalar(arma::find(bleed_ind_i == 1));
              arma::uvec b_i_time = arma::find(pr_B == 2);
              if(arma::any(b_i_time <= pos_bleed)) {
                  valid_prop = true;
              }
          } else {
              // clinic=1, rbc=0 -> NEED S2, *no* restriction on time of S2
              if(b_i_rule) {
                  valid_prop = true;
              }
          }
      } else if(clinic_rule == 0) {
          if(rbc_rule == 1) {
              // clinic=0, rbc=1 -> NEED S2, *yes* restriction on time of S2
              int pos_bleed = arma::as_scalar(arma::find(bleed_ind_i == 1));
              arma::uvec b_i_time = arma::find(pr_B == 2);
              if(arma::any(b_i_time <= pos_bleed)) {
                  valid_prop = true;
              }
          } else {
              // clinic=0, rbc=0 -> No restrictions, consider all state seq.
              valid_prop = true;
          }
      } else {
          // clinic=-1, rbc=1 -> evaluate likelihood anyways because S1,S4,S5
          // clinic=-1, rbc=0 -> evaluate likelihood anyways because S1,S4,S5
          valid_prop = true; 
      }
      
      // If the proposed state sequence is the same, then we do not need to 
      // evaluate the likelihood. Thus valid_prop = false can be set.
      if(arma::accu(pr_B == B_temp) == pr_B.n_elem) {
          valid_prop = false;
      }
      // -----------------------------------------------------------------------
      
      if(valid_prop) {
          
        arma::vec log_prev_vec = log_f_i_cpp(i, ii, t_pts, par, par_index,A_temp,
                                             B_temp,Y_temp,z_temp,Dn_temp,Xn_temp,
                                             Dn_omega_temp, W_temp);
          
        double log_target_prev = log_prev_vec(0);
    
        arma::vec twos(pr_B.n_elem, arma::fill::zeros);
        arma::vec threes = twos; // THREE STATE
        arma::vec fours = twos;
        arma::vec fives = twos;
        
        twos.elem(arma::find(pr_B == 2)) += 1;
        threes.elem(arma::find(pr_B == 3)) += 1; // THREE STATE
        fours.elem(arma::find(pr_B == 4)) += 1;
        fives.elem(arma::find(pr_B == 5)) += 1;
        
        arma::vec ones(pr_B.n_elem, arma::fill::ones);
        
        arma::mat bigB = arma::join_rows(ones, arma::cumsum(twos));
        bigB = arma::join_rows(bigB, arma::cumsum(threes)); // THREE STATE
        bigB = arma::join_rows(bigB, arma::cumsum(fours));
        bigB = arma::join_rows(bigB, arma::cumsum(fives));
        
        arma::mat I = arma::eye(4,4);
        for(int jj = 0; jj < n_i; jj++) {
            pr_Dn(jj) = arma::kron(I, bigB.row(jj));
        }
        
        arma::vec log_target_vec = log_f_i_cpp(i,ii,t_pts,par,par_index,A_temp,
                                     pr_B,Y_temp,z_temp,pr_Dn,Xn_temp,
                                     Dn_omega_temp, W_temp);
        double log_target = log_target_vec(0);
        
        // Note that the proposal probs cancel in the MH ratio
        double diff_check = log_target - log_target_prev;
        double min_log = log(arma::randu(arma::distr_param(0,1)));
        if(diff_check > min_log){
          B_temp = pr_B;
          Dn_temp = pr_Dn;
        }
      }
    }
    B_return(ii) = B_temp;
    Dn_return(ii) = Dn_temp;
  }
  List B_Dn = List::create(B_return, Dn_return);
  
  return B_Dn;
}

// [[Rcpp::export]]
Rcpp::List update_Dn_Xn_cpp( const arma::vec EIDs, arma::field <arma::vec> &B, 
                             const arma::mat &Y, const arma::vec &par, 
                             const arma::field<arma::uvec> &par_index,
                             const arma::vec &x, int n_cores) {

  // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
  //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
  // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
  // "i" is the numeric EID number
  // "ii" is the index of the EID

  arma::vec eids = Y.col(0);
  arma::field<arma::field<arma::mat>> Dn(EIDs.n_elem);
  arma::field<arma::field<arma::mat>> Xn(EIDs.n_elem);
  
  // omp_set_num_threads(n_cores);
  // # pragma omp parallel for
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
    arma::vec fours = twos;
    arma::vec fives = twos;
    
    twos.elem(arma::find(b_i == 2)) += 1;
    threes.elem(arma::find(b_i == 3)) += 1; // THREE STATE
    fours.elem(arma::find(b_i == 4)) += 1;
    fives.elem(arma::find(b_i == 5)) += 1;
    
    arma::vec ones(b_i.n_elem, arma::fill::ones);
    
    arma::mat bigB = arma::join_rows(ones, arma::cumsum(twos));
    bigB = arma::join_rows(bigB, arma::cumsum(threes)); // THREE STATE
    bigB = arma::join_rows(bigB, arma::cumsum(fours));
    bigB = arma::join_rows(bigB, arma::cumsum(fives));
    
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
arma::field <arma::vec> update_alpha_i_cpp( const arma::vec &EIDs, const arma::vec &par, 
                                            const arma::field<arma::uvec> &par_index,
                                            const arma::mat &Y, arma::field<arma::field<arma::mat>> &Dn, 
                                            const arma::field <arma::field<arma::mat>> &Xn,
                                            const arma::field<arma::field<arma::mat>> &Dn_omega, 
                                            const arma::field <arma::vec> &W, 
                                            arma::field <arma::vec> &B, int n_cores){

  // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
  //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
  // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
  // "i" is the numeric EID number
  // "ii" is the index of the EID

  arma::vec eids = Y.col(0);

  arma::vec vec_R = par.elem(par_index(4) - 1);
  arma::mat R = arma::reshape(vec_R, 4, 4);
  arma::mat invR = arma::inv_sympd(R);

  arma::vec vec_alpha_tilde = par.elem(par_index(1) - 1);

  arma::vec vec_beta = par.elem(par_index(0) - 1);

  arma::vec sigma_upsilon_vec = par.elem(par_index(2) - 1);
  arma::mat Upsilon = arma::reshape(sigma_upsilon_vec, 20, 20); // THREE STATE
  arma::mat inv_Upsilon = arma::inv_sympd(Upsilon);

  arma::field<arma::vec> A(EIDs.n_elem);
  
  arma::vec vec_A_total = par.elem(par_index(3) - 1);
  arma::vec vec_A_scale = { exp(vec_A_total(0)) / (1+exp(vec_A_total(0))),
                            exp(vec_A_total(1)) / (1+exp(vec_A_total(1))),
                            exp(vec_A_total(2)) / (1+exp(vec_A_total(2))),
                            exp(vec_A_total(3)) / (1+exp(vec_A_total(3))),
                            exp(vec_A_total(4)) / (1+exp(vec_A_total(4))),
                            exp(vec_A_total(5)) / (1+exp(vec_A_total(5))),
                            exp(vec_A_total(6)) / (1+exp(vec_A_total(6))),
                            exp(vec_A_total(7)) / (1+exp(vec_A_total(7))),
                            exp(vec_A_total(8)) / (1+exp(vec_A_total(8))),
                            exp(vec_A_total(9)) / (1+exp(vec_A_total(9))),
                            exp(vec_A_total(10)) / (1+exp(vec_A_total(10))),
                            exp(vec_A_total(11)) / (1+exp(vec_A_total(11))),
                            exp(vec_A_total(12)) / (1+exp(vec_A_total(12))),
                            exp(vec_A_total(13)) / (1+exp(vec_A_total(13))),
                            exp(vec_A_total(14)) / (1+exp(vec_A_total(14))),
                            exp(vec_A_total(15)) / (1+exp(vec_A_total(15))),
                            exp(vec_A_total(16)) / (1+exp(vec_A_total(16))),
                            exp(vec_A_total(17)) / (1+exp(vec_A_total(17))),
                            exp(vec_A_total(18)) / (1+exp(vec_A_total(18))),
                            exp(vec_A_total(19)) / (1+exp(vec_A_total(19)))};
  arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); // THREE STATE
  
  // omp_set_num_threads(n_cores);
  // # pragma omp parallel for
  for (int ii = 0; ii < EIDs.n_elem; ii++) {
      int i = EIDs(ii);

      arma::uvec sub_ind = arma::find(eids == i);
      int n_i = sub_ind.n_elem;
      
      arma::vec b_i = B(ii);

      arma::mat Y_temp = Y.rows(sub_ind);
      arma::mat Y_i = Y_temp.cols(1, 4);
      Y_i = Y_i.t();
      
      arma::field<arma::mat> Dn_alpha_full = Dn(ii);
      arma::field<arma::mat> Dn_omega_full = Dn_omega(ii);
      arma::field<arma::mat> Xn_full = Xn(ii);
      arma::vec vec_omega_ii = W(ii);
      
      arma::mat interm_W(20,20,arma::fill::zeros); // THREE STATE
      arma::vec interm_V(20,arma::fill::zeros); // THREE STATE
      for(int k=1; k < n_i; k++) {
          arma::mat D_k = Dn_alpha_full(k);
          arma::mat D_k_1 = Dn_alpha_full(k-1);
          
          arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
          arma::mat A_1_k = arma::diagmat(vec_A_k);
          
          arma::mat diff_hold_d = D_k - A_1_k * D_k_1;
          interm_W = interm_W + diff_hold_d.t() * invR * diff_hold_d;
          
          arma::mat X_k = Xn_full(k);
          arma::mat X_k_1 = Xn_full(k-1);
          arma::mat diff_hold_x = X_k - A_1_k * X_k_1;
          
          arma::mat D_k_o = Dn_omega_full(k);
          arma::mat D_k_o_1 = Dn_omega_full(k-1);
          arma::mat diff_hold_d_o = D_k_o - A_1_k * D_k_o_1;
          
          arma::vec m_k = Y_i.col(k) - A_1_k * Y_i.col(k-1) - 
                            diff_hold_x * vec_beta - diff_hold_d_o*vec_omega_ii;
          interm_V = interm_V + diff_hold_d.t() * invR * m_k;
      }
      
      
      
      arma::vec vec_A = A_all_state.col(b_i(0) - 1);
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
          Dn_alpha_full(0).t() * inv_Gamma * (Y_i.col(0) - Xn_full(0) * vec_beta 
                                                  - Dn_omega_full(0) * vec_omega_ii) + interm_V;

      arma::vec mu = W_i * V_i;

      arma::mat alpha_i = rmvnorm(1, mu, W_i);
      arma::vec vec_alpha_i = alpha_i.t();
      
      int count_while_loop = 0;
      int count_while_loop_big = 0;
      while(vec_alpha_i(1) > 0) {
          alpha_i = rmvnorm(1, mu, W_i);
          vec_alpha_i = alpha_i.t();
          
          count_while_loop += 1;
          if(count_while_loop > 10000) {
              count_while_loop_big += 1;
              Rcpp::Rcout << "stuck in alpha, i = " << ii << ", " << count_while_loop_big << std::endl;
              count_while_loop = 0;
          }
          if(count_while_loop_big > 1000) {
              break;
          }
      }

      A(ii) = vec_alpha_i;
  }
  
  return A;
  
}

// [[Rcpp::export]]
arma::field <arma::vec> update_omega_i_cpp( const arma::vec &EIDs, const arma::vec &par, 
                                            const arma::field<arma::uvec> &par_index,
                                            const arma::mat &Y, arma::field<arma::field<arma::mat>> &Dn, 
                                            const arma::field <arma::field<arma::mat>> &Xn,
                                            const arma::field<arma::field<arma::mat>> &Dn_omega, 
                                            const arma::field <arma::vec> &A, 
                                            arma::field <arma::vec> &B, int n_cores){

    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    arma::vec eids = Y.col(0);
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat invR = arma::inv_sympd(R);
    
    arma::vec vec_omega_tilde = par.elem(par_index(7) - 1);
    
    arma::vec vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_upsilon_omega = par.elem(par_index(8) - 1);
    arma::vec vec_upsilon_omega_inv = 1 / exp(vec_upsilon_omega);
    arma::mat inv_Upsilon = arma::diagmat(vec_upsilon_omega_inv);
    
    arma::field<arma::vec> W(EIDs.n_elem);
    
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = { exp(vec_A_total(0)) / (1+exp(vec_A_total(0))),
                              exp(vec_A_total(1)) / (1+exp(vec_A_total(1))),
                              exp(vec_A_total(2)) / (1+exp(vec_A_total(2))),
                              exp(vec_A_total(3)) / (1+exp(vec_A_total(3))),
                              exp(vec_A_total(4)) / (1+exp(vec_A_total(4))),
                              exp(vec_A_total(5)) / (1+exp(vec_A_total(5))),
                              exp(vec_A_total(6)) / (1+exp(vec_A_total(6))),
                              exp(vec_A_total(7)) / (1+exp(vec_A_total(7))),
                              exp(vec_A_total(8)) / (1+exp(vec_A_total(8))),
                              exp(vec_A_total(9)) / (1+exp(vec_A_total(9))),
                              exp(vec_A_total(10)) / (1+exp(vec_A_total(10))),
                              exp(vec_A_total(11)) / (1+exp(vec_A_total(11))),
                              exp(vec_A_total(12)) / (1+exp(vec_A_total(12))),
                              exp(vec_A_total(13)) / (1+exp(vec_A_total(13))),
                              exp(vec_A_total(14)) / (1+exp(vec_A_total(14))),
                              exp(vec_A_total(15)) / (1+exp(vec_A_total(15))),
                              exp(vec_A_total(16)) / (1+exp(vec_A_total(16))),
                              exp(vec_A_total(17)) / (1+exp(vec_A_total(17))),
                              exp(vec_A_total(18)) / (1+exp(vec_A_total(18))),
                              exp(vec_A_total(19)) / (1+exp(vec_A_total(19)))};
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); // THREE STATE
    
    // omp_set_num_threads(n_cores);
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        arma::vec b_i = B(ii);
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1, 4);
        Y_i = Y_i.t();
        
        arma::field<arma::mat> Dn_alpha_full = Dn(ii);
        arma::field<arma::mat> Dn_omega_full = Dn_omega(ii);
        arma::field<arma::mat> Xn_full = Xn(ii);
        arma::vec vec_alpha_ii = A(ii);

        arma::mat interm_W(vec_omega_tilde.n_elem, vec_omega_tilde.n_elem, arma::fill::zeros);
        arma::vec interm_V(vec_omega_tilde.n_elem, arma::fill::zeros);
        for(int k=1; k < n_i; k++) {
            arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
            arma::mat A_1_k = arma::diagmat(vec_A_k);
            
            arma::mat D_k_o = Dn_omega_full(k);
            arma::mat D_k_o_1 = Dn_omega_full(k-1);
            arma::mat diff_hold_d_o = D_k_o - A_1_k * D_k_o_1;
            
            interm_W = interm_W + diff_hold_d_o.t() * invR * diff_hold_d_o;
            
            arma::mat X_k = Xn_full(k);
            arma::mat X_k_1 = Xn_full(k-1);
            arma::mat diff_hold_x = X_k - A_1_k * X_k_1;
            
            arma::mat D_k = Dn_alpha_full(k);
            arma::mat D_k_1 = Dn_alpha_full(k-1);
            arma::mat diff_hold_d = D_k - A_1_k * D_k_1;
            
            arma::vec m_k = Y_i.col(k) - A_1_k * Y_i.col(k-1) - 
                                diff_hold_x * vec_beta - diff_hold_d * vec_alpha_ii;
            interm_V = interm_V + diff_hold_d_o.t() * invR * m_k;
        }
        
        arma::vec vec_A = A_all_state.col(b_i(0) - 1);
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
        
        arma::mat W_i_inv = inv_Upsilon + Dn_omega_full(0).t() * inv_Gamma * Dn_omega_full(0) + interm_W;
        
        arma::mat W_i = inv(W_i_inv);
        
        arma::mat V_i = inv_Upsilon * vec_omega_tilde + 
                            Dn_omega_full(0).t() * inv_Gamma * 
                                (Y_i.col(0) - Xn_full(0) * vec_beta - Dn_alpha_full(0) * vec_alpha_ii) + interm_V;
        
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

    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta, 
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // The prior mean for vec_alpha_tilde
    arma::vec vec_alpha_tilde_0 = {9.57729783,          -1,        0.1, 0, 0,
                                  88.69780576,  5.04150472,         -4, 0, 0,
                                  79.74903940, -5.04150472,          4, 0, 0,
                                    5.2113319,   0.5360813, -0.6866748, 0, 0}; // THREE STATE
    
    // The prior PRECISION matrix for vec_alpha_tilde
    // arma::vec inv_Sigma_alpha_diag = {1, 1, 1, 0.0025, 0.01, 0.01, 0.0025, 0.01, 0.01, 1, 1, 1};
    arma::vec inv_Sigma_alpha_diag = { 0.1, 0.3, 0.5, 0.05, 0.05,
                                      0.05, 0.5, 0.5, 0.05, 0.05,
                                      0.05, 0.5, 0.5, 0.05, 0.05,
                                       0.1, 0.3, 0.5, 0.05, 0.05}; // THREE STATE
    
    arma::mat inv_Sigma_alpha = arma::diagmat(inv_Sigma_alpha_diag);
    
    arma::vec sigma_upsilon_vec = par.elem(par_index(2) - 1);
    arma::mat sigma_upsilon = arma::reshape(sigma_upsilon_vec, 20, 20); // THREE STATE
    
    arma::mat Upsilon = sigma_upsilon;
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
    arma::vec alpha_tilde_temp = arma::mvnrnd(mu, U, 1);
    
    // // Ensuring that the parent means for HR and MAP obey the directional assumption
    // int count_while_loop = 0;
    // int count_while_loop_big = 0;
    // while(alpha_tilde_temp(1) > alpha_tilde_temp(2)) {
    //     alpha_tilde_temp = arma::mvnrnd(mu, U, 1);
    //     count_while_loop += 1;
    //     if(count_while_loop > 1000) {
    //         count_while_loop_big += 1;
    //         Rcpp::Rcout << "stuck in alpha tilde " << count_while_loop_big << std::endl;
    //         // Rcpp::Rcout << "hr bleed: " << alpha_tilde_temp(4) << ", hr recov: " << alpha_tilde_temp(5) <<
    //         //     ", MAP bleed: " << alpha_tilde_temp(7) << ", MAP recov: " << alpha_tilde_temp(8) << std::endl;
    //         count_while_loop = 0;
    //     }
    // }
    
    par.elem(vec_alpha_tilde_ind - 1) = alpha_tilde_temp;
    
    return par;
}

// [[Rcpp::export]]
arma::vec update_omega_tilde_cpp( const arma::vec EIDs, arma::vec par, 
                                  const arma::field<arma::uvec> par_index,
                                  const arma::field <arma::vec> W, const arma::mat Y){

    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // The prior mean for vec_omega_tilde
    arma::vec vec_omega_tilde_0 = {-1, -1,  1, -1,  1,  1, -1, -1, -1,  1,  1,  1,  1,
                                   -1,  1,  1,  1, -1,  1, -1,  1, -1, -1, -1, -1, -1,
                                   -1, -1, -1, -1,  1, -1,  1, -1, -1,  1,  1, -1,  1,
                                   -1, -1,  1,  1, -1, -1, -1, -1, -1, -1,  1, -1,  1,
                                    1, -1,  1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,
                                   -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1, -1,
                                   -1, -1,  1, -1, -1, -1, -1, -1, -1,  1};
    vec_omega_tilde_0 = 3 * vec_omega_tilde_0;
    
    // The prior PRECISION matrix for vec_omega_tilde
    arma::vec inv_sig_omega_vec(vec_omega_tilde_0.n_elem, arma::fill::ones);
    inv_sig_omega_vec = 0.1 * inv_sig_omega_vec;
    arma::mat inv_Sigma_omega = arma::diagmat(inv_sig_omega_vec);
    
    arma::vec vec_upsilon_omega = par.elem(par_index(8) - 1);
    arma::vec vec_upsilon_omega_inv = 1 / exp(vec_upsilon_omega);
    arma::mat inv_Upsilon = arma::diagmat(vec_upsilon_omega_inv);
    
    int length_EIDs = EIDs.n_elem;
    
    arma::mat inv_U = inv_Upsilon * length_EIDs + inv_Sigma_omega;
    arma::mat U = inv(inv_U);
    
    arma::mat total_omega = W(0);
    for (int i = 1; i < W.n_elem; i++) {
        total_omega = arma::join_horiz(total_omega, W(i));
    }
    
    arma::mat sum_omega = arma::sum(total_omega, 1);
    
    arma::mat hold = inv_Sigma_omega * vec_omega_tilde_0 + inv_Upsilon * sum_omega;
    
    arma::vec mu = U * hold;
    
    arma::uvec vec_omega_tilde_ind = par_index(7);
    par.elem(vec_omega_tilde_ind - 1) = rmvnorm(1, mu, U);
    
    return par;
}

// [[Rcpp::export]]
arma::vec update_beta_Upsilon_R_cpp( const arma::vec &EIDs, arma::vec par, 
                                     const arma::field<arma::uvec> &par_index,
                                     const arma::field <arma::vec> &A, const arma::mat &Y,
                                     const arma::field<arma::field<arma::mat>> &Dn, 
                                     const arma::field <arma::field<arma::mat>> &Xn, 
                                     const arma::field<arma::field<arma::mat>> &Dn_omega, 
                                     const arma::field <arma::vec> &W, arma::field <arma::vec> &B,
                                     int n_cores) {
    // Conjugate updates for beta and sigma_upsilon
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta, 
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID

    // The prior info for vec_beta --------------------------------------------
    arma::uvec vec_beta_ind = par_index(0);
    arma::vec vec_beta_0 = {0.25, -1, 2, -0.25};
    // We know one unit of RBC leads to 1 unit increase in hemo in 1 hour

    arma::vec scalar_mult(vec_beta_ind.n_elem, arma::fill::ones);
    scalar_mult.fill(0.01);
    // Make the prior tighter for RBC effect on hemoglobin
    scalar_mult(0) = 4;
    arma::mat inv_Sigma_beta = arma::diagmat(scalar_mult);

    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat invR = arma::inv_sympd(R);

    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = { exp(vec_A_total(0)) / (1+exp(vec_A_total(0))),
                                exp(vec_A_total(1)) / (1+exp(vec_A_total(1))),
                                exp(vec_A_total(2)) / (1+exp(vec_A_total(2))),
                                exp(vec_A_total(3)) / (1+exp(vec_A_total(3))),
                                exp(vec_A_total(4)) / (1+exp(vec_A_total(4))),
                                exp(vec_A_total(5)) / (1+exp(vec_A_total(5))),
                                exp(vec_A_total(6)) / (1+exp(vec_A_total(6))),
                                exp(vec_A_total(7)) / (1+exp(vec_A_total(7))),
                                exp(vec_A_total(8)) / (1+exp(vec_A_total(8))),
                                exp(vec_A_total(9)) / (1+exp(vec_A_total(9))),
                                exp(vec_A_total(10)) / (1+exp(vec_A_total(10))),
                                exp(vec_A_total(11)) / (1+exp(vec_A_total(11))),
                                exp(vec_A_total(12)) / (1+exp(vec_A_total(12))),
                                exp(vec_A_total(13)) / (1+exp(vec_A_total(13))),
                                exp(vec_A_total(14)) / (1+exp(vec_A_total(14))),
                                exp(vec_A_total(15)) / (1+exp(vec_A_total(15))),
                                exp(vec_A_total(16)) / (1+exp(vec_A_total(16))),
                                exp(vec_A_total(17)) / (1+exp(vec_A_total(17))),
                                exp(vec_A_total(18)) / (1+exp(vec_A_total(18))),
                                exp(vec_A_total(19)) / (1+exp(vec_A_total(19)))};
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); // THREE STATE

    arma::vec vec_beta = par.elem(vec_beta_ind - 1);

    // The prior info for sigma_upsilon ----------------------------------------
    int nu_Upsilon = 100;
    
    // arma::vec scalar_mult2(20, arma::fill::ones); // THREE STATE
    arma::vec scalar_mult2 = {  9,  2,  2,  2,  2, 
                              400, 16, 16, 16, 16, 
                              400, 16, 16, 16, 16, 
                                9,  2,  2,  2,  2};
    scalar_mult2 = (nu_Upsilon - 20 - 1) * scalar_mult2;
    
    arma::mat Psi_Upsilon = arma::diagmat(scalar_mult2);

    arma::uvec vec_sigma_upsilon_ind = par_index(2);

    arma::uvec vec_alpha_tilde_ind = par_index(1);
    arma::vec vec_alpha_tilde = par.elem(vec_alpha_tilde_ind - 1);

    arma::vec eids = Y.col(0);

    arma::field<arma::mat> in_V_main(EIDs.n_elem);
    arma::field<arma::mat> in_inv_W_main(EIDs.n_elem);
    arma::field<arma::mat> in_Upsilon_cov_main(EIDs.n_elem);
    // arma::field<arma::mat> in_Upsilon_omega(EIDs.n_elem);

    // omp_set_num_threads(n_cores);
    // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;

        arma::vec b_i = B(ii);
        arma::vec vec_alpha_i = A(ii);
        arma::vec vec_omega_i = W(ii);
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1, 4);
        Y_i = Y_i.t();
        
        arma::field<arma::mat> Dn_alpha_full = Dn(ii);
        arma::field<arma::mat> Dn_omega_full = Dn_omega(ii);
        arma::field<arma::mat> Xn_full = Xn(ii);
        
        arma::mat interm_W(4,4,arma::fill::zeros);
        arma::vec interm_V(4,arma::fill::zeros);
        for(int k=1; k < n_i; k++) {
            arma::mat X_k = Xn_full(k);
            arma::mat X_k_1 = Xn_full(k-1);
            arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
            arma::mat A_1_k = arma::diagmat(vec_A_k);
            arma::mat diff_hold_x = X_k - A_1_k * X_k_1;
            interm_W = interm_W + diff_hold_x.t() * invR * diff_hold_x;
            
            arma::mat D_k = Dn_alpha_full(k);
            arma::mat D_k_1 = Dn_alpha_full(k-1);
            arma::mat diff_hold_d = D_k - A_1_k * D_k_1;
            
            arma::mat D_k_o = Dn_omega_full(k);
            arma::mat D_k_o_1 = Dn_omega_full(k-1);
            arma::mat diff_hold_d_o = D_k_o - A_1_k * D_k_o_1;
            
            arma::vec m_k = Y_i.col(k) - A_1_k * Y_i.col(k-1) 
                                - diff_hold_d * vec_alpha_i - diff_hold_d_o * vec_omega_i;
            interm_V = interm_V + diff_hold_x.t() * invR * m_k;
        }
        
        
        arma::vec vec_A = A_all_state.col(b_i(0) - 1);
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
        
        arma::mat V_i = Xn_full(0).t() * inv_Gamma * 
                (Y_i.col(0) - Dn_alpha_full(0)*vec_alpha_i - Dn_omega_full(0)*vec_omega_i) + interm_V;
        
        arma::mat hold2 = vec_alpha_i - vec_alpha_tilde;
        arma::mat in_Upsilon_cov = hold2 * hold2.t();
        
        in_V_main(ii) = V_i;
        in_inv_W_main(ii) = W_i_inv;
        in_Upsilon_cov_main(ii) = in_Upsilon_cov;
    }

    arma::mat sum_in_V = in_V_main(0);
    arma::mat sum_in_inv_W = in_inv_W_main(0);
    arma::mat sum_in_Upsilon_cov = in_Upsilon_cov_main(0);
    for(int ii = 1; ii < EIDs.n_elem; ii++) {
        sum_in_V = sum_in_V + in_V_main(ii);
        sum_in_inv_W = sum_in_inv_W + in_inv_W_main(ii);
        sum_in_Upsilon_cov = sum_in_Upsilon_cov + in_Upsilon_cov_main(ii);
    }

    arma::mat V = inv_Sigma_beta * vec_beta_0 + sum_in_V;

    arma::mat inv_W = inv_Sigma_beta + sum_in_inv_W;
    arma::mat W_b = arma::inv_sympd(inv_W);

    arma::mat Upsilon_cov = Psi_Upsilon + sum_in_Upsilon_cov;

    int n_sub = EIDs.n_elem;

    par.elem(vec_beta_ind - 1) = arma::mvnrnd(W_b * V, W_b);
    par.elem(vec_sigma_upsilon_ind - 1) = arma::vectorise(riwish(nu_Upsilon + n_sub, Upsilon_cov));

    return par;
}

// [[Rcpp::export]]
arma::mat update_Y_i_cpp( const arma::vec &EIDs, const arma::vec &par, 
                          const arma::field<arma::uvec> &par_index, 
                          const arma::field <arma::vec> &A, arma::mat Y,
                          arma::field <arma::field<arma::mat>> &Dn, 
                          const arma::field <arma::field<arma::mat>> &Xn, const arma::mat &otype,
                          const arma::field<arma::field<arma::mat>> &Dn_omega,
                          const arma::field <arma::vec> &W, arma::field <arma::vec> &B,
                          int n_cores) {

    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta, 
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID

    arma::mat newY(Y.n_rows, 4); 
    arma::vec eids = Y.col(0);
    
    arma::vec vec_beta = par.elem(par_index(0) - 1);
    
    arma::vec vec_R = par.elem(par_index(4) - 1);
    arma::mat R = arma::reshape(vec_R, 4, 4);
    arma::mat invR = arma::inv_sympd(R);

    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = { exp(vec_A_total(0)) / (1+exp(vec_A_total(0))),
                                exp(vec_A_total(1)) / (1+exp(vec_A_total(1))),
                                exp(vec_A_total(2)) / (1+exp(vec_A_total(2))),
                                exp(vec_A_total(3)) / (1+exp(vec_A_total(3))),
                                exp(vec_A_total(4)) / (1+exp(vec_A_total(4))),
                                exp(vec_A_total(5)) / (1+exp(vec_A_total(5))),
                                exp(vec_A_total(6)) / (1+exp(vec_A_total(6))),
                                exp(vec_A_total(7)) / (1+exp(vec_A_total(7))),
                                exp(vec_A_total(8)) / (1+exp(vec_A_total(8))),
                                exp(vec_A_total(9)) / (1+exp(vec_A_total(9))),
                                exp(vec_A_total(10)) / (1+exp(vec_A_total(10))),
                                exp(vec_A_total(11)) / (1+exp(vec_A_total(11))),
                                exp(vec_A_total(12)) / (1+exp(vec_A_total(12))),
                                exp(vec_A_total(13)) / (1+exp(vec_A_total(13))),
                                exp(vec_A_total(14)) / (1+exp(vec_A_total(14))),
                                exp(vec_A_total(15)) / (1+exp(vec_A_total(15))),
                                exp(vec_A_total(16)) / (1+exp(vec_A_total(16))),
                                exp(vec_A_total(17)) / (1+exp(vec_A_total(17))),
                                exp(vec_A_total(18)) / (1+exp(vec_A_total(18))),
                                exp(vec_A_total(19)) / (1+exp(vec_A_total(19)))};
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); // THREE STATE
    
      // omp_set_num_threads(n_cores);
      // # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {	
        
        int i = EIDs(ii);
        
        arma::uvec sub_ind = arma::find(eids == i);

        arma::vec b_i = B(ii);

        arma::field<arma::mat> Dn_ii = Dn(ii);
        arma::field<arma::mat> Dn_omega_ii = Dn_omega(ii);
        arma::field<arma::mat> Xn_ii = Xn(ii);
        arma::vec vec_alpha_ii = A(ii);
        arma::vec vec_omega_ii = W(ii);

        // Index of observed versus missing data
        // 1 = observed, 0 = missing
        arma::mat otype_i = otype.rows(sub_ind);
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1,4);
        Y_i = Y_i.t();
        otype_i = otype_i.t();
        arma::mat Y_i_new = Y_i;

        for(int k = 0; k < Y_i.n_cols; k++) {
            if(all(otype_i.col(k) == 1)) {
                Y_i_new.col(k) = Y_i.col(k);
            } else {
                if(k == 0) {
                    arma::vec vec_A = A_all_state.col(b_i(k) - 1);
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
                    
                    arma::vec vec_A_p1 = A_all_state.col(b_i(k+1) - 1);
                    arma::mat A_p1 = arma::diagmat(vec_A_p1);
                    
                    arma::vec nu_k = Dn_ii(k) * vec_alpha_ii + Dn_omega_ii(k) * vec_omega_ii + Xn_ii(k) * vec_beta;
                    arma::vec nu_k_p1 = Dn_ii(k+1) * vec_alpha_ii + Dn_omega_ii(k+1) * vec_omega_ii + Xn_ii(k+1) * vec_beta;
                    
                    arma::vec y_val_kp1 = Y_i_new.col(k+1);

                    arma::mat inv_W_i = inv_Gamma + A_p1.t() * invR * A_p1;
                    arma::mat W_i = inv(inv_W_i);
                    arma::vec V_i = inv_Gamma*nu_k + A_p1.t()*invR*(y_val_kp1 - nu_k_p1 + A_p1*nu_k);

                    arma::vec y_i_mean = W_i * V_i;

                    arma::vec new_value = arma::mvnrnd(y_i_mean, W_i, 1);
                    arma::vec update_value = Y_i_new.col(k);
                    arma::uvec ind_replace = arma::find(otype_i.col(k) == 0);
                    update_value.elem(ind_replace) = new_value.elem(ind_replace);
                    
                    // Prevent negatives
                    int count_while_loop = 0;
                    int count_while_loop_big = 0;
                    while(arma::any(update_value <= 0)) {
                            new_value = arma::mvnrnd(y_i_mean, W_i, 1);
                            update_value = Y_i_new.col(k);
                            update_value.elem(ind_replace) = new_value.elem(ind_replace);

                            count_while_loop += 1;
                            if(count_while_loop > 10000) {
                                count_while_loop_big += 1;
                                Rcpp::Rcout << "stuck in impute, i = " << ii << ", " << count_while_loop_big << std::endl;
                                count_while_loop = 0;
                            }
                            if(count_while_loop_big > 1000) {
                                break;
                            }
                    }

                    Y_i_new.col(k) = update_value;
                } else if(k == Y_i.n_cols - 1) {
                    
                    arma::vec vec_A = A_all_state.col(b_i(k) - 1);
                    arma::mat A_k = arma::diagmat(vec_A);
                    
                    arma::vec nu_k = Dn_ii(k) * vec_alpha_ii + Dn_omega_ii(k) * vec_omega_ii + Xn_ii(k) * vec_beta;
                    arma::vec nu_k_m1 = Dn_ii(k-1) * vec_alpha_ii + Dn_omega_ii(k-1) * vec_omega_ii +  Xn_ii(k-1) * vec_beta;
                    
                    arma::vec y_val_km1 = Y_i_new.col(k-1);

                    arma::vec y_i_mean = nu_k + A_k * (y_val_km1 - nu_k_m1);

                    arma::vec new_value = arma::mvnrnd(y_i_mean, R, 1);
                    arma::vec update_value = Y_i_new.col(k);
                    arma::uvec ind_replace = arma::find(otype_i.col(k) == 0);
                    update_value.elem(ind_replace) = new_value.elem(ind_replace);
                    
                    // Prevent negatives
                    int count_while_loop = 0;
                    int count_while_loop_big = 0;
                    while(arma::any(update_value <= 0)) {
                        new_value = arma::mvnrnd(y_i_mean, R, 1);
                        update_value = Y_i_new.col(k);
                        update_value.elem(ind_replace) = new_value.elem(ind_replace);

                        count_while_loop += 1;
                        if(count_while_loop > 10000) {
                            count_while_loop_big += 1;
                            Rcpp::Rcout << "stuck in impute, i = " << ii << ", " << count_while_loop_big << std::endl;
                            count_while_loop = 0;
                        }
                        if(count_while_loop_big > 1000) {
                            break;
                        }
                    }

                    Y_i_new.col(k) = update_value;
                } else {

                    arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
                    arma::vec vec_A_p1 = A_all_state.col(b_i(k+1) - 1);
                    arma::mat A_k = arma::diagmat(vec_A_k);                    
                    arma::mat A_p1 = arma::diagmat(vec_A_p1);

                    arma::vec nu_k    = Dn_ii(k) * vec_alpha_ii + Dn_omega_ii(k) * vec_omega_ii + Xn_ii(k) * vec_beta;
                    arma::vec nu_k_m1 = Dn_ii(k-1) * vec_alpha_ii + Dn_omega_ii(k-1) * vec_omega_ii + Xn_ii(k-1) * vec_beta;
                    arma::vec nu_k_p1 = Dn_ii(k+1) * vec_alpha_ii + Dn_omega_ii(k+1) * vec_omega_ii + Xn_ii(k+1) * vec_beta;
                    
                    arma::vec y_val_km1 = Y_i_new.col(k-1);
                    arma::vec y_val_kp1 = Y_i_new.col(k+1);

                    arma::mat inv_W_i = invR + A_p1.t() * invR * A_p1;
                    arma::mat W_i = inv(inv_W_i);

                    arma::vec V_i = invR * (nu_k + A_k * (y_val_km1 - nu_k_m1)) + 
                                        A_p1.t() * invR * (y_val_kp1 - nu_k_p1 + A_p1 * nu_k);

                    arma::vec y_i_mean = W_i * V_i;

                    arma::vec new_value = arma::mvnrnd(y_i_mean, W_i, 1);
                    arma::vec update_value = Y_i_new.col(k);
                    arma::uvec ind_replace = arma::find(otype_i.col(k) == 0);
                    update_value.elem(ind_replace) = new_value.elem(ind_replace);
                    
                    // Prevent negatives
                    int count_while_loop = 0;
                    int count_while_loop_big = 0;
                    while(arma::any(update_value <= 0)) {
                        new_value = arma::mvnrnd(y_i_mean, W_i, 1);
                        update_value = Y_i_new.col(k);
                        update_value.elem(ind_replace) = new_value.elem(ind_replace);

                        count_while_loop += 1;
                        if(count_while_loop > 10000) {
                            count_while_loop_big += 1;
                            Rcpp::Rcout << "stuck in impute, i = " << ii << ", " << count_while_loop_big << std::endl;
                            count_while_loop = 0;
                        }
                        if(count_while_loop_big > 1000) {
                            break;
                        }
                    }

                    Y_i_new.col(k) = update_value;
                }
            }
        }

        Y_i_new = Y_i_new.t();
        newY.rows(sub_ind) = Y_i_new;
    }
    
    Y.cols(1,4) = newY;
    return Y;
}

// [[Rcpp::export]]
Rcpp::List proposal_R_cpp(const int nu_R, const arma::mat psi_R, 
                          const arma::mat &Y, arma::field<arma::field<arma::mat>> &Dn, 
                          const arma::field<arma::field<arma::mat>> &Xn, const arma::field <arma::vec> A, 
                          const arma::vec par, const arma::field<arma::uvec> par_index, 
                          const arma::vec EIDs, arma::field <arma::vec> B,
                          const arma::field<arma::field<arma::mat>> Dn_omega,
                          const arma::field <arma::vec> W){
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    arma::vec eids = Y.col(0);
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = { exp(vec_A_total(0)) / (1+exp(vec_A_total(0))),
                                exp(vec_A_total(1)) / (1+exp(vec_A_total(1))),
                                exp(vec_A_total(2)) / (1+exp(vec_A_total(2))),
                                exp(vec_A_total(3)) / (1+exp(vec_A_total(3))),
                                exp(vec_A_total(4)) / (1+exp(vec_A_total(4))),
                                exp(vec_A_total(5)) / (1+exp(vec_A_total(5))),
                                exp(vec_A_total(6)) / (1+exp(vec_A_total(6))),
                                exp(vec_A_total(7)) / (1+exp(vec_A_total(7))),
                                exp(vec_A_total(8)) / (1+exp(vec_A_total(8))),
                                exp(vec_A_total(9)) / (1+exp(vec_A_total(9))),
                                exp(vec_A_total(10)) / (1+exp(vec_A_total(10))),
                                exp(vec_A_total(11)) / (1+exp(vec_A_total(11))),
                                exp(vec_A_total(12)) / (1+exp(vec_A_total(12))),
                                exp(vec_A_total(13)) / (1+exp(vec_A_total(13))),
                                exp(vec_A_total(14)) / (1+exp(vec_A_total(14))),
                                exp(vec_A_total(15)) / (1+exp(vec_A_total(15))),
                                exp(vec_A_total(16)) / (1+exp(vec_A_total(16))),
                                exp(vec_A_total(17)) / (1+exp(vec_A_total(17))),
                                exp(vec_A_total(18)) / (1+exp(vec_A_total(18))),
                                exp(vec_A_total(19)) / (1+exp(vec_A_total(19)))};
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); // THREE STATE

    arma::mat psi_prop_R_interm(4, 4, arma::fill::zeros);

    for (int ii = 0; ii < EIDs.n_elem; ii++) {

        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1,4);
        Y_i = Y_i.t();
        arma::vec vec_Y_i = arma::vectorise(Y_i);
        
        arma::vec b_i = B(ii);
        arma::vec vec_alpha_i = A(ii);
        arma::vec vec_omega_i = W(ii);
        arma::vec vec_beta = par.elem(par_index(0) - 1);
        arma::field<arma::mat> Xn_i = Xn(ii);
        arma::field<arma::mat> Dn_i = Dn(ii);
        arma::field<arma::mat> Dn_omega_i = Dn_omega(ii);
        
        arma::mat M(4, Y_i.n_cols - 1, arma::fill::zeros);
        for(int k = 1; k < Y_i.n_cols; k++) {
            arma::vec nu_k = Dn_i(k) * vec_alpha_i + Dn_omega_i(k) * vec_omega_i + Xn_i(k) * vec_beta;
            arma::vec nu_k_1 = Dn_i(k-1) * vec_alpha_i + Dn_omega_i(k-1) * vec_omega_i + Xn_i(k-1) * vec_beta;
            arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
            arma::mat A_1_k = arma::diagmat(vec_A_k);
            
            arma::vec m_temp = nu_k + A_1_k * (Y_i.col(k-1) - nu_k_1);
            M.col(k-1) = m_temp;
        }
        
        arma::mat Y_i_1 = Y_i.cols(1, Y_i.n_cols - 1);
        
        arma::mat hold = Y_i_1 - M;
        arma::mat hold2 = hold * hold.t();

        psi_prop_R_interm += hold2;
    }

    arma::mat psi_prop_R = psi_prop_R_interm + psi_R;
    int nu_prop_R = Y.n_rows + nu_R - EIDs.n_elem;

    List nu_psi_R = List::create(psi_prop_R, nu_prop_R);
    return nu_psi_R;
}

// [[Rcpp::export]]
Rcpp::List proposal_R_cpp_new(const int nu_R, const arma::mat psi_R, arma::mat curr_R,
                              const arma::mat &Y, arma::field<arma::field<arma::mat>> &Dn, 
                              const arma::field<arma::field<arma::mat>> &Xn, 
                              const arma::field <arma::vec> A, const arma::vec par, 
                              const arma::field<arma::uvec> par_index, 
                              const arma::vec EIDs, arma::field <arma::vec> B,
                              const arma::field<arma::field<arma::mat>> Dn_omega,
                              const arma::field <arma::vec> W){
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    arma::vec eids = Y.col(0);
    arma::vec vec_A_total = par.elem(par_index(3) - 1);
    arma::vec vec_A_scale = { exp(vec_A_total(0)) / (1+exp(vec_A_total(0))),
                                exp(vec_A_total(1)) / (1+exp(vec_A_total(1))),
                                exp(vec_A_total(2)) / (1+exp(vec_A_total(2))),
                                exp(vec_A_total(3)) / (1+exp(vec_A_total(3))),
                                exp(vec_A_total(4)) / (1+exp(vec_A_total(4))),
                                exp(vec_A_total(5)) / (1+exp(vec_A_total(5))),
                                exp(vec_A_total(6)) / (1+exp(vec_A_total(6))),
                                exp(vec_A_total(7)) / (1+exp(vec_A_total(7))),
                                exp(vec_A_total(8)) / (1+exp(vec_A_total(8))),
                                exp(vec_A_total(9)) / (1+exp(vec_A_total(9))),
                                exp(vec_A_total(10)) / (1+exp(vec_A_total(10))),
                                exp(vec_A_total(11)) / (1+exp(vec_A_total(11))),
                                exp(vec_A_total(12)) / (1+exp(vec_A_total(12))),
                                exp(vec_A_total(13)) / (1+exp(vec_A_total(13))),
                                exp(vec_A_total(14)) / (1+exp(vec_A_total(14))),
                                exp(vec_A_total(15)) / (1+exp(vec_A_total(15))),
                                exp(vec_A_total(16)) / (1+exp(vec_A_total(16))),
                                exp(vec_A_total(17)) / (1+exp(vec_A_total(17))),
                                exp(vec_A_total(18)) / (1+exp(vec_A_total(18))),
                                exp(vec_A_total(19)) / (1+exp(vec_A_total(19)))};
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5); // THREE STATE

    arma::mat psi_prop_R_interm(4, 4, arma::fill::zeros);

    for (int ii = 0; ii < EIDs.n_elem; ii++) {

        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat Y_i = Y_temp.cols(1,4);
        Y_i = Y_i.t();
        
        arma::vec b_i = B(ii);
        arma::vec vec_alpha_i = A(ii);
        arma::vec vec_omega_i = W(ii);
        arma::vec vec_beta = par.elem(par_index(0) - 1);
        arma::field<arma::mat> Xn_i = Xn(ii);
        arma::field<arma::mat> Dn_i = Dn(ii);
        arma::field<arma::mat> Dn_omega_i = Dn_omega(ii);
        
        for(int k = 0; k < Y_i.n_cols; k++) {
            if(k == 0) {
                arma::vec vec_A = A_all_state.col(b_i(k) - 1);
                arma::mat Gamma     = {{curr_R(0,0) / (1 - vec_A(0) * vec_A(0)), 
                                        curr_R(0,1) / (1 - vec_A(0) * vec_A(1)), 
                                        curr_R(0,2) / (1 - vec_A(0) * vec_A(2)), 
                                        curr_R(0,3) / (1 - vec_A(0) * vec_A(3))},
                                       {curr_R(1,0) / (1 - vec_A(1) * vec_A(0)), 
                                        curr_R(1,1) / (1 - vec_A(1) * vec_A(1)), 
                                        curr_R(1,2) / (1 - vec_A(1) * vec_A(2)), 
                                        curr_R(1,3) / (1 - vec_A(0) * vec_A(3))},
                                       {curr_R(2,0) / (1 - vec_A(2) * vec_A(0)), 
                                        curr_R(2,1) / (1 - vec_A(2) * vec_A(1)), 
                                        curr_R(2,2) / (1 - vec_A(2) * vec_A(2)), 
                                        curr_R(2,3) / (1 - vec_A(0) * vec_A(3))},
                                       {curr_R(3,0) / (1 - vec_A(3) * vec_A(0)), 
                                        curr_R(3,1) / (1 - vec_A(3) * vec_A(1)), 
                                        curr_R(3,2) / (1 - vec_A(3) * vec_A(2)), 
                                        curr_R(3,3) / (1 - vec_A(0) * vec_A(3))}};
                arma::mat inv_Gamma = arma::inv_sympd(Gamma);
                
                arma::mat curr_R_sqrt = arma::sqrtmat_sympd(curr_R);
                arma::mat inv_gamma_sqrt = arma::sqrtmat_sympd(inv_Gamma);
                
                arma::vec nu_1 = Dn_i(k) * vec_alpha_i + Dn_omega_i(k) * vec_omega_i + Xn_i(k) * vec_beta;
                arma::vec y_diff = Y_i.col(k) - nu_1;
                y_diff = curr_R_sqrt * inv_gamma_sqrt * y_diff;
                
                psi_prop_R_interm = psi_prop_R_interm + y_diff * y_diff.t();
            } else {
                arma::vec nu_k = Dn_i(k) * vec_alpha_i + Dn_omega_i(k) * vec_omega_i + Xn_i(k) * vec_beta;
                arma::vec nu_k_1 = Dn_i(k-1) * vec_alpha_i + Dn_omega_i(k-1) * vec_omega_i + Xn_i(k-1) * vec_beta;
                arma::vec vec_A_k = A_all_state.col(b_i(k) - 1);
                arma::mat A_1_k = arma::diagmat(vec_A_k);
                
                arma::vec m_temp = nu_k + A_1_k * (Y_i.col(k-1) - nu_k_1);
                
                arma::vec y_diff = Y_i.col(k) - m_temp;
                psi_prop_R_interm = psi_prop_R_interm + y_diff * y_diff.t();
            }
        }
    }

    arma::mat psi_prop_R = psi_prop_R_interm + psi_R;
    int nu_prop_R = Y.n_rows + nu_R;

    List nu_psi_R = List::create(psi_prop_R, nu_prop_R);
    return nu_psi_R;
}

// [[Rcpp::export]]
void test_fnc() {
    
    int nu_R = 100;
    
    // arma::mat Psi_R(4,4,arma::fill::eye);
    arma::vec scalar_vec_R = {4.58, 98.2, 101.3, 7.6};
    arma::mat Psi_R = arma::diagmat(scalar_vec_R);
    Psi_R = (nu_R - 4 - 1) * Psi_R;
    
    Rcpp::Rcout << Psi_R << std::endl;
    
    for(int i = 0; i < 10; i++) {
        int a = 0;
        Rcpp::Rcout << i << std::endl;
        
        while(a < 5) {
            a += 1;
            
            if(a == 3) {
                break;
            }
        }
        Rcpp::Rcout << "a " << a << std::endl;
    }
    // int N = 5;
    // Rcpp::Rcout << "Case (c) Full" << std::endl;
    // for(int w=0; w < N; w++) {
    //   Rcpp::Rcout << "() -> () -> " << w+1 << std::endl;
    //   Rcpp::Rcout << Omega_List_GLOBAL(0)(w) << std::endl;
    // }
    // 
    // Rcpp::Rcout << "Case (b) Full" << std::endl;
    // for(int i = 0; i < N; i++) {
    //   for(int j = 0; j < N; j++) {
    //     Rcpp::Rcout << i+1 << "-->" << j+1 << std::endl;
    //     Rcpp::Rcout << Omega_List_GLOBAL(1)(i, j) << std::endl;
    //   }
    // }
    // 
    // Rcpp::Rcout << "Case (a) Full" << std::endl;
    // for(int w=0; w < N; w++) {
    //   Rcpp::Rcout << w + 1 << " -> () -> ()" << std::endl;
    //   Rcpp::Rcout << Omega_List_GLOBAL(2)(w) << std::endl;
    // }
    // 
    // 
    // Rcpp::Rcout << "Case (c) Sub" << std::endl;
    // for (int w = 0; w < N; w++) {
    //   Rcpp::Rcout << "() -> () -> " << w + 1 << std::endl;
    //   Rcpp::Rcout << Omega_List_GLOBAL_sub(0)(w) << std::endl;
    // }
    // 
    // Rcpp::Rcout << "Case (b) Sub" << std::endl;
    // for (int i = 0; i < N; i++) {
    //   for (int j = 0; j < N; j++) {
    //     Rcpp::Rcout << i + 1 << "-->" << j + 1 << std::endl;
    //     Rcpp::Rcout << Omega_List_GLOBAL_sub(1)(i, j) << std::endl;
    //   }
    // }
    // 
    // Rcpp::Rcout << "Case (a) Sub" << std::endl;
    // for (int w = 0; w < N; w++) {
    //   Rcpp::Rcout << w + 1 << " -> () -> ()" << std::endl;
    //   Rcpp::Rcout << Omega_List_GLOBAL_sub(2)(w) << std::endl;
    // }
    // 
    // bool valid_prop = false;
    // 
    // int clinic_rule = 0;
    // int rbc_rule = 1;
    // arma::vec pr_B = {1,1,1,1,1,1,2,2,2,2,3,3,3,3,3,1,1,1,2,3};
    // arma::vec old_B = pr_B;
    // old_B(4) = 3;
    // arma::vec bleed_indicator = {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0};
    // 
    // int pos_bleed = arma::as_scalar(arma::find(bleed_indicator == 1));
    // Rcpp::Rcout << pos_bleed << std::endl;
    // 
    // bool b_i_rule = arma::any(arma::vectorise(pr_B)==2);
    // Rcpp::Rcout << b_i_rule << std::endl;
    // 
    // arma::uvec b_i_time = arma::find(pr_B == 2);
    // Rcpp::Rcout << b_i_time << std::endl;
    // 
    // bool b_i_time_rule = arma::any(b_i_time <= pos_bleed);
    // Rcpp::Rcout << b_i_time_rule << std::endl;
    // 
    // if(clinic_rule >= 0) {
    //     if (clinic_rule == 1) {
    //         if(b_i_rule) {valid_prop = true;}
    //     } else {
    //         if (rbc_rule == 0 || (rbc_rule == 1 && b_i_rule)) {valid_prop = true;}
    //     }
    // } else {
    //     valid_prop = true; // evaluate likelihood anyways because S1&S3
    // }
    // 
    // if(arma::accu(pr_B == old_B) == pr_B.n_elem) {
    //     Rcpp::Rcout << "same" << std::endl;
    // } 
}
