#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec state_prob_dist(const int k, const int n_i, int t_pt_length, 
                          const arma::vec &par, const arma::field<arma::uvec> &par_index, 
                          const arma::mat &y_i, const arma::mat &z_i, const arma::vec &b_i,
                          const arma::vec &alpha_i, const arma::field<arma::mat> &x_i,
                          const arma::field<arma::mat> &Dn_omega, const arma::vec &w_i,
                          arma::mat omega_set) {
    
    // Parameter initialization ------------------------------------------------
    arma::mat vec_beta = par.elem(par_index(0) - 1);
    
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
    arma::mat A_all_state = arma::reshape(vec_A_scale, 4, 5);
    
    arma::vec vec_zeta_content = par.elem(par_index(5) - 1);
    arma::mat zeta = arma::reshape(vec_zeta_content, 2, 12); 
    
    arma::vec vec_init_content = par.elem(par_index(6) - 1);
    arma::vec init_logit = {1, exp(vec_init_content(0)), exp(vec_init_content(1)),
                            exp(vec_init_content(2)), exp(vec_init_content(3))}; 
    arma::vec P_init = init_logit / arma::accu(init_logit); 
    
    arma::vec prob_dist(omega_set.n_rows, arma::fill::ones);
    
    // ------------------------------------------------------------------------
    for(int j = 0; j < omega_set.n_rows; j++) {
        
        // Initialize the new state sequence and the points to evaluate likelihood
        arma::vec t_pts;
        arma::vec ss_j = b_i;
        if (k == 1) {
            // () -> () -> 1-5
            ss_j.rows(0, omega_set.n_cols - 1) = omega_set.row(j).t();
            t_pts = arma::linspace(0, n_i - 1, n_i);
        
        } else if (k <= n_i - t_pt_length) {
            // 1-5 -> () -> () -> 1-5
            ss_j.rows(k - 1, k + t_pt_length - 2) = omega_set.row(j).t();
            t_pts = arma::linspace(k - 1, n_i - 1, n_i - k + 1);
            
        } else if (k == n_i - t_pt_length + 1) {
            // 1-5 -> () -> ()
            ss_j.rows(k - 1, k + t_pt_length - 2) = omega_set.row(j).t();
            t_pts = arma::linspace(k - 1, n_i - 1, n_i - k + 1);
        }
        
        // Create new design matrix -------------------------------------------
        arma::field<arma::mat> Dn_ss_j(ss_j.n_elem);
        
        arma::vec twos(ss_j.n_elem, arma::fill::zeros);
        arma::vec threes = twos; 
        arma::vec fours = twos;
        arma::vec fives = twos;
        
        twos.elem(arma::find(ss_j == 2)) += 1;
        threes.elem(arma::find(ss_j == 3)) += 1; 
        fours.elem(arma::find(ss_j == 4)) += 1;
        fives.elem(arma::find(ss_j == 5)) += 1;
        
        arma::vec ones(ss_j.n_elem, arma::fill::ones);
        
        arma::mat bigB = arma::join_rows(ones, arma::cumsum(twos));
        bigB = arma::join_rows(bigB, arma::cumsum(threes));
        bigB = arma::join_rows(bigB, arma::cumsum(fours));
        bigB = arma::join_rows(bigB, arma::cumsum(fives));
        
        arma::mat I = arma::eye(4,4);
        for(int jj = 0; jj < n_i; jj++) {
            Dn_ss_j(jj) = arma::kron(I, bigB.row(jj));
        } 
        
        // Likelihood computations ---------------------------------------------
        double like_comp = 0;
        for(int jj = 0; jj < t_pts.n_elem; jj++) {
            
            int t_j = t_pts(jj);
            
            // (1) Transition probabilities ------------------------------------
            if(jj < 3) {
                if(t_j == 0) {
                    like_comp = like_comp + log(P_init(ss_j(t_j) - 1));
                } else{
                    // State space component
                    double q1_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(0));
                    double q1 = exp(q1_sub);
                    double q2_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(1));
                    double q2 = exp(q2_sub);
                    double q3_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(2));
                    double q3 = exp(q3_sub);
                    double q4_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(3));
                    double q4 = exp(q4_sub);
                    
                    double q5_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(4));
                    double q5 = exp(q5_sub);
                    double q6_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(5));
                    double q6 = exp(q6_sub);
                    double q7_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(6));
                    double q7 = exp(q7_sub);
                    double q8_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(7));
                    double q8 = exp(q8_sub);
                    
                    double q9_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(8));
                    double q9 = exp(q9_sub);
                    double q10_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(9));
                    double q10 = exp(q10_sub);
                    double q11_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(10));
                    double q11 = exp(q11_sub);
                    double q12_sub = arma::as_scalar(z_i.row(t_j) * zeta.col(11));
                    double q12 = exp(q12_sub);
                    
                    arma::mat Q = { {   1,   q1,  0,  q2,  0},
                                    {   0,    1, q3,  q4,  0},
                                    {  q5,   q6,  1,  q7,  0},
                                    {   0,   q8,  0,   1, q9},
                                    { q10,  q11,  0, q12,  1}}; 
                    
                    arma::vec q_row_sums = arma::sum(Q, 1);
                    arma::mat P_i = Q.each_col() / q_row_sums;
                    int b_k_1 = ss_j(t_j-1);
                    int b_k = ss_j(t_j);
                    
                    like_comp = like_comp + log(P_i(b_k_1 - 1, b_k - 1));
                }
            }
            
            // (2) Response outcome --------------------------------------------
            if(t_j == 0) {
                arma::vec vec_A = A_all_state.col(ss_j(t_j) - 1);
                
                arma::mat Gamma = { {R(0,0) / (1 - vec_A(0) * vec_A(0)), 
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
                
                arma::vec y_1 = y_i.col(t_j);
                arma::vec nu_1 = Dn_ss_j(t_j) * alpha_i 
                                + Dn_omega(t_j) * w_i 
                                + x_i(t_j) * vec_beta;
                arma::vec like_y = dmvnorm(y_1.t(), nu_1, Gamma, true);
                
                like_comp = like_comp + arma::as_scalar(like_y);
                
            } else {
                arma::vec vec_A = A_all_state.col(ss_j(t_j) - 1);
                arma::mat A_1 = arma::diagmat(vec_A);
                
                arma::vec y_k_1 = y_i.col(t_j - 1);
                arma::vec y_k = y_i.col(t_j);
                arma::vec nu_k_1 = Dn_ss_j(t_j - 1) * alpha_i 
                                    + Dn_omega(t_j - 1) * w_i 
                                    + x_i(t_j - 1) * vec_beta;
                arma::vec nu_k = Dn_ss_j(t_j) * alpha_i 
                                    + Dn_omega(t_j) * w_i 
                                    + x_i(t_j) * vec_beta;

                arma::vec mean_k = nu_k + A_1 * (y_k_1 - nu_k_1);
                arma::vec like_y_k = dmvnorm(y_k.t(), mean_k, R, true);
                
                like_comp = like_comp + arma::as_scalar(like_y_k);
            }
        }
        
        // Converting out of the "log" domain
        like_comp = exp(like_comp);
        
        prob_dist(j) = like_comp;
        // Rcpp::Rcout << ss_j.t() << std::endl;
        // Rcpp::Rcout << "prob = " << prob_dist(j) << std::endl;
    }
    
    prob_dist= (1/arma::accu(prob_dist)) * prob_dist;
    // Rcpp::Rcout << prob_dist << std::endl;
    
    return prob_dist;
}

// [[Rcpp::export]]
Rcpp::List update_b_i_up(const arma::vec EIDs, const arma::vec &par, 
                         const arma::field<arma::uvec> &par_index, 
                         const arma::field <arma::vec> &A, 
                         arma::field <arma::vec> &B, 
                         const arma::mat &Y, const arma::mat &z, 
                         arma::field<arma::field<arma::mat>> &Dn, 
                         const arma::field <arma::field<arma::mat>> &Xn, 
                         const arma::field<arma::field<arma::mat>> &Dn_omega, 
                         const arma::field <arma::vec> &W,
                         const arma::vec &bleed_indicator, int n_cores,
                         int t_pt_length) {
    
    // par_index KEY: (0) beta, (1) alpha_tilde, (2) sigma_upsilon, (3) vec_A, (4) R, (5) zeta,
    //                (6) init, (7) omega_tilde, (8) vec_upsilon_omega
    // Y key: (0) EID, (1) hemo, (2) hr, (3) map, (4) lactate, (5) RBC, (6) clinic
    // "i" is the numeric EID number
    // "ii" is the index of the EID
    
    // In previous iterations, t_pt_length = 2
    arma::vec eids = Y.col(0); 
    arma::vec rbc_rule_vec = Y.col(5);
    arma::vec clinic_rule_vec = Y.col(6); 
    
    arma::field<arma::vec> B_return(EIDs.n_elem);
    arma::field<arma::field<arma::mat>> Dn_return(EIDs.n_elem);
    
    omp_set_num_threads(n_cores);
    # pragma omp parallel for
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        // Subject-specific information ----------------------------------------
        int i = EIDs(ii);
        
        arma::uvec sub_ind = arma::find(eids == i);
        int n_i = sub_ind.n_elem;
        
        int rbc_rule = rbc_rule_vec(sub_ind.min());
        int clinic_rule = clinic_rule_vec(sub_ind.min());
        
        // Sub-setting fields --------------------------------------------------
        arma::vec B_temp = B(ii);
        arma::vec A_temp = A(ii);
        arma::vec W_temp = W(ii);
        arma::vec bleed_ind_i = bleed_indicator.elem(sub_ind);
        
        arma::field<arma::mat> Dn_temp = Dn(ii);
        arma::field<arma::mat> Dn_omega_temp = Dn_omega(ii);
        arma::field<arma::mat> Xn_temp = Xn(ii);
        
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat z_temp = z.rows(sub_ind);
        
        arma::mat y_i = Y_temp.cols(1, 4);
        y_i = y_i.t();
        
        // Looping through subject state space ---------------------------------
        for (int k = 0; k < n_i - (t_pt_length - 1); k++) {
            
            // All possible state transitions given current time point ---------
            arma::mat Omega_set;
            if (clinic_rule >= 0) {
                Omega_set = Omega_fun_cpp_new_multi(k + 1, n_i, B_temp, false, 
                                                    t_pt_length);
            } else { 
                Omega_set = Omega_fun_cpp_new_multi(k + 1, n_i, B_temp, true, 
                                                    t_pt_length);
            } 
            
            // Learn the proposal distribution ---------------------------------
            arma::vec ss_prob = state_prob_dist(k+1, n_i, t_pt_length, par,
                                                par_index, y_i, z_temp, B_temp,
                                                A_temp, Xn_temp, Dn_omega_temp,
                                                W_temp, Omega_set);
            
            arma::vec x_sample = arma::linspace(1, Omega_set.n_rows, Omega_set.n_rows);
            arma::vec row_ind = RcppArmadillo::sample(x_sample, 1, false, ss_prob);
            
            // Gibbs update ----------------------------------------------------
            B_temp.rows(k, k+t_pt_length-1) = Omega_set.row(row_ind(0)-1).t();
            
            arma::vec twos(B_temp.n_elem, arma::fill::zeros);
            arma::vec threes = twos; // THREE STATE
            arma::vec fours = twos;
            arma::vec fives = twos;
            
            twos.elem(arma::find(B_temp == 2)) += 1;
            threes.elem(arma::find(B_temp == 3)) += 1; // THREE STATE
            fours.elem(arma::find(B_temp == 4)) += 1;
            fives.elem(arma::find(B_temp == 5)) += 1;
            
            arma::vec ones(B_temp.n_elem, arma::fill::ones);
            
            arma::mat bigB = arma::join_rows(ones, arma::cumsum(twos));
            bigB = arma::join_rows(bigB, arma::cumsum(threes)); // THREE STATE
            bigB = arma::join_rows(bigB, arma::cumsum(fours));
            bigB = arma::join_rows(bigB, arma::cumsum(fives));
            
            arma::mat I = arma::eye(4,4);
            for(int jj = 0; jj < n_i; jj++) {
                Dn_temp(jj) = arma::kron(I, bigB.row(jj));
            } 
        }
        B_return(ii) = B_temp;
        Dn_return(ii) = Dn_temp;
    } 
    List B_Dn = List::create(B_return, Dn_return);
    
    return B_Dn;
} 
