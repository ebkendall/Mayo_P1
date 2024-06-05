#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

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

// Exported functions -------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List update_b_i_online(const arma::vec EIDs, const arma::vec &par, 
                          const arma::field<arma::uvec> &par_index, 
                          const arma::field <arma::vec> &A, 
                          arma::field <arma::vec> &B, 
                          const arma::mat &Y, const arma::mat &z, 
                          arma::field<arma::field<arma::mat>> &Dn, 
                          const arma::field <arma::field<arma::mat>> &Xn, 
                          const arma::field<arma::field<arma::mat>> &Dn_omega, 
                          const arma::field <arma::vec> &W, int n_cores,
                          const arma::vec &x) {
    
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
    
    for (int ii = 0; ii < EIDs.n_elem; ii++) {
        int i = EIDs(ii);
        arma::uvec sub_ind = arma::find(eids == i);
        
        int n_i = sub_ind.n_elem;
        
        // Subsetting fields
        arma::vec B_temp = B(ii);
        arma::vec A_temp = A(ii);
        arma::vec W_temp = W(ii);
        
        arma::field<arma::mat> Dn_temp = Dn(ii);
        arma::field<arma::mat> Dn_omega_temp = Dn_omega(ii);
        arma::field<arma::mat> Xn_temp = Xn(ii);
        
        // Subsetting the remaining data
        arma::mat Y_temp = Y.rows(sub_ind);
        arma::mat z_temp = z.rows(sub_ind);
        
        int clinic_rule = clinic_rule_vec(sub_ind.min());
        
        arma::vec n_RBC_ii = x.elem(sub_ind);
        int rbc_rule = 0;
        int twelve_hr_ind = 0;
        int twoFour_hr_ind = 0;
        
        if(n_i > 48) {
            twelve_hr_ind = n_i - 48;
        }
        if(n_i > 96) {
            twoFour_hr_ind = n_i - 96;
        }
        if(n_RBC_ii(n_i-1) - n_RBC_ii(twelve_hr_ind) >= 3) {
            rbc_rule = 1;
        }
        if(n_RBC_ii(n_i-1) - n_RBC_ii(twoFour_hr_ind) >= 6) {
            rbc_rule = 1;
        }
        
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
                if(b_i_rule) {
                    valid_prop = true;
                }
            } else if(clinic_rule == 0) {
                if(rbc_rule == 1) {
                    if(b_i_rule) {
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
    Rcpp::List B_Dn = Rcpp::List::create(B_return, Dn_return);
    
    return B_Dn;
}