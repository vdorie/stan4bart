functions {
  vector make_theta_L(int len_theta_L, int[] p, real dispersion,
                      vector tau, vector scale, vector zeta,
                      vector rho, vector z_T)
  {
    vector[len_theta_L] theta_L;
    int zeta_mark = 1;
    int rho_mark = 1;
    int z_T_mark = 1;
    int theta_L_mark = 1;

    // each of these is a diagonal block of the implicit Cholesky factor
    for (i in 1:size(p)) { 
      int nc = p[i];
      if (nc == 1) { // "block" is just a standard deviation
        theta_L[theta_L_mark] = tau[i] * scale[i] * dispersion;
        // unlike lme4, theta[theta_L_mark] includes the dispersion term in it
        theta_L_mark += 1;
      }
      else { // block is lower-triangular               
        matrix[nc,nc] T_i; 
        real std_dev;
        real T21;
        real trace_T_i = square(tau[i] * scale[i] * dispersion) * nc;
        vector[nc] pi = segment(zeta, zeta_mark, nc); // gamma(zeta | shape, 1)
        pi /= sum(pi);                            // thus dirichlet(pi | shape)
        
        // unlike lme4, T_i includes the dispersion term in it
        zeta_mark += nc;
        std_dev = sqrt(pi[1] * trace_T_i);
        T_i[1,1] = std_dev;
        
        // Put a correlation into T_i[2,1] and scale by std_dev
        std_dev = sqrt(pi[2] * trace_T_i);
        T21 = 2.0 * rho[rho_mark] - 1.0;
        rho_mark += 1;
        T_i[2,2] = std_dev * sqrt(1.0 - square(T21));
        T_i[2,1] = std_dev * T21;
        
        for (r in 2:(nc - 1)) { // scaled onion method to fill T_i
          int rp1 = r + 1;
          vector[r] T_row = segment(z_T, z_T_mark, r);
          real scale_factor = sqrt(rho[rho_mark] / dot_self(T_row)) * std_dev;
          z_T_mark += r;
          std_dev = sqrt(pi[rp1] * trace_T_i);
          for(c in 1:r) T_i[rp1,c] = T_row[c] * scale_factor;
          T_i[rp1,rp1] = sqrt(1.0 - rho[rho_mark]) * std_dev;
          rho_mark += 1;
        }
        
        // now vech T_i
        for (c in 1:nc) for (r in c:nc) {
          theta_L[theta_L_mark] = T_i[r,c];
          theta_L_mark += 1;
        }
      }
    }
    return theta_L;
  }
  
  vector make_b(vector z_b, vector theta_L, int[] p, int[] l)
  {
    vector[rows(z_b)] b;
    int b_mark = 1;
    int theta_L_mark = 1;
    for (i in 1:size(p)) {
      int nc = p[i];
      if (nc == 1) {
        real theta_L_start = theta_L[theta_L_mark];
        for (s in b_mark:(b_mark + l[i] - 1)) 
          b[s] = theta_L_start * z_b[s];
        b_mark += l[i];
        theta_L_mark += 1;
      }
      else {
        matrix[nc,nc] T_i = rep_matrix(0, nc, nc);
        for (c in 1:nc) {
          T_i[c,c] = theta_L[theta_L_mark];
          theta_L_mark += 1;
          for(r in (c+1):nc) {
            T_i[r,c] = theta_L[theta_L_mark];
            theta_L_mark += 1;
          }
        }
        for (j in 1:l[i]) {
          vector[nc] temp = T_i * segment(z_b, b_mark, nc);
          b_mark -= 1;
          for (s in 1:nc) b[b_mark + s] = temp[s];
          b_mark += nc + 1;
        }
      }
    }
    return b;
  }
  
  real decov_lp(vector z_b, vector z_T, vector rho, vector zeta, vector tau,
                real[] regularization, real[] delta, vector shape,
                int t, int[] p)
  {
    int pos_reg = 1;
    int pos_rho = 1;
    target += normal_lpdf(z_b | 0, 1);
    target += normal_lpdf(z_T | 0, 1);
    for (i in 1:t) if (p[i] > 1) {
      vector[p[i] - 1] shape1;
      vector[p[i] - 1] shape2;
      real nu = regularization[pos_reg] + 0.5 * (p[i] - 2);
      pos_reg += 1;
      shape1[1] = nu;
      shape2[1] = nu;
      for (j in 2:(p[i]-1)) {
        nu -= 0.5;
        shape1[j] = 0.5 * j;
        shape2[j] = nu;
      }
      target += beta_lpdf(rho[pos_rho:(pos_rho + p[i] - 2)] | shape1, shape2);
      pos_rho += p[i] - 1;
    }
    target += gamma_lpdf(zeta | delta, 1);
    target += gamma_lpdf(tau  | shape, 1);
    return target();
  }
  
  vector hs_prior(vector z_beta, real[] global, vector[] local, 
                  real global_prior_scale, real error_scale, real c2) {
    int K = rows(z_beta);
    vector[K] lambda = local[1] .* sqrt(local[2]);
    real tau = global[1] * sqrt(global[2]) * global_prior_scale * error_scale;
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt( c2 * lambda2 ./ (c2 + square(tau) * lambda2) );
    return z_beta .* lambda_tilde * tau;
  }
  
  vector hsplus_prior(vector z_beta, real[] global, vector[] local, 
                      real global_prior_scale, real error_scale, real c2) {
    int K = rows(z_beta);
    vector[K] lambda = local[1] .* sqrt(local[2]);
    vector[K] eta = local[3] .* sqrt(local[4]);
    real tau = global[1] * sqrt(global[2]) * global_prior_scale * error_scale;
    vector[K] lambda_eta2 = square(lambda .* eta);
    vector[K] lambda_tilde = sqrt( c2 * lambda_eta2 ./ 
                                 ( c2 + square(tau) * lambda_eta2) );
    return z_beta .* lambda_tilde * tau;
  }
  
  real CFt(real z, real df) {
    real z2 = square(z);
    real z3 = z2 * z;
    real z5 = z2 * z3;
    real z7 = z2 * z5;
    real z9 = z2 * z7;
    real df2 = square(df);
    real df3 = df2 * df;
    real df4 = df2 * df2;
    return z + (z3 + z) / (4 * df) + (5 * z5 + 16 * z3 + 3 * z) / (96 * df2)
           + (3 * z7 + 19 * z5 + 17 * z3 - 15 * z) / (384 * df3)
           + (79 * z9 + 776 * z7 + 1482 * z5 - 1920 * z3 - 945 * z) / (92160 * df4);
  }
  
  vector csr_matrix_times_vector2(int m, int n, vector w, 
                                  int[] v, int[] u, vector b);
  
  vector pw_gauss(vector y, vector eta, real sigma) {
    return -0.5 * log(6.283185307179586232 * sigma * sigma) - 0.5 * square(y - eta) / (sigma * sigma);
  }
}
data {
  // dimensions
  int<lower=0> N; // number of observations
  int<lower=0> K; // number of mean-parameters
  
  // data
  matrix[N,K] X;  // dense centered predictor matrix
    
  int<lower=0> len_y;
  real lb_y; // lower bound on y
  real<lower=lb_y> ub_y; // upper bound on y
  vector<lower=lb_y, upper=ub_y>[len_y] y;
  
  // intercept
  int<lower=0,upper=1> has_intercept;  // 1 = yes
  int<lower=0,upper=1> is_binary; // 1 = yes
  
  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus, 
  //   5 = laplace, 6 = lasso, 7 = product_normal
  int<lower=0,upper=7> prior_dist;
  int<lower=0,upper=2> prior_dist_for_intercept;
  
  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
  int<lower=0,upper=3> prior_dist_for_aux;
  
    
  int<lower=0,upper=1> has_weights;  // 0 = No, 1 = Yes
  vector[has_weights ? N : 0] weights;
  
  vector[N] offset_;
  
  // hyperparameter values are set to 0 if there is no prior
  vector<lower=0>[K] prior_scale;
  real<lower=0> prior_scale_for_intercept;
  real<lower=0> prior_scale_for_aux;
  
  vector[K] prior_mean;
  real prior_mean_for_intercept;
  real<lower=0> prior_mean_for_aux;

  vector<lower=0>[K] prior_df;
  real<lower=0> prior_df_for_intercept;
  real<lower=0> prior_df_for_aux;
  
  real<lower=0> global_prior_df;     // for hs priors only
  real<lower=0> global_prior_scale;  // for hs priors only
  real<lower=0> slab_df;     // for hs prior only
  real<lower=0> slab_scale;  // for hs prior only
  int<lower=2> num_normals[prior_dist == 7 ? K : 0];
    
  // glmer stuff, see table 3 of
  // https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  int<lower=0> t;               // num. terms (maybe 0) with a | in the glmer formula
  int<lower=1> p[t];            // num. variables on the LHS of each |
  int<lower=1> l[t];            // num. levels for the factor(s) on the RHS of each |
  int<lower=0> q;               // conceptually equals \sum_{i=1}^t p_i \times l_i
  int<lower=0> len_theta_L;     // length of the theta_L vector

  // hyperparameters for glmer stuff; if t > 0 priors are mandatory
  vector<lower=0>[t] shape; 
  vector<lower=0>[t] scale;
  int<lower=0> len_concentration;
  real<lower=0> concentration[len_concentration];
  int<lower=0> len_regularization;
  real<lower=0> regularization[len_regularization];
  
    
  int<lower=0> num_non_zero;  // number of non-zero elements in the Z matrix
  vector[num_non_zero] w;     // non-zero elements in the implicit Z matrix
  int<lower=0, upper=q-1> v[num_non_zero];               // column indices for w
  int<lower=0, upper=rows(w) + 1> u[N + 1];  // where the non-zeros start in each row
}
transformed data {
  int<lower=0> len_z_T = 0;
  int<lower=0> len_var_group = sum(p);
  int<lower=0> len_rho = sum(p) - t;
  int<lower=1> pos = 1;
  real<lower=0> delta[len_concentration];
  int<lower=0> hs;
  if (prior_dist <= 2) hs = 0;
  else if (prior_dist == 3) hs = 2;
  else if (prior_dist == 4) hs = 4;
  else hs = 0;
  
  for (i in 1:t) {
    if (p[i] > 1) {
      for (j in 1:p[i]) {
        delta[pos] = concentration[j];
        pos += 1;
      }
    }
    for (j in 3:p[i]) len_z_T += p[i] - 1;
  }
}
parameters {
  real<lower=negative_infinity(),upper=positive_infinity()> gamma[has_intercept];
  
  // from parameters_glm.stan
  vector[prior_dist == 7 ? sum(num_normals) : K] z_beta;
  real<lower=0> global[hs];
  vector<lower=0>[K] local[hs];
  real<lower=0> caux[hs > 0];
  vector<lower=0>[K] mix[prior_dist == 5 || prior_dist == 6];
  real<lower=0> one_over_lambda[prior_dist == 6];
  vector[q] z_b;
  vector[len_z_T] z_T;
  vector<lower=0,upper=1>[len_rho] rho;
  vector<lower=0>[len_concentration] zeta;
  vector<lower=0>[t] tau;
  
  real<lower=0> aux_unscaled[!is_binary];
}
transformed parameters {
  real aux[!is_binary];
 // tparameters_glm.stan
  vector[K] beta;
  vector[q] b; 
  vector[len_theta_L] theta_L;    
  
  if (!is_binary) {
    aux[1] = prior_dist_for_aux == 0 ? aux_unscaled[1] : (prior_dist_for_aux <= 2 ? 
             prior_scale_for_aux * aux_unscaled[1] + prior_mean_for_aux :
             prior_scale_for_aux * aux_unscaled[1]);
  }
  
  
  if      (prior_dist == 0) beta = z_beta;
  else if (prior_dist == 1) beta = z_beta .* prior_scale + prior_mean;
  else if (prior_dist == 2) for (k in 1:K) {
    beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
  }
  else if (prior_dist == 3) {
    real c2 = square(slab_scale) * caux[1];
    beta = hs_prior(z_beta, global, local, global_prior_scale, aux[1], c2);
  }
  else if (prior_dist == 4) {
    real c2 = square(slab_scale) * caux[1];
    beta = hsplus_prior(z_beta, global, local, global_prior_scale, aux[1], c2);
  }
  else if (prior_dist == 5) // laplace
    beta = prior_mean + prior_scale .* sqrt(2 * mix[1]) .* z_beta;
  else if (prior_dist == 6) // lasso
    beta = prior_mean + one_over_lambda[1] * prior_scale .* sqrt(2 * mix[1]) .* z_beta;
  else if (prior_dist == 7) { // product_normal
    int z_pos = 1;
    for (k in 1:K) {
      beta[k] = z_beta[z_pos];
      z_pos += 1;
      for (n in 2:num_normals[k]) {
        beta[k] *= z_beta[z_pos];
        z_pos += 1;
      }
      beta[k] *= prior_scale[k] ^ num_normals[k];
      beta[k] += prior_mean[k];
    }
  }


  if (!is_binary) {
    if (prior_dist_for_aux == 0) // none
      aux[1] = aux_unscaled[1];
    else {
      aux[1] = prior_scale_for_aux * aux_unscaled[1];
      if (prior_dist_for_aux <= 2) // normal or student_t
        aux[1] += prior_mean_for_aux;
    }
    
    theta_L = make_theta_L(len_theta_L, p, 
                           aux[1], tau, scale, zeta, rho, z_T);
  } else {
    theta_L = make_theta_L(len_theta_L, p, 
                           1.0, tau, scale, zeta, rho, z_T);
  }

  
  b = make_b(z_b, theta_L, p, l);
}
model {
  vector[N] eta;
  real dummy;
  
  real actual_aux = is_binary ? 1.0 : aux[1];
  
  eta = offset_;
  if (K > 0)
    eta += X * beta;
  eta += csr_matrix_times_vector2(N, q, w, v, u, b);
  
  if (has_intercept == 1)
    eta += gamma[1];

  if (has_weights == 0) {
    target += normal_lpdf(y | eta, actual_aux);
  } else {
    //vector[N] summands;
    //summands = pw_gauss(y, eta, actual_aux);
    //target += dot_product(weights, summands);
    // ignores sum of log of weights, target += 0.5 * sum(log(weights))
    target += -0.5 * N * log(6.283185307179586232 * actual_aux * actual_aux) - 0.5 * dot_product(weights, square((y - eta))) / (actual_aux * actual_aux);
  }

  if (!is_binary && prior_dist_for_aux > 0 && prior_scale_for_aux > 0) {
    real log_half = -0.693147180559945286;
    if (prior_dist_for_aux == 1)
      target += normal_lpdf(aux_unscaled[1] | 0, 1) - log_half;
    else if (prior_dist_for_aux == 2)
      target += student_t_lpdf(aux_unscaled[1] | prior_df_for_aux, 0, 1) - log_half;
    else 
     target += exponential_lpdf(aux_unscaled[1] | 1);
  }
  
  // priors_glm.stan
  // Log-priors for coefficients
       if (prior_dist == 1) target += normal_lpdf(z_beta | 0, 1);
  else if (prior_dist == 2) target += normal_lpdf(z_beta | 0, 1); // Student t via Cornish-Fisher expansion
  else if (prior_dist == 3) { // hs
    real log_half = -0.693147180559945286;
    target += normal_lpdf(z_beta | 0, 1);
    target += normal_lpdf(local[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
    target += normal_lpdf(global[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
    target += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);
  }
  else if (prior_dist == 4) { // hs+
    real log_half = -0.693147180559945286;
    target += normal_lpdf(z_beta | 0, 1);
    target += normal_lpdf(local[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
    target += normal_lpdf(local[3] | 0, 1) - log_half;
    // unorthodox useage of prior_scale as another df hyperparameter
    target += inv_gamma_lpdf(local[4] | 0.5 * prior_scale, 0.5 * prior_scale);
    target += normal_lpdf(global[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
    target += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);
  }
  else if (prior_dist == 5) { // laplace
    target += normal_lpdf(z_beta | 0, 1);
    target += exponential_lpdf(mix[1] | 1);
  }
  else if (prior_dist == 6) { // lasso
    target += normal_lpdf(z_beta | 0, 1);
    target += exponential_lpdf(mix[1] | 1);
    target += chi_square_lpdf(one_over_lambda[1] | prior_df[1]);
  }
  else if (prior_dist == 7) { // product_normal
    target += normal_lpdf(z_beta | 0, 1);
  }
  /* else prior_dist is 0 and nothing is added */
  
  // Log-prior for intercept  
  if (has_intercept == 1) {
    if (prior_dist_for_intercept == 1)  // normal
      target += normal_lpdf(gamma | prior_mean_for_intercept, prior_scale_for_intercept);
    else if (prior_dist_for_intercept == 2)  // student_t
      target += student_t_lpdf(gamma | prior_df_for_intercept, prior_mean_for_intercept, 
                               prior_scale_for_intercept);
    /* else prior_dist is 0 and nothing is added */
  }
  
  dummy = decov_lp(z_b, z_T, rho, zeta, tau, 
                   regularization, delta, shape, t, p);
}
