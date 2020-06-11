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
  
  vector csr_matrix_times_vector2(int m, int n, vector w, 
                                  int[] v, int[] u, vector b);
  
  vector pw_gauss(vector y, vector eta, real sigma) {
    return -0.5 * log(6.283185307179586232 * sigma) - 0.5 * square((y - eta)) / sigma;
  }
}
data {
  // dimensions
  int<lower=0> N; // number of observations
    
  int<lower=0> len_y;
  real lb_y; // lower bound on y
  real<lower=lb_y> ub_y; // upper bound on y
  vector<lower=lb_y, upper=ub_y>[len_y] y;
  
  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
  int<lower=0,upper=3> prior_dist_for_aux;
  
    
  int<lower=0,upper=1> has_weights;  // 0 = No, 1 = Yes
  vector[has_weights ? N : 0] weights;
  
  vector[N] offset_;
  
  real<lower=0> prior_scale_for_aux;
  real<lower=0> prior_mean_for_aux;
  real<lower=0> prior_df_for_aux;
    
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
  vector[q] z_b;
  vector[len_z_T] z_T;
  vector<lower=0,upper=1>[len_rho] rho;
  vector<lower=0>[len_concentration] zeta;
  vector<lower=0>[t] tau;
  
  real<lower=0> aux_unscaled;
}
transformed parameters {
  real aux = prior_dist_for_aux == 0 ? aux_unscaled : (prior_dist_for_aux <= 2 ? 
             prior_scale_for_aux * aux_unscaled + prior_mean_for_aux :
             prior_scale_for_aux * aux_unscaled);

  vector[q] b;
  vector[len_theta_L] theta_L;    

  if (prior_dist_for_aux == 0) {
    aux = aux_unscaled;
  } else {
    aux = prior_scale_for_aux * aux_unscaled;
    if (prior_dist_for_aux <= 2) // normal or student_t
      aux += prior_mean_for_aux;
  }

  theta_L = make_theta_L(len_theta_L, p, 
                         aux, tau, scale, zeta, rho, z_T);
  b = make_b(z_b, theta_L, p, l);
}
model {
  vector[N] eta;
  real dummy;
  
  eta = offset_;
  eta += csr_matrix_times_vector2(N, q, w, v, u, b);

  if (has_weights == 0) {
    target += normal_lpdf(y | eta, aux);
  } else {
    vector[N] summands;
    summands = pw_gauss(y, eta, aux);
    target += dot_product(weights, summands);
  }

  if (prior_dist_for_aux > 0 && prior_scale_for_aux > 0) {
    real log_half = -0.693147180559945286;
    if (prior_dist_for_aux == 1)
      target += normal_lpdf(aux_unscaled | 0, 1) - log_half;
    else if (prior_dist_for_aux == 2)
      target += student_t_lpdf(aux_unscaled | prior_df_for_aux, 0, 1) - log_half;
    else 
     target += exponential_lpdf(aux_unscaled | 1);
  }
  
  dummy = decov_lp(z_b, z_T, rho, zeta, tau, 
                   regularization, delta, shape, t, p);
}
