#define TMB_LIB_INIT R_init_rrGDM
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{

  // DATA
  DATA_VECTOR(Y);      // pairwise dissimilarities (num)
  DATA_VECTOR(Y_den);  // pairwise dissimilarities (den)
  DATA_MATRIX(X);   // environmental predictors (n_samples x n_env)
  DATA_MATRIX(W);   // predictors on sigma
  DATA_IMATRIX(D);   // mapping of sample pairs
  DATA_SPARSE_MATRIX(Z);  // other random effects
  DATA_INTEGER(has_sample_re);  // flag for re (sample-level)
  DATA_INTEGER(has_random);   // flag for re (other)
  DATA_INTEGER(family);    // family (0 = normal, 1 = binomial, 2 = beta), 3 = betabinomial)
  DATA_INTEGER(link);    // link (0 = identity, 1 = neg_exp, 2 = neg_gaus)
  DATA_INTEGER(constraint);  // flag for constraint (0 = none, 1 = lambda, 2 = site re)
  DATA_IVECTOR(map_re);  // mapping of RE columns to vars
  DATA_INTEGER(n_lv);    // n latent dimensions (k)

  // Parameters
  PARAMETER_MATRIX(lambda);  // env. loadings
  PARAMETER_MATRIX(u_sample);   // U (sample-level error)
  PARAMETER(log_scale);     // Scale parameter (dispersion)
  PARAMETER_VECTOR(beta);   // effects on sigma (dispersion)
  PARAMETER_MATRIX(u);  // additional random effects
  PARAMETER_VECTOR(log_sigma_re);  // vector of SD for each RE term

  // Derived variables
  int n_pairs = Y.size();   // number of pairs (dissimilarities)
  int n_samples = X.rows();  // number of samples (dissimilarities)
  int n_env = X.cols();     // number of env. variables

  vector<Type> sigma_re = exp(log_sigma_re); // sigma of sample error
  Type SD = 0;     // SD of Gaussian family (pair error)
  Type phi = 0;    // phi of beta (and beta binomial) distribution (pair error)


  // ---------- 3. LATENT COMPOSITIONS ---------
  matrix<Type> L_m = X * lambda;
  matrix<Type> L = L_m;

  if (has_sample_re) { // add sample level error
    L += u_sample;
  }

  if (has_random && u.size() > 0) {   // add other random effects
    L += Z * u;
  }

  //  SD model
  vector<Type> l_sample_sd = W * beta.matrix();

  // U likelihood contr. to likelihood
  Type nll = 0.0;
  if (has_sample_re) {
    for (int i = 0; i < n_samples; i++) {
      Type sample_sd_i = exp(l_sample_sd(i));
      for (int j = 0; j < n_lv; j++) {
        if (constraint == 2 && j > i) continue;  // Skip upper triangle if constraint == 2
        nll -= dnorm(u_sample(i, j), Type(0), sample_sd_i, true);
  }}}

  // Other random effects: contribution to NLL
  if (has_random) {
    for (int i = 0; i < u.rows(); i++) {
      for (int k = 0; k < n_lv; k++) {
        nll -= dnorm(u(i, k), Type(0), sigma_re(map_re[i]), true);
  }}}


  // Distribution-specific scale parameters
  if (family == 0) {
    SD = exp(log_scale);
  } else if (family == 2 || family == 3) {
    phi = exp(log_scale);
  }

  //----------PAIRWISE LOOP (START)----------
  // !!! Could use matrix multiplication instead, left as loop to allow parallel comp.
  vector<Type> mu_vec(n_pairs);
  for (int i = 0; i < n_pairs; i++) {

    // Select two samples
    int s1 = D(i, 0);
    int s2 = D(i, 1);

    // Calculate Euclidean distances
    Type dist_2 = 0.0;
    for(int k = 0; k < n_lv; k++) {
      Type diff = L(s1, k) - L(s2, k);
      dist_2 += diff * diff;
    }

    dist_2 += Type(1e-10);   // avoid exact 0
    Type dist = sqrt(dist_2);

    // Inverse link
    Type mu;
    switch (link) {
    case 1:
      mu = Type(1) - exp(-(dist));
      break;
    case 2:
      mu = Type(1) - exp(-(dist_2));
      break;
    default:
      mu = dist;
    break;
    }

    if (family != 0) mu = Type(1e-9) + (Type(1) - Type(2e-9)) * mu;
    mu_vec(i) = mu;
  }

  //---------- END PAIRWISE LOOP ----------

  // likelihood
  switch (family) {
    case 0: // Gaussian
      nll -= sum(dnorm(Y, mu_vec, SD, true));
      break;
    case 1: // Binomial
      nll -= sum(dbinom(Y, Y_den, mu_vec, true));
      break;
    case 2: { // Beta
        vector<Type> a = mu_vec * phi;
        vector<Type> b = (Type(1.0) - mu_vec) * phi;
        nll -= sum(dbeta(Y, a, b, true));
        break;
      }
    case 3: { // Beta-Binomial
      for(int i = 0; i < n_pairs; i++) {
      Type a = mu_vec(i) * phi;
      Type b = (Type(1.0) - mu_vec(i)) * phi;
      nll -= lgamma(Y_den(i) + Type(1)) - lgamma(Y(i) + Type(1)) - lgamma(Y_den(i) - Y(i) + Type(1))
        + lgamma(Y(i) + a) + lgamma(Y_den(i) - Y(i) + b) - lgamma(Y_den(i) + a + b)
        - lgamma(a) - lgamma(b) + lgamma(a + b);
    }
      break;
  }}

  // Derived quantities:
  vector<Type> lambda_mod(n_env); // total env. effects
  for(int i = 0; i < n_env; i++){
    lambda_mod(i) = sqrt(lambda.row(i).array().square().sum());
  }

  REPORT(lambda);
  REPORT(lambda_mod);
  REPORT(beta);
  REPORT(sigma_re);
  REPORT(L);
  REPORT(L_m);

  if (has_sample_re) REPORT(u_sample);
  if (has_random) REPORT(u);

  ADREPORT(lambda);
  ADREPORT(beta);
  ADREPORT(sigma_re);
  if (family == 0) ADREPORT(SD);
  if (family == 2 || family == 3) ADREPORT(phi);

  return nll;
}
