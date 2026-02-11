#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // ---------- 1. DATA AND PARAMETERS ---------

  // Data
  DATA_VECTOR(Y);               // pairwise dissimilarities ()[0,1] for Beta distribution)
  DATA_VECTOR(Y_den);           // vector of total amount of species
  DATA_MATRIX(X);               // environmental var (n_samples x n_env)
  DATA_MATRIX(W);               // effects on sigma
  DATA_MATRIX(D);               // matrix mapping Z to sample pairs (n_pairs x 2)
  DATA_SPARSE_MATRIX(Z);        // extra random effects
  DATA_INTEGER(has_sample_re);  // flag for re (sample-level)
  DATA_INTEGER(has_random);     // flag for re (other)
  DATA_INTEGER(family);         // family (0 = normal, 1 = binomial, 2 = beta)
  DATA_INTEGER(link);           // link function (0 = identity, 1 = neg_exp, 2 = neg_gaus)
  DATA_INTEGER(constraint);     // flag for constraint (0 = no constraint, 1 = constraint on beta, 2 = constraint on site re)
  DATA_IVECTOR(map_re);         // mapping from RE columns
  DATA_INTEGER(n_lv);           // n latent dimensions

  // Parameters
  PARAMETER_VECTOR(lambda_diag);          // lambda diagonal(env. loadings)
  PARAMETER_VECTOR(lambda_upper);         // lambda, upper triangle
  PARAMETER_VECTOR(lambda_lower);         // lambda, lower triangle
  PARAMETER_MATRIX(u_sample);           // U (sample error)
  PARAMETER(log_scale);                 // Scale parameter
  PARAMETER_VECTOR(beta);             // effects on sigma
  PARAMETER_MATRIX(u);                  // additional random effects ()intercepts)
  PARAMETER_VECTOR(log_sigma_re);       // log of SD for each RE term

  // Derived variables
  int n_pairs = Y.size();               // number of pairs (dissimilarities)
  int n_samples = X.rows();             // number of samples (dissimilarities)
  int n_env = X.cols();                 // number of env. variables
  vector<Type> sigma_re = exp(log_sigma_re);

  Type SD = 0;                          // SD of Gaussian family
  Type phi = 0;                         // phi of beta (and beta binomial) distribution


  // ---------- 2. LOADINGS AND CONSTRAINS ---------

  // CONSTRUCT LAMBDA MATRIX (environmental loadings)
  matrix<Type> lambda(n_env, n_lv);
  lambda.setZero();

  // Diagonal
  int dim_lambda = (n_env < n_lv) ? n_env : n_lv;
  for(int i = 0; i < dim_lambda; i++){
    if(constraint == 1){
      lambda(i, i) = exp(lambda_diag(i));
    } else {
      lambda(i, i) = lambda_diag(i);
    }
  }

  // Lower
  int idx = 0;
  for(int i = 0; i < n_env; i++){
    for(int j = 0; j < n_lv; j++){
      if(i > j){
        lambda(i, j) = lambda_lower(idx++);
      }
    }
  }

  // Upper
  if (constraint != 1) { // Upper tri remains 0s and ignored through 'map' in model fitting function when constraint = 2
    idx = 0;
    for(int i = 0; i < n_env; i++){
      for(int j = 0; j < n_lv; j++){
        if(i < j){
          lambda(i, j) = lambda_upper(idx++);
        }
      }
    }
  }


  // ---------- 3. LATENT (COM)POSITIONS ---------
  matrix<Type> L_m = X * lambda;         // expected position
  matrix<Type> L = L_m + u_sample;     // sample-level error
  if (has_random && u.size() > 0) {    // additional random effects
    L += Z * u;
  }

  // SD sample-level error (distributional model)
  vector<Type> l_sample_sd = W * beta.matrix();

  // initiate negative log-likelihood
  Type nll = 0.0;

  // Sample-level error contribution to nll
  if (has_sample_re) {
    for (int i = 0; i < n_samples; i++) {
      Type sample_sd_i = exp(l_sample_sd(i));
      for (int j = 0; j < n_lv; j++) {
        // Skip upper triangle when constraint == 2
        if (constraint == 2 && j > i) continue;
        nll -= dnorm(u_sample(i, j), Type(0), sample_sd_i, true);
      }
    }
  }

  // Other random effect contribution to nll
  if (has_random) {
    for (int i = 0; i < u.rows(); i++) {
      for (int k = 0; k < n_lv; k++) {
        nll -= dnorm(u(i, k), Type(0), sigma_re(map_re[i]), true);
      }
    }
  }

  // ---------- 4. COMPUTE PAIRWISE DISTANCES AND LIKELIHOOD OF DISSIMILARITIES ---------

  // Distribution-specific scale parameters (if any)
  if (family == 0) {
    SD = exp(log_scale); // Gaussian SD
  } else if (family == 2) {
    phi = exp(log_scale); // Beta scale
  }

  // Loop through each pair of sites
  for (int i = 0; i < n_pairs; i++) {
    // Convert matrix elements to integers
    int s1 = CppAD::Integer(D(i, 0));
    int s2 = CppAD::Integer(D(i, 1));

    // Calculate squared Euclidean distance between latent positions (L)
    Type dist_2 = 0.0;
    for(int k = 0; k < n_lv; k++) {
      Type diff = L(s1, k) - L(s2, k);
      dist_2 += diff * diff;
    }

    // Apply inverse link
    Type mu;
    switch (link) {
    case 1:
      mu = Type(1) - exp(-sqrt(dist_2));
      break;
    case 2:
      mu = Type(1) - exp(-dist_2);
      break;
    default:
      mu = sqrt(dist_2);
    break;
    }

    // likelihood
    switch (family) {
    case 0: // Gaussian
      nll -= dnorm(Y(i), mu, SD, true);
      break;
    case 1: // Binomial
      nll -= dbinom(Y(i), Y_den(i), mu, true);
      break;
    case 2: { // Beta
        Type a = mu * phi;
        Type b = (Type(1) - mu) * phi;
        nll -= dbeta(Y(i), a, b, true);
        break;
      }
    }

  }

  // Report intermediate variables
  REPORT(lambda);                            // environmental loadings
  REPORT(beta);                              // Effects on dispersion (SD of sample-level error)
  REPORT(sigma_re);                          // SD of other random effects
  REPORT(L);                                 // Site loadings (conditional) for plot
  REPORT(L_m);                               // Site loadings (marginal) for plot
  if (has_sample_re) REPORT(u_sample);       // U: Sample-level error loadings

  // report SE
  ADREPORT(lambda);                        // B: Env. loadings
  ADREPORT(beta);                          // Effects on dispersion (SD sample error)
  ADREPORT(sigma_re);                      // SD of other random effects

  if (family == 0) ADREPORT(SD);
  if (family == 2) ADREPORT(phi);

  return nll;
}
