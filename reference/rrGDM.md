# Fit a reduced-rank generalised dissimilarity model

Fit a reduced-rank Generalised Dissimilarity Model using TMB

## Usage

``` r
rrGDM(
  Y,
  X = NULL,
  k = 2,
  binary = FALSE,
  formula = ~.,
  formula_sd = ~1,
  dist_method = "bray",
  scale_diss = NULL,
  family = "gaussian",
  link = NULL,
  se = TRUE,
  sample_error = TRUE,
  constraint = 1,
  trans = function(x) x,
  sparse = NULL,
  seed = NULL,
  init_jitter = 0.1,
  n_starts = 1,
  previous_best = TRUE,
  init_method = "dbRDA",
  start_fit = NULL,
  trace = FALSE,
  control = list()
)
```

## Arguments

- Y:

  A matrix or data frame containing the response variable (usually
  species abundance matrix).

- X:

  A data frame containing the predictor variables.

- k:

  An integer specifying the number of latent variables to be used in the
  model. Default is 2.

- binary:

  Logical indicating whether the response variable should be treated as
  binary (e.g., presence-absence). Default is FALSE.

- formula:

  Formula specifying the model to be fitted.

- formula_sd:

  Formula specifying the model for the dispersion parameter. Default is
  ~1 (intercept only).

- dist_method:

  A character string specifying the distance method to be used (e.g.,
  "euclidean", "bray"). Default is "bray".

- scale_diss:

  A numeric vector of length two indicating a new scale for
  dissimilarities. Useful for distributions bounded between 0 and 1
  (e.g., beta family). Ignored, with a warning, for the `"binomial"` and
  `"beta.binomial"` families.

- family:

  A character string specifying the family of the response variable:
  `"gaussian"`, `"binomial"`, `"beta"` or `"beta.binomial"`. Default is
  `"gaussian"`.

- link:

  Character string specifying the link function to be used:
  `"identity"`, `"negexp"` or `"neggaus"`. Default is 'identity' for
  Gaussian distribution, and 'negexp' for beta and binomial
  distributions.

- se:

  Logical indicating whether to calculate standard errors (TRUE) or not
  (FALSE). Default is TRUE.

- sample_error:

  Logical. If TRUE, includes site-level random effects to account for
  variation not explained by predictors. Default is TRUE.

- constraint:

  Integer. Specifies identifiability constraints. 0: no constraint; 1:
  Constraint applied to environmental loadings matrix (\\\Lambda\\); 2:
  Constraint applied to site random effects matrix (U). Default is 1.

- trans:

  Function used to transform the response matrix before computing
  dissimilarities (e.g., log, sqrt). Default is no transformation.

- sparse:

  Optional object of class `sparse_diss` (see
  [`sparse_diss`](https://dhercar.github.io/rrGDM/reference/sparse_diss.md)).
  If supplied, the model is fitted using only the pairs of samples
  listed in `sparse$edges`, instead of all pairs. It must be created
  from the same samples as `Y`. Default is NULL (all pairs).

- seed:

  An integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) at the start of the
  fit, so it changes the state of R's random number generator. Default
  is NULL (no seed is set).

- init_jitter:

  Numeric. Standard deviation of the Gaussian noise added to the
  starting values of every start after the first. Only used when
  `n_starts > 1`. Default is 0.1.

- n_starts:

  Integer. Number of optimisation starts. The first start uses the
  values from `init_method` (or `start_fit`); each further start adds
  Gaussian noise (SD `init_jitter`) to the best solution so far
  (`previous_best = TRUE`) or to the initial values (`FALSE`). With more
  than one start, only fits that converged and have a positive-definite
  Hessian are eligible as the final fit. Default is 1.

- previous_best:

  Logical. Controls where the jitter for multiple starts after the first
  is applied. If `TRUE` (default), each start jitters around the best
  solution found so far (perturb-and-reoptimise). If `FALSE`, every
  start jitters around the initial values (from `init_method`), so
  starts do not depend on earlier fits. Only relevant when
  `n_starts > 1`.

- init_method:

  A character string specifying a method used to obtain initial values
  for optimisation. Currently, only `'PCOA'` and `'dbRDA'` are
  available; `NULL` starts from small random values. Default is
  `'dbRDA'`. Ignored if `start_fit` is supplied.

- start_fit:

  Optional fitted `rrGDM` model used to initialise the optimisation
  (e.g., a fast sparse fit before refitting with more pairs). It must
  have the same samples, `k`, `family`, `link`, `constraint`, predictors
  and random-effect structure as the new model; the set of pairs
  (`sparse`) may differ. If supplied, `init_method` is ignored. Default
  is NULL.

- trace:

  Logical. If `TRUE`, TMB prints the optimisation progress and the
  objective of each start is reported. Default is `FALSE`.

- control:

  A list of control parameters passed to
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html). `iter.max` defaults
  to 1000.

## Value

An object of class `"rrGDM"`: a list with the components below.

**Estimation**

- `opt`:

  Result of [`nlminb`](https://rdrr.io/r/stats/nlminb.html) for the
  fixed-effect parameters

- `obj`:

  The TMB object created by
  [`MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html)

- `best_par`:

  Numeric vector with estimated parameters at the optimum, including the
  modes of random effects

- `sd_rep`:

  An [`sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html) object with
  standard errors and the covariance matrix of the estimates, or `NULL`
  if `se = FALSE` or the standard errors could not be computed. See
  [`include_se`](https://dhercar.github.io/rrGDM/reference/include_se.md)

- `pd_hess`:

  `TRUE` if the Hessian at the optimum is positive definite and not
  nearly singular, `FALSE` if not, and `NA` if there is no `sd_rep`

**Data**

- `Y`, `X`:

  The community matrix and the predictor data as supplied

- `Y_diss`, `Y_den`:

  Observed values for the fitted pairs. For the `"gaussian"` and
  `"beta"` families, `Y_diss` is the dissimilarity (rescaled to
  `scale_diss` if given) and `Y_den` is 1. For `"binomial"` and
  `"beta.binomial"`, `Y_diss` is the numerator of the dissimilarity and
  `Y_den` the denominator, so that `Y_diss / Y_den` is the
  dissimilarity.

- `pairs`:

  Two-column matrix with the row numbers of `Y` of the two samples in
  each fitted pair, in the same order as `Y_diss`.

- `sparse`:

  The
  [`sparse_diss`](https://dhercar.github.io/rrGDM/reference/sparse_diss.md)
  object supplied in `sparse`, or `NULL` if all pairs were used.

- `obs_full`:

  Only for fits with `sparse`: observed dissimilarities of *all* pairs
  of samples, in the order of
  [`dist`](https://rdrr.io/r/stats/dist.html)

**Design matrices**

- `X_design`:

  Design matrix of the predictors in `formula` (samples x predictors,
  without intercept; zero columns if there are no predictors).

- `W_design`:

  Design matrix of `formula_sd` (samples x terms, with intercept)

- `Z_design`:

  Sparse indicator matrix for the random intercepts in `formula`
  (samples x total number of groups).

- `map_re`:

  Integer vector giving, for each column of `Z_design`, the
  random-effect term it belongs to (index into `re_vars`). `integer(0)`
  if there are no random effects.

- `re_vars`:

  Names of the grouping variables with random intercepts. `character(0)`
  if there are none.

**Model specification**

- `k`, `family`, `link`:

  Number of latent dimensions, family and link used

- `inv_link`:

  Function that maps distances in the latent space to expected
  dissimilarities (identity, \\1 - e^{-d}\\ or \\1 - e^{-d^2}\\).

- `constraint`:

  Identifiability constraint used

- `formula`, `formula_sd`:

  The formulas supplied.

- `formula_fix`:

  `formula` without the random-effect terms (`~0` if there are no
  predictors).

- `terms_fix`, `terms_sd`:

  [`terms`](https://rdrr.io/r/stats/terms.html) objects of the
  fixed-effect and `formula_sd` models.

- `xlevels_fix`, `xlevels_sd`, `contrasts_fix`, `contrasts_sd`:

  Factor levels and contrasts used to build the design matrices. Used to
  build the design matrices for new data.

- `call`:

  The matched call.

## Details

With `constraint = 1` (or `2`) the upper triangle of the loadings
\\\Lambda\\ (or of the sample-level errors \\U\\) is fixed at 0 so that
the model is identifiable. A positive diagonal is currently not imposed,
so the sign of a latent axis is arbitrary.

## See also

[`summary.rrGDM`](https://dhercar.github.io/rrGDM/reference/summary.rrGDM.md),
[`anova.rrGDM`](https://dhercar.github.io/rrGDM/reference/anova.rrGDM.md),
[`wald_test`](https://dhercar.github.io/rrGDM/reference/wald_test.md),
[`get_LV`](https://dhercar.github.io/rrGDM/reference/get_LV.md),
[`predict.rrGDM`](https://dhercar.github.io/rrGDM/reference/predict.rrGDM.md),
[`biplot.rrGDM`](https://dhercar.github.io/rrGDM/reference/biplot.rrGDM.md),
[`stress`](https://dhercar.github.io/rrGDM/reference/stress.md),
[`cor_r2`](https://dhercar.github.io/rrGDM/reference/cor_r2.md),
[`include_se`](https://dhercar.github.io/rrGDM/reference/include_se.md),
[`sparse_diss`](https://dhercar.github.io/rrGDM/reference/sparse_diss.md)

## Examples

``` r
  # Toy community with 25 sites, 25 species, 2 environmental gradients
  set.seed(2)
  n <- 25 # sites
  p <- 25 # species

  env <- data.frame(temp = rnorm(n), moist = rnorm(n))
  eta <- as.matrix(env) %*% matrix(rnorm(2 * p), 2, p)
  Y <- matrix(rpois(n * p, exp(1 + 0.8 * eta)), n, p)

  fit <- rrGDM(Y, env, k = 2,
               formula = ~ temp + moist,
               family = "binomial")
#> ℹ Preparing data
#> ✔ Preparing data ... done
#> 
#> ℹ Finding initial parameters (dbRDA)
#> ✔ Finding initial parameters (dbRDA) ... done
#> 
#> ℹ optimisation
#> ℹ Start 1/1: objective = 1670.1661
#> ℹ optimisation
#> ✔ optimisation ... done
#> 
  summary(fit)
#> Call:
#> rrGDM(Y = Y, X = env, k = 2, formula = ~temp + moist, family = "binomial")
#> 
#> Family: binomial (link: negexp); 2 latent dimensions; 25 samples (300 pairs)
#> Log-likelihood: -1670.17 (df = 4); AIC: 3348.33
#> 
#> ── Effects on composition (location) ──
#> 
#>   term   Size    Wald Df Pr_Chisq Sig
#>   temp 0.5957 603.881  1   <1e-06 ***
#>  moist 0.6527 644.000  2   <1e-06 ***
#> ── Effects on dispersion (log of sample-level SD) ──
#> 
#>         term Estimate     SE      z   Pr_z Sig
#>  (Intercept)   -2.030 0.1088 -18.66 <1e-06 ***
#> 
#> Sample-level SD at baseline: tau = exp(Intercept) = 0.1313
#> ── Fit of the latent space (conditional) ──
#> 
#>   Stress (metric, response scale):    0.06025
#>   Stress (non-metric, link scale):    0.09938
#>   corR2:                              0.9550
#>   Stress and corR2 show how well the latent space captures dissimilarities,
#>   not the explanatory power of the predictors.
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

  # Model the sample-level dispersion too
  fit2 <- rrGDM(Y, env, k = 2,
                formula = ~ temp + moist,
                formula_sd = ~ temp,
                family = "binomial")
#> ℹ Preparing data
#> ✔ Preparing data ... done
#> 
#> ℹ Finding initial parameters (dbRDA)
#> ✔ Finding initial parameters (dbRDA) ... done
#> 
#> ℹ optimisation
#> ℹ Start 1/1: objective = 1669.7663
#> ℹ optimisation
#> ✔ optimisation ... done
#> 
  summary(fit2)
#> Call:
#> rrGDM(Y = Y, X = env, k = 2, formula = ~temp + moist, formula_sd = ~temp, 
#>     family = "binomial")
#> 
#> Family: binomial (link: negexp); 2 latent dimensions; 25 samples (300 pairs)
#> Log-likelihood: -1669.77 (df = 5); AIC: 3349.53
#> 
#> ── Effects on composition (location) ──
#> 
#>   term   Size    Wald Df Pr_Chisq Sig
#>   temp 0.5879 564.089  1   <1e-06 ***
#>  moist 0.6442 589.395  2   <1e-06 ***
#> ── Effects on dispersion (log of sample-level SD) ──
#> 
#>         term Estimate     SE      z   Pr_z Sig
#>  (Intercept)   -2.086 0.1224 -17.04 <1e-06 ***
#>         temp   0.1297 0.1440 0.9005  0.368    
#> 
#> Sample-level SD at baseline: tau = exp(Intercept) = 0.1241
#> ── Fit of the latent space (conditional) ──
#> 
#>   Stress (metric, response scale):    0.06027
#>   Stress (non-metric, link scale):    0.09932
#>   corR2:                              0.9550
#>   Stress and corR2 show how well the latent space captures dissimilarities,
#>   not the explanatory power of the predictors.
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
