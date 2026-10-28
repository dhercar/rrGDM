# Wald Tests for Environmental Loadings and Dispersion Effects

Wald tests for sets of predictor effects in a fitted `rrGDM` model:

- **Environmental loadings** (`lambda`): for each predictor, a
  multivariate (multi-df) Wald test of the null that all free loadings
  of that predictor (one per latent dimension) are jointly zero. This
  accounts for the potential covariance between loadings within a
  predictor. The `Size` column gives the modulus (Euclidean norm) of the
  loading vector across the `k` latent dimensions, a rotation-invariant
  measure of the predictor's total effect on community composition.

- **Dispersion effects** (`beta`): for each predictor in the dispersion
  model, a univariate Wald z-test, since each coefficient is a scalar
  effect on \\\log(\tau)\\.

## Usage

``` r
wald_test(object, include_intercept = FALSE)
```

## Arguments

- object:

  A fitted model of class `rrGDM`, fitted with `se = TRUE`.

- include_intercept:

  Logical. Should the intercept term of `formula_sd` be included in the
  `beta` output? Default is `FALSE`.

## Value

A list of class `wald_rrGDM` with two elements:

- `lambda`: a `data.frame` with one row per predictor in `X_design`
  (zero rows if the mean model has no predictors), with columns `term`,
  `Size`, `Wald`, `Df` (see Details), `Pr_Chisq`, `Sig`.

- `beta`: a `data.frame` with one row per predictor in `W_design`
  (excluding the intercept unless `include_intercept = TRUE`), with
  columns `term`, `Estimate`, `SE`, `z`, `Pr_z`, `Sig`.

## Details

This test is preferred over testing each latent dimension separately
(e.g., with per-dimension p-values corrected for multiple comparisons),
because it explicitly accounts for the covariance between loadings. If
the covariance matrix \\\boldsymbol{\Sigma}\_p\\ of a predictor's
loadings is singular or ill-conditioned (judged on the correlation
scale, so independently of the units of the coefficients; reciprocal
condition number below `1e-8`), the Wald statistic and p-value for that
predictor are set to `NA` with a warning. This typically indicates
unidentified or nearly collinear coefficients, for which the test is not
meaningful. If the Hessian at the optimum was not positive definite, the
function stops with an error.

The degrees of freedom of the `lambda` test depend on
`object$constraint`:

- Under `constraint == 1`, the loadings matrix is constrained to be
  lower triangular: predictor (row) `p` only has free loadings in latent
  dimensions `1:min(p, k)`, the remainder being structural zeros that
  are fixed (not estimated) when the model is fit. The Wald test for
  predictor `p` therefore uses `Df = min(p, k)`, and only the free
  loadings (and their covariance sub-block) enter the test statistic.

- Under `constraint %in% c(0, 2)`, `lambda` is unconstrained, so every
  predictor is tested with `Df = k`. Note that under `constraint == 0`
  the model is not identified under rotation of the latent space, so
  individual loadings (and therefore this test) should be interpreted
  with caution.

## See also

[`anova.rrGDM`](https://dhercar.github.io/rrGDM/reference/anova.rrGDM.md)
for tests by variable,
[`summary.rrGDM`](https://dhercar.github.io/rrGDM/reference/summary.rrGDM.md)

## Examples

``` r
set.seed(2)
# simulate data
n <- 25 # sites
p <- 25 # species

env <- data.frame(temp = rnorm(n), moist = rnorm(n))
eta <- as.matrix(env) %*% matrix(rnorm(2 * p), 2, p)
Y <- matrix(rpois(n * p, exp(1 + 0.8 * eta)), n, p,
            dimnames = list(paste0("site", 1:n), paste0("sp", 1:p)))

# fit model
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

wald_test(fit)
#> 
#> ── Wald tests: environmental loadings ──
#> 
#>   term   Size    Wald Df Pr_Chisq Sig
#>   temp 0.5957 603.881  1   <1e-06 ***
#>  moist 0.6527 644.000  2   <1e-06 ***
#> ── Wald tests: dispersion effects ──
#> 
#> ℹ No dispersion terms.
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Wald tests on loadings can be anti-conservative with few samples.
```
