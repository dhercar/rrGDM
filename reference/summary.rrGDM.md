# Summary method for rrGDM models

Summarises a fitted `rrGDM` object. The summary reports:

- the model call and the (approximate marginal) log-likelihood;

- for each predictor in the mean model, the effect size (`Size`, the
  modulus of its loadings across the latent dimensions) and a
  multivariate Wald test (see
  [`wald_test`](https://dhercar.github.io/rrGDM/reference/wald_test.md));

- the intercept and slopes of the dispersion model (effects on
  \\\log(\tau)\\, where \\\tau\\ is the sample-level SD), with
  univariate Wald z-tests;

- the standard deviation of additional random effects, if any;

- the observation-level (pairwise) parameter of the response
  distribution, if the family has one (`phi` for beta and beta-binomial,
  `SD` for Gaussian)

- the conditional stress and corR2 of the model.

## Usage

``` r
# S3 method for class 'rrGDM'
summary(object, ...)
```

## Arguments

- object:

  A fitted model of class `rrGDM`.

- ...:

  Currently unused.

## Value

A list of class `summary.rrGDM`, with elements including `call`,
`loglik`, `df`, `aic`, `lambda` (see
[`wald_test`](https://dhercar.github.io/rrGDM/reference/wald_test.md)),
`beta`, `tau0` (sample-level SD when all dispersion predictors are at
their reference value), `sigma_re`, `scale_par`, `stress_metric`,
`stress_nonmetric` and `cor_r2`.

## Details

The log-likelihood is the approximate marginal log-likelihood obtained
with the Laplace approximation, and `df` is the number of (non-mapped)
fixed parameters estimated.

Two stress values are reported, both computed from the conditional
latent coordinates (see
[`stress`](https://dhercar.github.io/rrGDM/reference/stress.md)): the
metric stress on the dissimilarity (response) scale, and Kruskal's
non-metric stress on the latent distances (link scale), which is the one
comparable with NMDS software. The corR2 (see
[`cor_r2`](https://dhercar.github.io/rrGDM/reference/cor_r2.md)) is the
conditional squared correlation between observed and fitted
dissimilarities. Both indicate how well the latent space captures the
dissimilarities, but not the explanatory power of the predictors.

If the model was fitted without standard errors (`se = FALSE`), the
Hessian at the optimum is not positive definite, or the covariance
matrix of the estimates is missing, a warning is issued and only point
estimates are reported.

If the stress or corR2 cannot be computed, a warning is given and the
value is `NA`.

## See also

[`wald_test`](https://dhercar.github.io/rrGDM/reference/wald_test.md),
[`stress`](https://dhercar.github.io/rrGDM/reference/stress.md),
[`cor_r2`](https://dhercar.github.io/rrGDM/reference/cor_r2.md)
