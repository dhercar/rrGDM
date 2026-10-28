# Wald analysis of variance table for rrGDM models

Joint Wald chi-square tests for each *variable* in a fitted `rrGDM`
model, however many terms (design-matrix columns) the variable
contributes. Where
[`wald_test`](https://dhercar.github.io/rrGDM/reference/wald_test.md)
gives one test per column of the design matrix,
[`anova()`](https://rdrr.io/r/stats/anova.html) groups every coefficient
that belongs to the same variable into a single test, so that:

- a factor with several levels gets one test rather than one per level

- terms such as `x + I(x^2)`, `poly(x, 2)` or `splines::ns(x, 3)` are
  tested together as "the effect of `x`"

- a variable that enters both the location model (`formula`) and the
  dispersion model (`formula_sd`) can be tested on both at once, or on
  each separately.

## Usage

``` r
# S3 method for class 'rrGDM'
anova(object, ..., interactions = TRUE, separate = FALSE)
```

## Arguments

- object:

  A fitted model of class `rrGDM`, fitted with `se = TRUE`.

- ...:

  Not used. Only a single model can be supplied (comparison of several
  models is not supported).

- interactions:

  Logical. Controls how terms involving more than one variable (e.g.,
  `x:f`) are handled. If `TRUE` (default), a variable's test includes
  its main-effect terms *and* every interaction term it takes part in,
  i.e., it tests whether the variable has any effect on composition (or
  dispersion) at all. If `FALSE`, interaction terms are omitted.
  Variables that only appear in interactions are then omitted.

- separate:

  Logical. A variable can appear in both the location model (`formula`)
  and the dispersion model (`formula_sd`). If `FALSE` (default), these
  are tested jointly in a single row
  (`effect = "location + dispersion"`), accounting for the covariance
  between location and dispersion coefficients. If `TRUE`, the location
  and dispersion effects of each variable get a row each.

## Value

A `data.frame` of class `anova_rrGDM`, with one row per variable (or per
variable and effect type if `separate = TRUE`) and columns:

- `variable`: name of the variable.

- `effect`: `"location"`, `"dispersion"` or `"location + dispersion"`.

- `terms`: the model terms grouped in the test.

- `Df`: degrees of freedom (number of free coefficients tested).

- `Wald`: the Wald chi-squared statistic.

- `Pr_Chisq`: p-value.

- `Sig`: significance code.

## Details

Intercepts and random effects are excluded.

For each row, the null hypothesis is that all of the coefficients
involved are jointly zero. As in
[`wald_test`](https://dhercar.github.io/rrGDM/reference/wald_test.md),
structural zeros of the loadings matrix under `constraint == 1` are
excluded from both the statistic and the degrees of freedom. If the
joint covariance matrix of the coefficients being tested is singular or
ill-conditioned (judged on the correlation scale; reciprocal condition
number below `1e-8`), the Wald statistic and p-value for that row are
set to `NA` with a warning, since the test is not meaningful when
coefficients are unidentified or nearly collinear. If the Hessian at the
optimum was not positive definite, the function stops with an error.

With `interactions = TRUE` the test compares the model against one where
the variable is dropped entirely. With `interactions = FALSE` the test
of a variable that interacts with others refers to its effect when the
interacting variables are at their reference level (factors) or zero
(continuous predictors).

As with
[`wald_test`](https://dhercar.github.io/rrGDM/reference/wald_test.md),
tests on location effects can be anti-conservative with few samples.

## See also

[`wald_test`](https://dhercar.github.io/rrGDM/reference/wald_test.md)
for per-term tests, effect sizes and dispersion estimates.

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

anova(fit)
#> 
#> ── Wald tests by variable ──
#> 
#> • Interaction terms included in each variable's test
#> • Location and dispersion effects of the same variable tested jointly
#> 
#> variable  effect    Df     Wald  Pr_Chisq
#> temp      location   1  603.881    <1e-06  ***
#> moist     location   2  644.000    <1e-06  ***
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Tests on location effects can be anti-conservative with few samples.
anova(fit, separate = TRUE)
#> ── Wald tests by variable ──
#> 
#> • Interaction terms included in each variable's test
#> • Location and dispersion effects tested separately
#> 
#> variable  effect    Df     Wald  Pr_Chisq
#> temp      location   1  603.881    <1e-06  ***
#> moist     location   2  644.000    <1e-06  ***
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Tests on location effects can be anti-conservative with few samples.
```
