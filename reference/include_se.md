# Calculate standard errors

Calculate standard errors

## Usage

``` r
include_se(object, ...)
```

## Arguments

- object:

  A fitted `rrGDM` object.

- ...:

  Further arguments passed to
  [`sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html).

## Value

A model object, including the sd report

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

fit0 <- rrGDM(Y, env, k = 2, formula = ~ temp + moist, family = "binomial", se = FALSE)
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
fit0 <- include_se(fit0)
wald_test(fit0)
#> 
#> ── Wald tests: environmental loadings ──
#> 
#>   term   Size    Wald Df Pr_Chisq Sig
#>   temp 0.5957 603.857  1   <1e-06 ***
#>  moist 0.6526 643.979  2   <1e-06 ***
#> ── Wald tests: dispersion effects ──
#> 
#> ℹ No dispersion terms.
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Wald tests on loadings can be anti-conservative with few samples.
```
