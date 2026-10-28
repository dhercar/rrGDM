# Conditional correlation-based R2 (corR2) for a fitted rrGDM model

Computes the conditional `corR2` of a fitted `rrGDM` object: the squared
Pearson correlation between the observed dissimilarities used to fit the
model and the dissimilarities implied by the conditional latent
coordinates (i.e., including sample-level errors and any other random
effects), after applying the inverse link function.

## Usage

``` r
cor_r2(object, message = FALSE)
```

## Arguments

- object:

  A fitted object of class `rrGDM`.

- message:

  Logical. Should a message clarifying the interpretation of the corR2
  be printed? Default is `TRUE`.

## Value

A single numeric value: the conditional corR2, between 0 and 1.

## Details

The corR2 is useful to assess how well the latent space captures the
observed dissimilarities (e.g., to choose the number of latent
dimensions `k`), but it is **not** a measure of the explanatory power of
the predictors. Because the conditional coordinates include the
sample-level errors, the corR2 can be high even if the predictors
explain little. Like
[`stress`](https://dhercar.github.io/rrGDM/reference/stress.md), it
increases as `k` increases, so it is best used descriptively.

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

cor_r2(fit)
#> [1] 0.9549986
```
