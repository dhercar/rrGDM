# Calculate Stress for a fitted rrGDM model

Computes a goodness-of-fit "stress" statistic for a fitted `rrGDM`
object, comparing the observed dissimilarities used to fit the model
against a model-implied quantity computed from the conditional latent
coordinates (i.e., including sample-level errors): either the raw
Euclidean distance between fitted latent positions (`scale = "link"`),
or the inverse-link-transformed expected dissimilarity
(`scale = "response"`).

For `type = "nonmetric"`, this is Kruskal's stress-1, \$\$Stress =
\sqrt{\frac{\sum (\hat{d} - D)^2}{\sum D^2}}\$\$ where \\D\\ is the
chosen model-implied quantity and \\\hat{d}\\ ("disparities") is the
isotonic (monotone) regression of \\D\\ on the observed dissimilarities.

For `type = "metric"`, observed and model-implied dissimilarities are
compared directly, with no regression step: \$\$Stress =
\sqrt{\frac{\sum (y\_{ij} - \hat{y}\_{ij})^2}{\sum y\_{ij}^2}}\$\$ where
\\y\_{ij}\\ are the observed dissimilarities and \\\hat{y}\_{ij}\\ the
model-implied ones. Because it does not allow for a linear
recalibration, it also penalises a poorly calibrated link. This is not
directly comparable to the stress reported by NMDS software.

Stress is returned as a proportion between 0 (best) and 1 (worst) for
bounded dissimilarities.

## Usage

``` r
stress(object, type = c("nonmetric", "metric"), scale = c("link", "response"))
```

## Arguments

- object:

  A fitted object of class `rrGDM`.

- type:

  Character string: `"nonmetric"` (default) computes Kruskal's stress-1
  using isotonic regression, as in standard NMDS. `"metric"` compares
  observed and model-implied dissimilarities directly (see Description).

- scale:

  Character string: `"link"` (default) uses the raw Euclidean distance
  between fitted latent positions in the ordination. This is generally
  the only scale on which stress values from RRGDM and other methods
  (e.g.,
  [`vegan::metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.html),
  [`MASS::isoMDS`](https://rdrr.io/pkg/MASS/man/isoMDS.html)) should be
  compared directly. `"response"` instead uses the inverse-link
  transformed expected dissimilarity (`object$inv_link(dist)`), which is
  on the same scale as the observed dissimilarities used to fit the
  model. This is the appropriate scale for `type = "metric"`.

## Value

A single numeric value: the stress.

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

stress(fit)
#> [1] 0.09938138
stress(fit, type = "metric", scale = "response")
#> [1] 0.06024978
```
