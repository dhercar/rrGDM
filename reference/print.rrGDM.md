# Print method for rrGDM models

Prints an overview of a fitted `rrGDM` object.

## Usage

``` r
# S3 method for class 'rrGDM'
print(x, digits = 3, n_u = 5, ...)
```

## Arguments

- x:

  A fitted model of class `rrGDM`.

- digits:

  Number of significant digits used to print loadings. Default is 3.

- n_u:

  Number of rows to print for the sample-level error matrix `u` and, for
  each grouping variable, for the random-effect matrix (these have one
  row per sample or per group, so they are truncated). Default is 5.

- ...:

  Currently unused.

## Value

Invisibly returns `x`.

## Details

Convergence information includes the convergence code returned by
[`nlminb`](https://rdrr.io/r/stats/nlminb.html) (0 indicates successful
convergence), whether the Hessian at the optimum is positive definite
(only available if standard errors were requested) and the largest
absolute gradient of the fixed parameters. Loadings are printed as
estimated (i.e., unrotated). Use
[`get_LV`](https://dhercar.github.io/rrGDM/reference/get_LV.md) to
extract the full matrices.
