# Print method for sparse_diss objects

Prints an overview of a `sparse_diss` object: the number of samples and
pairs (full and included), the selection settings, the pairs kept by
type and the number of pairs per sample.

## Usage

``` r
# S3 method for class 'sparse_diss'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `sparse_diss`, as returned by
  [`sparse_diss`](https://dhercar.github.io/rrGDM/reference/sparse_diss.md).

- digits:

  Number of significant digits used to print dissimilarities. Default is
  3.

- ...:

  Currently unused.

## Value

Invisibly returns `x`.
