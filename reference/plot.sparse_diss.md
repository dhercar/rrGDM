# Plot method for sparse_diss objects

[`plot`](https://rdrr.io/r/graphics/plot.default.html) method for
objects of class `sparse_diss`. Shows the samples on two axes of a PCoA
of the full dissimilarity matrix, with the selected pairs drawn as
segments between samples.

## Usage

``` r
# S3 method for class 'sparse_diss'
plot(
  x,
  y = NULL,
  k = c(1, 2),
  edge_types = NULL,
  edge_col_by = "type",
  edge_alpha = 0.4,
  edge_linewidth = 0.3,
  add_points = TRUE,
  point_size = 1.5,
  ...
)
```

## Arguments

- x:

  An object of class `sparse_diss`, as returned by
  [`sparse_diss`](https://dhercar.github.io/rrGDM/reference/sparse_diss.md).

- y:

  Currently unused.

- k:

  Integer vector of length 2 indicating which PCoA axes to plot. Must be
  within the `k_pcoa` axes stored in `x`. Default is `c(1, 2)`.

- edge_types:

  Optional character vector with the type of pairs to draw (any of
  `'knn'`, `'random'`, `'zero'`). Default is `NULL` (all).

- edge_col_by:

  Character string: colour segments by pair `'type'` (default) or by
  dissimilarity (`'diss'`).

- edge_alpha:

  Numeric transparency of segments. Default is 0.4.

- edge_linewidth:

  Numeric width of segments. Default is 0.3.

- add_points:

  Logical indicating whether to add the default sample points layer via
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html).
  Default is `TRUE`. Set to `FALSE` to omit it and add a custom one
  instead.

- point_size:

  Numeric size of points. Default is 1.5.

- ...:

  Currently unused.

## Value

A `ggplot` object, so that it can be further customised (e.g.
`plot(sp) + ggplot2::theme_bw()`). The default data of the plot contains
the PCoA coordinates of each sample, in columns `PCoA.1` and `PCoA.2`
(the plotted axes, in order).

## Details

The PCoA is computed on the full dissimilarity matrix, on the link scale
(see `k_pcoa` and `link` in
[`sparse_diss`](https://dhercar.github.io/rrGDM/reference/sparse_diss.md)).
