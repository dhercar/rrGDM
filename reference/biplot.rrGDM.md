# Biplot for rrGDM objects

[`biplot`](https://rdrr.io/r/stats/biplot.html) method for objects of
class `rrGDM`, showing sample positions in two of the fitted latent
dimensions, optionally overlaid with environmental loadings and/or
species correlations with latent axes.

## Usage

``` r
# S3 method for class 'rrGDM'
biplot(
  x,
  k = c(1, 2),
  rotate = TRUE,
  type = c("conditional", "marginal"),
  col_by = NULL,
  shape_by = NULL,
  add_points = TRUE,
  env_load = TRUE,
  species_cor = FALSE,
  species_cor_n = 20,
  env_load_scale = 1,
  sp_load_scale = 1,
  text_size = 3,
  compare = FALSE,
  ...
)
```

## Arguments

- x:

  A fitted model of class `rrGDM`.

- k:

  Integer vector of length 2 indicating which latent dimensions to plot.
  Default is `c(1, 2)`.

- rotate:

  Logical indicating whether to apply varimax rotation (via
  [`varimax`](https://rdrr.io/r/stats/varimax.html)) before plotting.
  Default is `TRUE`. Varimax rotates the axes towards simple structure;
  it does not order them by variance explained. As in
  [`get_LV`](https://dhercar.github.io/rrGDM/reference/get_LV.md), the
  rotation is estimated from the conditional coordinates (`L`), whatever
  the `type` plotted.

- type:

  Character string: either `'conditional'` (site scores including
  sample-level error, `L`) or `'marginal'` (site scores as determined by
  environmental loadings only, `L_m`). Default is `'conditional'`.

- col_by:

  Optional character string naming a variable in `x$X` used to colour
  points. Ignored if `add_points = FALSE`.

- shape_by:

  Optional character string naming a variable in `x$X` used to shape
  points. Ignored if `add_points = FALSE`.

- add_points:

  Logical indicating whether to add the default sample points layer via
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html).
  Default is `TRUE`. Set to `FALSE` to omit this layer entirely and add
  a fully custom one instead. The returned plot's default data already
  contains the rotated latent coordinates together with all columns used
  to fit the model, so
  `biplot(model, add_points = FALSE) + ggplot2::geom_point(ggplot2::aes(fill = pH), shape = 21)`
  works directly, with more control over aesthetics, shapes, and
  legends.

- env_load:

  Logical indicating whether to add environmental loadings as arrows.
  Default is `TRUE`.

- species_cor:

  Logical indicating whether to add species labels positioned by their
  (Spearman) correlation with each plotted latent axis. Default is
  `FALSE`.

- species_cor_n:

  Optional integer. If `species_cor = TRUE`, only the `species_cor_n`
  species with the strongest correlation are labelled. Useful when the
  community matrix has hundreds or thousands of species. Default is
  `20`. Set to `NULL` to label all species.

- env_load_scale:

  Numeric scale factor applied to environmental loading arrows. Default
  is 1.

- sp_load_scale:

  Numeric scale factor applied to species-correlation labels. Default is
  1.

- text_size:

  Numeric text size used for species labels. Default is 3.

- compare:

  Logical indicating whether to overlay marginal scores on top of
  conditional scores, connected by a segment. Only valid when
  `type = 'conditional'`. Default is `FALSE`.

- ...:

  Currently unused.

## Value

A `ggplot` object. Unlike
[`biplot.default`](https://rdrr.io/r/stats/biplot.html), the object is
returned so it can be further customised (e.g.
`biplot(model) + ggplot2::theme_bw()`).

## Details

The function extracts fitted latent coordinates, optionally rotates them
via varimax rotation to improve interpretability, and plots the two
latent dimensions requested. Environmental loadings and species
correlations, if requested, are rotated using the same rotation matrix
so that they remain aligned with the sample scores.

## Examples

``` r
if (FALSE) { # \dontrun{
com <- gllvm::microbialdata$Y
env <- gllvm::microbialdata$Xenv
env$pH_s <- scale(env$pH)

m1 <- rrGDM(Y = com,
            X = env,
            k = 2,
            seed = 1,
            formula = ~ pH_s + as.factor(Region),
            family = 'binomial')

biplot(m1,
       shape_by = 'Region',
       col_by = 'pH',
       type = 'conditional')

# Custom:
biplot(m1, add_points = FALSE) +
  ggplot2::geom_point(ggplot2::aes(fill = pH), shape = 21, size = 3) +
  ggplot2::scale_fill_viridis_c()

# Add species correlation with latent axes
biplot(m1, species_cor = TRUE, species_cor_n = 20)
} # }
```
