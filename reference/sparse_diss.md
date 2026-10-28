# Select a sparse set of sample pairs for model fitting

Selects a subset of sample pairs to be used for fitting
[`rrGDM`](https://dhercar.github.io/rrGDM/reference/rrGDM.md), such that
every sample is in exactly the same number of pairs (`knn + n_random`).
Each sample is paired with its `knn` nearest neighbours, which inform
its local position, and with `n_random` random samples, which ensure
that the global structure (including distant samples) is represented.
Selection is based on the observed dissimilarities only.

## Usage

``` r
sparse_diss(
  Y,
  knn = 10,
  n_random = 3,
  max_zeros = Inf,
  dist_method = "bray",
  trans = function(x) x,
  binary = FALSE,
  k_pcoa = 2,
  link = "identity",
  seed = NULL
)
```

## Arguments

- Y:

  A community matrix or data frame (samples x species), or an object of
  class [`dist`](https://rdrr.io/r/stats/dist.html). If a `dist` object
  is supplied, `dist_method`, `trans` and `binary` are ignored.

- knn:

  Integer. Number of nearest neighbours per sample. Default is 10.

- n_random:

  Integer. Number of random pairs per sample. Default is 3.

- max_zeros:

  Integer or `Inf`. Maximum number of pairs with a dissimilarity of
  exactly 0 per sample. Default is `Inf`.

- dist_method:

  Dissimilarity index: a method of
  [`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  (e.g. `'bray'`, `'jaccard'`, `'euclidean'`); `'sorensen'` is
  Bray-Curtis on presence-absence data. Ignored if `Y` is a `dist`
  object. Default is `'bray'`.

- trans:

  Transformation applied to `Y` before calculating dissimilarities.
  Default is no transformation.

- binary:

  Logical indicating whether `Y` should be converted to presence-absence
  before calculating dissimilarities. Default is FALSE.

- k_pcoa:

  Integer (at least 2) or `NULL`. Number of PCoA axes of the full
  dissimilarity matrix stored for 2D representation. Use `NULL` to skip
  it. Default is 2.

- link:

  Character string: `'identity'` (default), `'negexp'` or `'neggaus'`.
  Link function of the model to be fitted (see
  [`rrGDM`](https://dhercar.github.io/rrGDM/reference/rrGDM.md)). The
  PCoA is computed on dissimilarities transformed to the scale of latent
  distances (e.g. `-log(1 - diss)` for `'negexp'`), after setting
  dissimilarities above 0.99 to 0.99. This does not affect which pairs
  are selected.

- seed:

  An integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) at the start, so it
  changes the state of R's random number generator. Default is NULL (no
  seed is set).

## Value

An object of class `sparse_diss`: a list including `edges`, a two-column
integer matrix (`s1`, `s2`) of the pairs kept, which
[`rrGDM`](https://dhercar.github.io/rrGDM/reference/rrGDM.md) uses to
subset the dissimilarities (argument `sparse`), and `pairs`, a data
frame with the pairs, their dissimilarity (`diss`) and their type
(`'knn'`, `'random'` or `'zero'`). Each pair is stored once, with
`s1 < s2`. `degree` gives the number of pairs per sample. The other
components are the sample names (`id`), the number of samples and pairs
(`n_samples`, `n_pairs_full`, `n_pairs`, `n_zero_pairs`), the settings
used (`knn`, `n_random`, `n_per_sample`, `max_zeros`, `link`,
`dist_method`, `binary`), the number of connected groups of samples
(`n_components`), the PCoA used by
[`plot.sparse_diss`](https://dhercar.github.io/rrGDM/reference/plot.sparse_diss.md)
(`pcoa`, `NULL` if `k_pcoa = NULL`) and the matched call (`call`).

## Details

Pairs are sorted by dissimilarity and kept if both samples have fewer
than `knn` pairs (and fewer than `max_zeros` for zero dissimilarities).
The remaining pairs up to `knn + n_random`, are filled with random pairs
among samples with free slots. Because the number of pairs per sample is
fixed, a sample may be paired with a sample that is not among its own
`knn` nearest neighbours. Pairs with a dissimilarity of exactly 0
(identical samples) are only added in the nearest-neighbour step, never
in the random step.

One sample has one pair fewer if both the number of samples and
`knn + n_random` are odd. A warning is given if some samples fall short
or if the pairs do not connect all samples (increase `n_random`).

## See also

[`plot.sparse_diss`](https://dhercar.github.io/rrGDM/reference/plot.sparse_diss.md),
[`print.sparse_diss`](https://dhercar.github.io/rrGDM/reference/print.sparse_diss.md)
