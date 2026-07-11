# Latent structure comparison: categorization, ordering or quantification

Fits the six latent structure models of Torres Irribarra & Diakow's
model selection framework to a binary response matrix and returns a
comparison table. The models progressively constrain the latent
structure, decomposing the order and scale assumptions:

## Usage

``` r
latent_structure_comparison(
  Y,
  nclass = NULL,
  item_order = c("auto", "columns"),
  n_starts = 5,
  control = list()
)
```

## Arguments

- Y:

  Binary response matrix (persons x items)

- nclass:

  Number of latent classes for the class-based models. The default
  `ceiling((ncol(Y) + 1) / 2)` is the smallest number at which the
  located class model is fit-equivalent to the semiparametric Rasch
  model (Lindsay et al. 1991).

- item_order:

  How to order items for the IIO/DM constraints: "auto" (default; by
  marginal proportion correct) or "columns" (the column order of `Y`)

- n_starts:

  Random starts for each class-based model (default 5)

- control:

  Control list passed to the fitters

## Value

An object of class `lca_structure_comparison`: a data frame with one row
per model (structure type, logLik, nominal parameter count, AIC, BIC)
plus the fitted models in `attr(, "fits")`.

## Details

- UN:

  Unconstrained latent class model - qualitative structure (differences
  of kind).

- MON:

  Ordered classes with class (person) monotonicity (Croon 1990).

- IIO:

  Ordered classes with invariant item ordering (item monotonicity).

- DM:

  Double monotonicity: both restrictions.

- LCR:

  Located latent classes (latent class Rasch; Lindsay, Clogg &
  Grego 1991) - a discrete quantitative structure.

- RM:

  The Rasch model with a normal latent distribution - a continuous
  quantitative structure.

Successive comparisons carry the framework's logic: UN vs MON/IIO asks
whether an ordering is tenable at all; the single-monotonicity models vs
DM asks whether persons and items share one proficiency progression; DM
vs LCR isolates the interval-scale (parameter separability) assumption;
LCR vs RM asks whether a continuous latent variable adds anything beyond
located classes. UN, MON, IIO, and DM share the same nominal parameter
count (inequality constraints do not reduce it), so their information
criteria differ only through fit; treat those comparisons descriptively
(chi-bar-square caveat).

## References

Torres Irribarra, D., & Diakow, R. Categorization, ordering or
quantification: selecting a latent variable model by comparing latent
structures.

## Examples

``` r
if (FALSE) { # \dontrun{
cmp <- latent_structure_comparison(resp, nclass = 4)
print(cmp)
} # }
```
