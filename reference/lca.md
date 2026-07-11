# Latent Class Family for Finite Mixture Models

Create a family object for fitting latent class models through the
unified
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
interface. The response is a matrix of binary manifest variables passed
as the first argument of
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md).

## Usage

``` r
lca(nclass = 2, ordering = "none")
```

## Arguments

- nclass:

  Number of latent classes (default 2)

- ordering:

  Class order restriction passed to
  [`fit_lca`](https://drjoshmcgrane.github.io/gllammr/reference/fit_lca.md):
  "none" (default), "increasing", or a list/matrix of class pairs
  defining a partial order

## Value

A family object of class `lca_family`

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- gllamm(indicator_matrix, family = lca(nclass = 3))
fit_ord <- gllamm(indicator_matrix,
                  family = lca(nclass = 3, ordering = "increasing"))
} # }
```
