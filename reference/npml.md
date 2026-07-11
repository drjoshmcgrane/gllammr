# Nonparametric maximum likelihood integration

Integration specification for
[`gllamm`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md):
replace the normal latent distribution with k estimated mass points and
masses (NPML; Aitkin 1999). Equivalent to
[`fit_npml`](https://drjoshmcgrane.github.io/gllammr/reference/fit_npml.md).

## Usage

``` r
npml(k = 2)
```

## Arguments

- k:

  Number of mass points (default 2)

## Value

An object of class `gllamm_integration`

## Examples

``` r
if (FALSE) { # \dontrun{
gllamm(y ~ x + (1 | g), data = d, family = binomial(),
       integration = npml(2))
} # }
```
