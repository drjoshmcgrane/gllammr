# Weighted isotonic regression (pool-adjacent-violators)

Returns the nondecreasing vector minimizing sum(w \* (y - x)^2), which
is the constrained M-step maximizer for binomial proportions and normal
means under a monotone-classes restriction (Croon 1990).

## Usage

``` r
.pava_weighted(y, w)
```
