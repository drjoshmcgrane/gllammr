# Weighted isotonic regression over a partial order

Minimizes sum(w \* (y - x)^2) subject to x\[a\] \<= x\[b\] for every row
(a, b) of `edges` - the constrained M-step maximizer for binomial
proportions and normal means under a partially ordered classes
restriction. A chain is solved exactly by pool-adjacent-violators; a
general DAG by Dykstra's cyclic projection algorithm (each constraint
set is a half-space whose weighted projection is a two-point pool),
which converges to the exact projection. Problems here are tiny
(length(y) = number of classes), so iteration cost is negligible.

## Usage

``` r
.isotonic_poset(y, w, edges)
```
