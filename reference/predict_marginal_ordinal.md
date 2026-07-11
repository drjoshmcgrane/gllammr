# Internal function for marginal ordinal predictions

Population-averaged category probabilities by Monte Carlo over the
estimated random-effects distributions of every term (fresh draws per
replicate; the same draw applies to all members of a group only in
expectation, which is what the marginal quantity requires).

## Usage

``` r
predict_marginal_ordinal(object, n_sim = 1000, newdata = NULL)
```
