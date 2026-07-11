# Expand a nested random-effects term into one term per level

(x \| a/b/c) becomes the lme4-equivalent (x \| a) + (x \| a:b) + (x \|
a:b:c). Each expanded term carries an interaction grouping over the
levels above it.

## Usage

``` r
expand_nested_random_term(parsed_term)
```

## Arguments

- parsed_term:

  Output of parse_random_term()

## Value

List of parsed terms (length 1 for non-nested input)
