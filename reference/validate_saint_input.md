# Validate a SAINT input list

Checks that `si` is a list with exactly the elements `inter`, `prey`,
`bait`, and that each is a data frame.

## Usage

``` r
validate_saint_input(si)
```

## Arguments

- si:

  A named list with elements `inter`, `prey`, `bait`.

## Value

`TRUE` invisibly on success; otherwise stops with a message.
