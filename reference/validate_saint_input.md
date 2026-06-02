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

## Examples

``` r
si <- list(
  inter = data.frame(ipId = "IP1", baitId = "BaitA", preyId = "Prey1"),
  prey = data.frame(preyId = "Prey1"),
  bait = data.frame(ipId = "IP1", baitId = "BaitA", CorT = "T")
)
validate_saint_input(si)
```
