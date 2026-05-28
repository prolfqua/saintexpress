# saintexpress

`saintexpress` provides a pure-R implementation of the SAINTexpress scoring
engines for spectral counts and intensities.

This package is intentionally focused on the R implementation. It does not ship
or run native SAINTexpress binaries. Native execution is handled by the separate
`saintexpressbin` package, and prolfqua integration is handled by
`prolfquasaint`.

## Usage

Prepare a SAINT input list with three data frames named `inter`, `prey`, and
`bait`:

```r
si <- list(
  inter = inter,
  prey = prey,
  bait = bait
)
```

Validate the input shape:

```r
saintexpress::validate_saint_input(si)
```

Run the spectral-count engine:

```r
result <- saintexpress::run_saint(si, mode = "spc")
```

Run the intensity engine:

```r
result <- saintexpress::run_saint(si, mode = "int")
```

`result` is a data frame in the SAINTexpress `list.txt` output shape.

## Optimizers

The default optimizer is base R:

```r
saintexpress::run_saint(si, mode = "spc", optimizer = "base")
```

The optional `nloptr` optimizer can be used when the `nloptr` package is
installed:

```r
saintexpress::run_saint(si, mode = "spc", optimizer = "nloptr")
```

## Testing Scope

This package uses small synthetic inputs for unit tests. TIP49 reference
fixtures and native-vs-R comparison tests live in `prolfquasaint`, where both
`saintexpress` and `saintexpressbin` are integrated.
