# Validate a SAINT input list

Checks that `si` is valid input for
[`run_saint()`](https://prolfqua.github.io/saintexpress/reference/run_saint.md):
if `validate_saint_input()` passes,
[`run_saint()`](https://prolfqua.github.io/saintexpress/reference/run_saint.md)
accepts the same `si` and `mode`. The checks are the ones listed in the
section "Input tables" below, applied after the columns have been
matched by name or position.

## Usage

``` r
validate_saint_input(si, mode = c("spc", "int"))
```

## Arguments

- si:

  A named list with the data frames `inter`, `prey` and `bait`, in any
  order. See the section "Input tables".

- mode:

  `"spc"` for spectral counts, `"int"` for intensities. Sets the allowed
  range of `quant`.

## Value

`TRUE` invisibly on success; otherwise stops with a message naming the
table and column at fault.

## Input tables

Columns are matched by name when all required names are present;
otherwise by position, in the order listed (as in the headerless native
input files). A column that already carries a required name must then be
at that name's position. Further columns are ignored.

`inter`, one row per prey detected in an IP (purification):

|  |  |  |
|----|----|----|
| Column | Type | Content |
| `ipId` | character | IP identifier, as in `bait$ipId` |
| `baitId` | character | Bait of the IP, as in `bait$baitId` |
| `preyId` | character | Prey identifier, as in `prey$preyId` |
| `quant` | numeric | Spectral count, whole number \>= 0 (`"spc"`), or intensity, \> 0 (`"int"`) |

`prey`, one row per prey:

|  |  |  |
|----|----|----|
| Column | Type | Content |
| `preyId` | character | Unique prey identifier |
| `preyLength` | numeric | Protein length (optional, not used in scoring) |
| `preyGeneId` | character | Gene name, reported as `PreyGene` (optional, defaults to `preyId`) |

The `prey` table is matched by name when it has both `preyId` and
`preyGeneId`. Otherwise the first column is `preyId`; with two columns
the second is `preyLength` if numeric, else `preyGeneId`; with three or
more, the second is `preyLength` and the third `preyGeneId`.

`bait`, one row per IP:

|          |           |                                                   |
|----------|-----------|---------------------------------------------------|
| Column   | Type      | Content                                           |
| `ipId`   | character | Unique IP identifier                              |
| `baitId` | character | Bait name; replicate IPs share it                 |
| `CorT`   | character | `"T"` for a test IP, `"C"` for a negative control |

Every `inter$ipId` must be in `bait`, every `inter$preyId` in `prey`,
each prey at most once per IP, and `inter$baitId` must agree with
`bait`. At least one test and one control IP are needed.

## Examples

``` r
si <- list(
  inter = data.frame(
    ipId = c("IP1", "IP2", "IP3"),
    baitId = c("BaitA", "BaitA", "Ctrl"),
    preyId = c("Prey1", "Prey1", "Prey1"),
    quant = c(12, 9, 1)
  ),
  prey = data.frame(preyId = "Prey1", preyLength = 500, preyGeneId = "GENE1"),
  bait = data.frame(
    ipId = c("IP1", "IP2", "IP3"),
    baitId = c("BaitA", "BaitA", "Ctrl"),
    CorT = c("T", "T", "C")
  )
)
validate_saint_input(si, mode = "spc")
```
