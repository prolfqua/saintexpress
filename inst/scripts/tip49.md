# TIP49 example data in inst/extdata/tip49

## Source

The data are a subset of TIP49, the example input shipped with SAINTexpress 3.6.3 (`SAINTexpress_v3.6.3__2018-03-09/example input files/TIP49/`, from the SAINTexpress distribution by Hyungwon Choi and colleagues). TIP49 is an AP-MS experiment on human protein complexes containing TIP49a/TIP49b (RUVBL1/RUVBL2), with spectral counts: 27 baits in 35 purifications, 9 negative-control purifications, 1207 preys.

Sardiu ME, Cai Y, Jin J, Swanson SK, Conaway RC, Conaway JW, Florens L, Washburn MP (2008). Probabilistic assembly of human protein interaction networks from label-free quantitative proteomics. PNAS 105(5), 1454-1459. doi:10.1073/pnas.0706983105

The subset keeps the package checks fast: 4 baits (ACTR8, ZNHIT1, INO80B, VPS72) in 7 purifications, all 9 controls, and the 78 preys detected in at least 3 of the 7 bait purifications.

## Files

Input (tab-separated, no header; every line is copied unchanged from the TIP49 files):

- `inter.dat`: IP id, bait id, prey id, spectral count (715 rows)
- `prey.dat`: prey id, protein length (78 rows)
- `bait.dat`: IP id, bait id, `T` (test) or `C` (control) (16 rows)

Output of the native SAINTexpress 3.6.3 binaries on these inputs, unchanged:

- `list_spc.txt`: `SAINTexpress-spc inter.dat prey.dat bait.dat`
- `list_int.txt`: `SAINTexpress-int inter.dat prey.dat bait.dat` (the spectral counts are used as intensities)

## How the subset was made

In R, from the TIP49 directory of the SAINTexpress distribution:

```r
rd <- function(f) read.delim(f, header = FALSE, colClasses = "character")
inter <- rd("inter.dat")
prey <- rd("prey.dat")
bait <- rd("bait.dat")
bait <- bait[bait$V2 %in% c("ACTR8", "ZNHIT1", "INO80B", "VPS72") | bait$V3 == "C", ]
inter <- inter[inter$V1 %in% bait$V1, ]
n_ips <- table(inter$V3[inter$V1 %in% bait$V1[bait$V3 == "T"]])
preys <- names(n_ips)[n_ips >= 3]
inter <- inter[inter$V3 %in% preys, ]
prey <- prey[prey$V1 %in% preys, ]
wr <- function(x, f) {
  write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
wr(inter, "inter.dat")
wr(prey, "prey.dat")
wr(bait, "bait.dat")
```

## How the native output was produced

The official precompiled Linux64 binaries from `SAINTexpress_v3.6.3__2018-03-09/Precompiled_binaries/Linux64/` were run with default options in a Debian container:

```sh
# Dockerfile: debian:bookworm-slim with libstdc++6, the two binaries copied to /usr/local/bin
docker build --platform linux/amd64 -t saintexpress:3.6.3 .
docker run --rm --platform linux/amd64 -v "$PWD:/work" -w /work saintexpress:3.6.3 \
  SAINTexpress-spc inter.dat prey.dat bait.dat   # writes list.txt, renamed to list_spc.txt
docker run --rm --platform linux/amd64 -v "$PWD:/work" -w /work saintexpress:3.6.3 \
  SAINTexpress-int inter.dat prey.dat bait.dat   # writes list.txt, renamed to list_int.txt
```

## MD5 checksums

```
e23b9042d734a434747c9b42184b8085  inter.dat
00a681219ba925bd27e5806f9d92c1d8  prey.dat
b1ceaff8ff29ba085c0876d927cd6e97  bait.dat
9b0a3533691c64870519c50111b7fcf0  list_spc.txt
fddee63756246230eb8be836f130b560  list_int.txt
```
