# TIP49 example data in inst/extdata/tip49

## Source

The TIP49 dataset is the example input shipped with SAINTexpress 3.6.3 (`SAINTexpress_v3.6.3__2018-03-09/example input files/TIP49/`, from the SAINTexpress distribution by Hyungwon Choi and colleagues). It is an AP-MS experiment on human protein complexes containing TIP49a/TIP49b (RUVBL1/RUVBL2), with spectral counts: 27 baits in 35 purifications, 9 negative-control purifications, 1207 preys.

Sardiu ME, Cai Y, Jin J, Swanson SK, Conaway RC, Conaway JW, Florens L, Washburn MP (2008). Probabilistic assembly of human protein interaction networks from label-free quantitative proteomics. PNAS 105(5), 1454-1459. doi:10.1073/pnas.0706983105

## Files

Input, copied unchanged from the SAINTexpress distribution (tab-separated, no header):

- `inter.dat`: IP id, bait id, prey id, spectral count (7855 rows)
- `prey.dat`: prey id, protein length (1207 rows)
- `bait.dat`: IP id, bait id, `T` (test) or `C` (control) (44 rows)

Output of the native SAINTexpress 3.6.3 binaries on these inputs, unchanged:

- `list_spc.txt`: `SAINTexpress-spc inter.dat prey.dat bait.dat`
- `list_int.txt`: `SAINTexpress-int inter.dat prey.dat bait.dat` (the spectral counts are used as intensities)

## How the native output was produced

The official precompiled Linux64 binaries from `SAINTexpress_v3.6.3__2018-03-09/Precompiled_binaries/Linux64/` were run with default options in a Debian container:

```sh
# Dockerfile: debian:bookworm-slim with libstdc++6, the two binaries copied to /usr/local/bin
docker build -t saintexpress:latest .
docker run --rm --platform linux/amd64 -v "$PWD:/work" -w /work saintexpress:latest \
  SAINTexpress-spc inter.dat prey.dat bait.dat   # writes list.txt, renamed to list_spc.txt
docker run --rm --platform linux/amd64 -v "$PWD:/work" -w /work saintexpress:latest \
  SAINTexpress-int inter.dat prey.dat bait.dat   # writes list.txt, renamed to list_int.txt
```

## MD5 checksums

```
b0a885df9188b3088bfd9b19dbaabe5f  inter.dat
0c6cc1e5869415fbf4c6407766be7368  prey.dat
e78aa0eb84efacbd87f305bbce3e67d5  bait.dat
7370f6b0b599730626f746a6cda2c387  list_spc.txt
31b5b206aa3b9defa24436dc50188d54  list_int.txt
```
