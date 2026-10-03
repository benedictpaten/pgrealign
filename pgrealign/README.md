# pgrealign

`pgrealign` realigns regions of a pangenome graph. It:
1. extracts each region's haplotype sequences, with each path's identity;
2. realigns them with a pluggable multiple aligner;
3. splices the graph induced from the new alignment back in, re-threading every path.

It was built to replace the tandem-repeat (STR/VNTR) alignments of the HPRC Minigraph-Cactus graph. The
extract, realign and replace stages do not know they are working on repeats. Give them any BED of
regions and any aligner.

The design and its validation are in [docs/pipeline-design.md](../docs/pipeline-design.md).

**Status:** under construction. Only `catalog` exists so far.

## Install

```
pip install -e .            # from the repository root; Python >= 3.9, standard library only
```

External tools, found on `PATH`:
- `vg` and `gbz-base`, for the graph stages;
- `abpoa` (>= 1.5), for the default aligner;
- a C compiler, which builds `trfind` on first use.

## Stages

| stage | what it does | status |
|---|---|---|
| `catalog` | tandem-repeat catalogue of the reference (RepeatMasker + trfind) | done |
| `regions` | anchor each target in the graph's own snarl tree; merge; size policy | planned |
| `extract` | one streaming pass: every path's run through every region, with identity and fragments | planned |
| `realign` | per region: alleles (+ fragments) -> MSA; default the full-panel medoid star | planned |
| `induce` | per region: MSA -> column graph + one walk per row | planned |
| `replace` | per contig: splice new regions in, re-thread every path, check spellings | planned |
| `merge` / `verify` | assemble the contigs into one GBZ; report | planned |

### catalog

```
pgrealign catalog --repeatmasker chm13v2.0_RepeatMasker.bed \
    --reference chr1.fa chr2.fa ... chrX.fa -j 6 -o chm13v2.0_tr.bed.gz
```

The catalogue merges RepeatMasker's Simple_repeat, Satellite and Low_complexity records with the
repeats that `trfind` finds in the reference: a seed-based finder for periods 7-2,000 bp. Records
within 50 bp of each other become one region. RepeatMasker alone is not enough, because its longest
Simple_repeat motif is 30 bp.

[catalogs/chm13v2.0_tr.bed.gz](../catalogs/chm13v2.0_tr.bed.gz) is the CHM13 v2.0 catalogue: 845,357
regions on chr1-22 and chrX, 265.6 Mb. It takes 26 s and 5.5 GB to build, and its `.manifest.json`
records exactly how it was built.

## Tests

```
python3 -m unittest discover -s tests -v
```
