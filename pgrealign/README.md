# pgrealign

`pgrealign` realigns regions of a pangenome graph. It:
1. extracts each region's haplotype sequences, with each path's identity;
2. realigns them with a pluggable multiple aligner;
3. splices the graph induced from the new alignment back in, re-threading every path.

It was built to replace the tandem-repeat (STR/VNTR) alignments of the HPRC Minigraph-Cactus graph. The
extract, realign and replace stages do not know they are working on repeats. Give them any BED of
regions and any aligner.

The design and its validation are in [docs/pipeline-design.md](../docs/pipeline-design.md).

**Status:** every stage works end to end on chr20 of the HPRC v2.1 eval graph. See
[results/full_graph_chr20.md](../results/full_graph_chr20.md): SV F1 +0.040 for short reads and +0.084 for ONT
(refined), significant per region, with small variants unchanged. WDL/Toil packaging is not written yet.

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
| `prepare` | cut a contig from a whole-genome GBZ, drop the gref cover, snarls, reference index | done |
| `regions` | anchor each target in the graph's own snarl tree; merge; size policy | done |
| `extract` | one streaming pass: every path's run through every region, with identity and fragments | done |
| `realign` | per region: alleles (+ fragments) -> MSA; default the full-panel medoid star | done |
| `induce` | per region: MSA -> column graph + one walk per row | done |
| `replace` | per contig: splice new regions in, re-thread every path, check spellings | done |
| `union` / `project` | realign on one panel, apply to a graph of other haplotypes | done |
| `verify` | every path spells what it did before | done |
| `merge` | assemble contigs into one GBZ (needs a global renumbering) | planned |

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
