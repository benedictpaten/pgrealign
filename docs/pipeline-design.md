# Realigning regions of a pangenome graph: pipeline design (draft)

**Status:** agreed 2026-10-02; the build has started.

**Decisions (2026-10-02):**
- **Scope:** every catalogued tandem repeat whose spanning alleles are not all identical.
- **Catalogue:** CHM13 RepeatMasker (Simple_repeat, Low_complexity, Satellite), rebuilt reproducibly
  with no truth input.
- **Disk:** the user frees space for the full-graph stages.
- **Packaging:** a Python package with a CLI, plus WDL for Toil and Terra.

## Goal

Take a complete pangenome graph: a GBZ with every haplotype as a path, such as the ~460-haplotype HPRC
v2.1 Minigraph-Cactus graph. Then:

1. find every region of a given kind (here every STR/VNTR);
2. extract each region's haplotype sequences, with each path's identity;
3. realign them by an established protocol (here the full-panel medoid star, results 4x);
4. replace each region's subgraph with the graph induced from the new alignment, re-threading every
   path;
5. write an updated graph that can be indexed and tested.

The extract, realign and replace steps must not know they are working on tandem repeats. A different
BED of regions, or a different aligner, should drop in.

## What exists, and what does not carry over

The chr20 and chr6 experiments (results/stage4_*.md) proved the protocol on a 34-haplotype graph.
Their tools do not scale to the full graph:

| step | today | problem on the full graph |
|---|---|---|
| catalogue | RepeatMasker + truth TRF spans, in a session scratchpad (lost) | must be reproducible, committed, and free of truth |
| anchors | chosen on the hap32 graph (truth rule or hap32 snarls) | they do not separate the full graph: 607/660 chr20 regions have walks touching only one anchor |
| extraction | one `gbz-base query` per region (47 s, 11 GB each), anonymous `hprc#k` records, non-spanning walks dropped | weeks genome-wide; no path identity, so nothing can be written back |
| realign | `iterate.py st_medoid`, inside a chr20 test harness | needs a clean API, fragment handling, and size limits with a fallback |
| replace | `patch_contig.py` on a per-contig hap32 GFA, paths matched by hap32 name | full-graph W lines, ~460 haplotypes, fragments, gref cover, IDs, disk |

## Design: six stages with file contracts

Each stage is a subcommand of one CLI and reads and writes files. That lets any stage be re-run,
swapped out or parallelised. Per-contig stages stream the graph once; per-region stages are
independent jobs.

```
targets.bed --> [1 regions] --> regions.tsv --> [2 extract] --> <region>/ packages
                (per contig)                    (per contig, one streaming pass)
<region>/alleles.fa --> [3 realign] --> msa.fa --> [4 induce] --> region.gfa + row walks
                        (per region, pluggable)   (per region, general MSA -> graph)
graph + all <region>/ results --> [5 replace] --> updated contig GFA/GBZ --> [6 merge/verify]
                                  (per contig, streaming)
```

0. **`catalog`** (tandem-repeat specific, optional). Builds `targets.bed` from an annotation
   (RepeatMasker Simple_repeat/Low_complexity/Satellite), with class, period and motif. Any other BED
   of target intervals works in its place.
1. **`regions`** (per contig).
   - Places each target on the reference path. Its anchors are the boundaries of the smallest run of
     consecutive snarls, in the full graph's own snarl tree, that encloses the target ± pad.
   - Merges overlapping regions and applies the size limit, recording why each target was dropped.
   - Output `regions.tsv`: id, anchors, reference span, member targets, status.
2. **`extract`** (per contig, one streaming pass over the walks, with no per-region queries).
   - For every region it records each walk's run between the anchors with its identity: sample,
     haplotype, contig, start and step range.
   - It classifies each run: spanning, prefix, suffix, internal, or visiting an anchor twice.
   - It writes a package `<region>/`:
     - `alleles.fa`: distinct spanning sequences, with weights;
     - `fragments.fa`;
     - `occurrences.tsv`: each run mapped to an allele or fragment;
     - `region.json`.
   - The selection policy (for example "the alleles differ") is applied here, from the measured
     lengths.
3. **`realign`** (per region, pluggable).
   - Input: `alleles.fa` with weights, plus `fragments.fa`. Output: `msa.fa`, one row per allele and
     per fragment.
   - The default method is the medoid star: centre by weighted 15-mer Jaccard medoid; abPOA global
     pairwise; left normalisation; insertion slots merged by abPOA.
   - Fragments are aligned to the centre with free end gaps and become rows with end gaps, so every
     path still spells its sequence.
   - Size tiers keep long regions from running out of memory. A region that fails or is too big is
     marked `keep`.
4. **`induce`** (per region, general). Turns the MSA into a column graph, unchops it, chops nodes to
   at most 1,024 bp, and writes each row's walk. It checks that every row spells its sequence.
5. **`replace`** (per contig, streaming). Writes the contig with each region's interior nodes
   replaced:
   - The interior nodes are the union of the runs' interior nodes, checked to be used by no other
     path.
   - New node IDs are fresh, above the genome-wide maximum, so contigs never clash.
   - Every run is rewritten to its row's walk, and its new spelling is checked against the old.
   - A region that fails any check is left as it was, with the reason recorded.
6. **`merge` / `verify`.**
   - Assembles the contigs into one GBZ, keeping the reference metadata.
   - Optionally checks every path's full sequence by md5.
   - Writes a report: regions replaced or kept and why, and nodes and edges before and after.

**Packaging.** A Python package in this repository (`pip install -e .`) with one console script.
- Each stage is a subcommand, so a stage can be run on its own.
- A `run` subcommand chains the stages over a contig list, with resume and `--jobs`.
- External tools (vg, abPOA) are found on PATH, and their versions are recorded in every output.
- Each stage gets unit tests on small fixture graphs. These cover fragments, reverse walks, nested
  regions and shared boundaries.

**Evaluation is not part of the pipeline.** HG002 calling, truvari and edit distances stay in
`tools/` and `work/`, and read the pipeline's outputs.

## Things the full graph forces us to decide

- **The gref_CHM13 cover** (261,708 paths in the eval graph). Its fragments end inside regions. Strip it
  before `replace` and rebuild it afterwards, or work on the non-gref graph.
- **Fragments.** Most regions have some (chr20: 586 of 636 have at least one non-spanning walk). They
  must be threaded, not refused, or most regions stay unpatched.
- **Huge regions** (centromeric and satellite arrays, merged chains over 100 kb). These are `keep` by
  policy and reported.
- **Disk.** The full chr6 W-line text is about 24 GB. `vg gbwt -G` needs the output GFA as a file.

## Validation on chr6 (held out)

- **V1. The protocol.** The medoid star has run only on chr20, where its gain over abPOA-on-hap32 was
  not significant on its own.
  - Run it on chr6 through the new `regions` + `extract` stages on the HG002-free full graph.
  - Project onto hap32, patch the chr6 hap32 graph, re-map, and call HG002.
  - Compare with the MC and abPOA-hap32 chr6 arms already scored (SV F1 0.5894 and 0.6425), using the
    per-region bootstrap.
- **V2. The pipeline on the full graph.** Run all six stages on chr6 of the full eval graph. Then:
  - check that every path spells its original sequence;
  - count regions replaced and kept, by reason;
  - check pair consistency on a sample of regions.
- **V3. End to end.** Sample HG002 haplotypes from the updated chr6 graph with `vg haplotypes`, then
  map and call. Compare with the same chain on the original graph. This needs distance, r-index and
  .hapl indexes and HG002 k-mer counts, so it comes after V2.

## Build order

1. **Package skeleton and `catalog`.** `pyproject.toml`, the `pgrealign` CLI, a provenance manifest in
   every output, and stdlib unittest. Then the CHM13 catalogue, committed with the script that builds
   it.
2. **Graph preparation.** A per-contig full graph without the gref cover (`vg gbwt -R gref_CHM13`,
   `vg chunk -C`), its snarls, and reference positions. Disk and time are measured on chr6.
3. **`regions` and `extract`.**
   - Anchors come from full-graph snarls.
   - Extraction is one streaming pass, keeping path identity and fragments.
   - Check: at the chr20 anchors used before, the multiset of spanning sequences must equal
     `hprc.fa.gz`.
4. **`realign` and `induce`.**
   - The medoid star moves out of `tools/iterate.py`, gaining fragment rows and size tiers.
   - Check: the chr20 `st_medoid` MSAs are reproduced from the same inputs.
5. **V1 on chr6.**
   - The full panel, plus the hap32 alleles at full-graph anchors, is realigned and projected onto
     hap32.
   - It is patched into the chr6 hap32 graph with `tools/patch_contig.py`, then mapped, called and
     scored against the MC and abPOA-hap32 arms.
6. **`replace`, `merge`, `verify` on the full chr6 graph (V2).**
7. **WDL and a container, and the user guide.**
8. **V3:** `vg haplotypes` end to end.
