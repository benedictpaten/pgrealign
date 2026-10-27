# pgrealign reference

This reference describes the `pgrealign` command line, the files its stages read and write, and the
Python functions through which another program, such as a graph builder, can call the realignment
method directly. [method.md](method.md) explains how each stage works and why; this document says
what each stage reads and writes, what its options control, and when to change them.
`pgrealign <stage> --help` lists the stage's options with most of their defaults, and
[whole-genome.md](whole-genome.md) describes how to build a realigned graph of a whole genome from
per-contig runs.

1. [Overview](#overview)
2. [The command line](#the-command-line)
3. [File formats](#file-formats)
4. [The Python API](#the-python-api)
5. [Known gaps](#known-gaps)

## Overview

pgrealign replaces chosen regions of a pangenome graph with graphs built from a new multiple
alignment of the same sequences. Every haplotype spells the same sequence afterwards as before,
apart from letter case ([Known gaps](#known-gaps)); only the way the sequences share nodes changes.
The stages are subcommands of one program, `pgrealign`, that pass their results to one another in
files, so each can be rerun, replaced or run in parallel on its own. Most stages work on one
contig, and `merge` joins the realigned contigs into one graph.

### Installation and tools

The [README](../README.md#install) gives the requirements: Python, vg, abPOA 1.5 or later, and a C
compiler for `catalog --reference` and `--strf`, unless `--trfind` and `--strf` name trfind and strf
binaries built elsewhere. The stages find `vg` and `abpoa` on `PATH` unless `--vg` and
`--abpoa` name them. pgrealign does not check vg's version; the vg in use must support the commands
pgrealign runs: `vg convert -f` (with `-H`), `vg chunk --gbz --contig`, `vg gbwt` (with `-Z`, `-R`,
`--set-reference`, `--tags`, `-x`, `-G`, `-f` and `--gbz-format`), `vg snarls -T -P`, `vg view -Rj`,
`vg stats -r`, `-z`, `-l`, and `vg paths -M` and `-F`.

### Graphs

pgrealign reads graphs as GFA. GFA is a text format with one record per line, the first character
naming its type: a node (an S line) has an integer id and a sequence, and an edge (an L line) joins
the ends of two nodes. GBZ is vg's compressed format, which holds the nodes and edges with a GBWT, a
compressed index of the graph's paths. `extract`, `replace` and `qc` read a file whose name ends in
`.gbz`, or in `.vg`, `.pg`, `.hg` or `.xg` for vg's other formats, through `vg convert -f`, which
writes it as GFA, and any other file as GFA, gzipped when its name ends in `.gz`. `verify` passes
every kind of file to vg.

### Walks

This section and those that follow it, up to
[The stages and their files](#the-stages-and-their-files), define pgrealign's terms;
[method.md](method.md) defines them at more length.

| term | meaning |
|---|---|
| sample, haplotype | one individual, and one of its copies of the genome |
| walk | one haplotype's route through the graph, a W line: steps such as `>12<13>14`, each a node visited forward (`>`) or in reverse (`<`); a haplotype assembled in pieces has one walk per piece |
| walk key | `sample#haplotype#contig#start`, which names a walk. vg keeps it unique in a GBZ; a GFA given to `extract` and `replace` must keep it unique too, since neither checks, and `replace` would confuse two walks that share a key |

### Paths and the reference

| term | meaning |
|---|---|
| path | in vg's sense, a walk or a P line, which is a named list of oriented nodes. `extract` and `replace` treat only W lines as walks; `verify` and `merge` compare paths of both kinds |
| reference sample, contig | the sample (`--ref-sample`, CHM13 unless set) whose walk on each of its sequences, or contigs, such as chr1, gives that contig's coordinates; `<ref>` stands for its name. In the graph a contig is the connected component that holds the reference's walk on it. Every coordinate in pgrealign's files is 0-based and half-open, as in BED |

### Snarls and targets

| term | meaning |
|---|---|
| snarl, chain | the part of the graph between two boundary nodes, in vg's snarl decomposition (a trivial snarl's boundary nodes are joined directly by an edge); and a series of snarls, each sharing a boundary node with the next, top-level when it lies inside no snarl |
| target | an interval of the reference to realign, from a catalogue or any BED file, under the contig names of the reference sample's W lines and in its coordinates |

### Regions and runs

| term | meaning |
|---|---|
| region, anchors | the part of the graph between two anchors, L and R: boundary nodes of one chain, on either side of one or more targets, each visited once, forward, by the reference, so that every walk that crosses the region passes through both. A region's id is its first target's |
| run | one walk's passage through one region: spanning (from one anchor to the other), prefix (from the walk's start, inside the region, to an anchor), suffix (from an anchor to the walk's end, inside the region), internal (an anchor-free walk, one that visits no anchor of any region that `regions` made `ok`, lying inside the region) or complex (any other, such as a U-turn) |
| interior | the nodes between a region's anchors that its recorded runs visit, which `replace` removes |

### Members

| term | meaning |
|---|---|
| allele | a distinct sequence that spanning runs spell, in the reference's direction; empty when L and R are adjacent |
| fragment, side | a distinct sequence of prefix, suffix or internal runs, together with the end of the region it touches: `left` (L), `right` (R) or `internal` |
| member, weight | an allele or a fragment, with an id `a1`, `a2`, ... or `f1`, `f2`, ... in the order `extract` meets them; and the number of runs that spell it |

### Packages and statuses

| term | meaning |
|---|---|
| package | what `extract` collects for one region: its anchors, its members and a status |
| status | what `regions`, `extract` and `realign` each say of a region. `replace` changes only a region that is `ok` in all three; every other keeps its sequence and topology, though under `replace`'s default `dense` numbering its nodes take new ids ([Node ids](#node-ids)) |

### Alignments

| term | meaning |
|---|---|
| MSA, column | a multiple sequence alignment: one row per member, its sequence with gap characters (`-`) inserted, every row the same length; and one position across the rows |
| unbanded, banded | a pairwise alignment fills a dynamic-programming matrix of one cell per pair of bases. Unbanded, it fills every cell and finds an optimal alignment under its scores; banded, it fills only a band of cells around its best path, which costs less but can miss the optimum |

### The medoid star

| term | meaning |
|---|---|
| star, centre | the way `realign` builds an MSA ([method.md](method.md#realign-the-star)): every other allele, and every fragment found in no allele (an unplaced fragment), is aligned to one allele, the centre, and the pairwise alignments are merged; a fragment found in an allele takes that allele's row |
| medoid | the centre `realign` chooses: the allele of positive weight whose weighted sum of k-mer distances to the other alleles is least, a k-mer being a substring of length k |

### The stages and their files

| stage | works on | reads | writes |
|---|---|---|---|
| `catalog` | a genome | RepeatMasker BED, reference FASTA | catalogue BED |
| `prepare` | a contig | whole-genome GBZ | `<contig>.gbz`, `<contig>.snarls.pb`, `<contig>.ref.tsv.gz` |
| `regions` | a contig | targets, snarls, reference index | `regions.tsv` |
| `extract` | a contig | `<contig>.gbz`, `regions.tsv` | `packages.jsonl.gz`, `occurrences.tsv.gz`, `interiors.tsv.gz`, `summary.json` |
| `realign` | a contig's packages, or one shard of them ([realign](#parallel-runs-and-resuming)) | `packages.jsonl.gz` | `msas.jsonl.gz` |
| `replace` | a contig | `<contig>.gbz`, `extract`'s directory, `msas.jsonl.gz` | realigned GFA, and GBZ if asked |
| `verify` | a contig | `<contig>.gbz` and the realigned GBZ | a JSON report |
| `qc` | a contig | `extract`'s directory, `msas.jsonl.gz` | `qc.tsv`, `qc.tsv.json` |
| `run` | a contig | whole-genome GBZ, targets | every contig stage's outputs |
| `merge` | a genome | per-contig GBZs | one GBZ |
| `union` | a contig | two `packages.jsonl.gz` | the union's `packages.jsonl.gz` |
| `project` | a contig | the union's packages, its `msas.jsonl.gz` | `msas.jsonl.gz` for another graph ([union](#union-and-project)) |

Each stage prints a JSON summary on standard output (`prepare` prints its output paths) and writes a
manifest (`verify` only with `-o`), a JSON record of what it ran ([format](#outputmanifestjson)). One
contig, stage by stage (`run` does the same in one command):

```sh
pgrealign prepare --gbz genome.gbz --contig chr1 -o work
pgrealign regions --targets catalogs/chm13v2.0_tr.bed.gz --contig chr1 \
    --snarls work/chr1.snarls.pb --ref work/chr1.ref.tsv.gz -o work/regions.tsv.gz
pgrealign extract --graph work/chr1.gbz --regions work/regions.tsv.gz -o work/extract
pgrealign realign --packages work/extract/packages.jsonl.gz -o work/msas.jsonl.gz -j 8
pgrealign replace --graph work/chr1.gbz --extract work/extract --msas work/msas.jsonl.gz \
    -o work/chr1.realigned.gfa --gbz work/chr1.realigned.gbz
pgrealign verify --before work/chr1.gbz --after work/chr1.realigned.gbz -o work/verify.json
pgrealign qc --extract work/extract --msas work/msas.jsonl.gz -o work/qc.tsv -j 8
```

The catalogue's path is relative to the repository's root. Its contig names and coordinates are
those of CHM13 v2.0, so it suits a graph whose reference sample is CHM13 v2.0, with contigs named
chr1, chr2 and so on.

## The command line

Every stage is run as `pgrealign <stage> [options]`, and `pgrealign --version` prints the version.
Options that control the same thing in several stages:

| option | stages | controls |
|---|---|---|
| `--vg` | prepare, regions, extract, replace, verify, run, merge | the vg binary |
| `--abpoa` | realign, qc, run | the abPOA binary |
| `-t`, `--threads` | prepare, extract, replace, run | threads for vg |
| `-j`, `--jobs` | catalog, realign, qc, run | the most work done in parallel |
| `--max-memory` | realign, qc, run | the memory, in GB, that the regions worked on at once may need in all; it decides only how many run at once ([Memory](#memory)) |
| `--max-pair-memory` | qc, run | the largest pair `qc` compares, by its memory estimate ([qc](#qc)) |
| `--ref-sample` | prepare, extract, qc, run | the reference sample, CHM13 unless set |
| `--drop-sample`, `--keep-all-samples` | prepare, run | the samples `prepare` removes ([Dropped samples](#dropped-samples)) |
| `--scores` | realign, qc, run | abPOA's scores |
| `--max-cells` | realign, qc, run | the largest pairwise alignment run unbanded ([Limits](#limits), [qc](#qc)) |
| `--centre`, `--max-excess` | realign, run | the allele each star is built around: the reference allele, unless too much of the haplotypes' sequence lies beyond its length, or the medoid ([realign](#realign)) |
| `--max-work`, `--max-region-memory`, `--max-address-space` | realign, run | the largest region `realign` aligns, by its work, memory estimate and abPOA's reservation, each decided from its package ([Limits](#limits)) |
| `--timeout` | realign, run | `realign`'s time limit per region, none by default ([Limits](#limits)) |
| `--min-share` | realign, run | the least share of an allele's 15-mers with the centre for its region to be aligned ([Limits](#limits)) |
| `--max-internal-bp` | extract, run | the longest anchor-free walk `extract` records ([Statuses](#statuses)) |
| `--id-start` | replace, run | the first node id `replace` gives ([Node ids](#node-ids)) |

`--scores` takes abPOA's own scoring flags as one argument. Leaving it out uses the method's scores,
`star.SCORES` ([Scores](#scores-and-scores)). A flag left out of `--scores` takes abPOA's default, so
`--scores ''` gives abPOA's own defaults.

### catalog

`catalog` builds a catalogue of the tandem repeats in a reference genome, the usual list of targets.
No other stage depends on the targets being tandem repeats, so any BED of intervals can take the
catalogue's place. It reads one or both of:

- `--repeatmasker`: a RepeatMasker annotation as BED, with at least seven tab-separated columns.
  `catalog` uses the name (column 4) and the repeat class (column 7), and keeps the Simple_repeat,
  Satellite and Low_complexity records.
- `--reference`: reference FASTA files, plain or gzipped, given one after another after a single
  `--reference`, as in `--reference chr1.fa chr2.fa`. `catalog` runs trfind, its own tandem-repeat
  finder, on them, on N files at once with `-j N`, so with one file per contig the contigs run in
  parallel.

- `--trf`: a TRF binary, as a path or a name on `PATH`, to run on the `--reference` files as well, the
  way RepeatMasker 4.1 finds its Simple_repeat records, in place of a RepeatMasker annotation the
  reference lacks ([method](method.md#trf)). `-j N` runs N TRF processes at once.
- `--strf`: runs strf, pgrealign's reimplementation of the part of TRF that `--trf` uses
  (`pgrealign/c/strf.c`, [method](method.md#strf)), in TRF's place: on the same batches, with the same
  stages, `-j N` processes at once. Alone, `--strf` compiles strf on first use, as trfind is compiled;
  `--strf BINARY` names a strf binary built elsewhere from `pgrealign/c/strf.c`, as a path or a name on
  `PATH`. `--trf` and `--strf` cannot both be given.

Sources that hold no record on the selected contigs, such as a reference with no tandem repeat, give
a catalogue with a header and no entry; giving neither source is an error.

`--contigs` takes a comma-separated list of contig names, or `all`. Without it, `catalog` keeps a
fixed list, chr1 to chr22 and chrX, and drops every other contig, such as chrY, chrM, or the records
of a FASTA named 1, 2 and so on; such contigs need `--contigs` or `all`.

It writes the catalogue to `-o`, gzipped when the name ends in `.gz`
([format](#target-beds-and-the-catalogue)). Its manifest holds each input's MD5 as well as its size,
so that someone rebuilding the catalogue can confirm they have the same inputs. The repository's
`catalogs/chm13v2.0_tr.bed.gz` is the CHM13 v2.0 catalogue, and its manifest holds the command that
built it.

`catalog` holds each FASTA file's selected sequences in memory while trfind runs on them, and
trfind holds its own copy, for N files at once under `-j N`. It also writes an uncompressed copy of
each file's selected sequences to the system's temporary directory for trfind to read.

The summary holds the number of `records` kept (RepeatMasker records of the three kept classes, and
trfind's records, on the selected contigs), the number of `regions` (catalogue entries) written,
`records_by_source_class` (keys such as `rm/STR` and `trfind/VNTR`), `regions_by_class`, the
entries' total length in `bp`, and, when trfind ran, `trfind`: the hashes, compiler flags and
platform that name its build, or, with `--trfind`, `binary`, the path of the binary that ran, and,
with `--trf`, `trf`: the TRF binary's path and the parameters of each stage. With `--strf`, `strf`
holds the parameters of each stage too, along with strf's build, or with `--strf BINARY` its `binary`,
as `trfind` does.

#### Records, trfind and merging

[method.md](method.md#the-catalogue) gives the rules by which `catalog` gives each record a class,
`STR`, `VNTR`, `SAT` or `LC`, finds repeats with trfind, and merges records into entries.

- trfind (`pgrealign/c/trfind.c`) finds tandem repeats from the sequence alone. Each 10-mer proposes
  as a period the distance to its previous copy, when that lies between `--min-period` and
  `--max-period`, both included; seeds of about the same period chain into a repeat. It reports the
  smallest divisor of the chain's period that explains the repeat about as well, so a reported
  period can fall below `--min-period`.
- `catalog` compiles trfind, and strf with `--strf`, on first use into `$PGREALIGN_CACHE`, or else
  `~/.cache/pgrealign`, with `$CC`, or else the first of `cc`, `clang` and `gcc` on `PATH`. The
  cached binary's name holds the tool's name and a hash of its source, the compiler flags and the
  platform, so a changed source is compiled again.
- `--trfind` names a trfind binary built elsewhere from `pgrealign/c/trfind.c`, as a path or a name
  on `PATH`, for a system with no C compiler, such as a container, and `--strf BINARY` names a strf
  binary built from `pgrealign/c/strf.c` the same way. `catalog` then compiles neither, and cannot
  check that the binary was built from this version's source.
- `catalog` merges the records of one contig that overlap or lie within `--gap` bp of each other
  into one entry. Raise `--gap` to join repeats separated by short stretches into one entry, and so
  into one target.

### prepare

`prepare` cuts one contig out of a whole-genome GBZ and builds the two indexes `regions` needs. It
reads `--gbz` and `--contig`, and writes to the directory `-o`:

- `<contig>.gbz`: the contig's connected component, from `vg chunk --gbz --contig`, without the
  dropped samples. Node ids stay those of the whole-genome graph. `vg chunk` must write exactly one
  component, or `prepare` stops.
- `<contig>.snarls.pb`: the snarl decomposition from `vg snarls -T -P <ref>#0#<contig>`, trivial
  snarls included, in vg's binary format.
- `<contig>.ref.tsv.gz`: the reference path's nodes and their positions ([format](#contigreftsvgz)).
- `<contig>.manifest.json`.

The reference path is the reference sample's W line on the contig, or a P line named
`<ref>#0#<contig>` or `<ref>#<contig>`. `prepare` counts the W lines and P lines together, and stops
unless it finds exactly one, so it refuses a reference split into several walks on the contig. A
P-line reference serves `prepare` and `regions`, but `extract` records runs of W lines only, so its
packages name no reference allele, and `replace` stops at a P line that visits a removed node
([Known gaps](#known-gaps)).

`prepare` writes each output under a temporary name and renames it into place when it is complete.
It skips each output that exists, so running it again after a failure resumes. It does not check the
options an existing output was made with: after changing `--ref-sample` or the dropped samples,
delete the outputs first.

#### Dropped samples

By default `prepare` removes the sample `gref_CHM13`, the graph's reference path cover: a copy of
the reference path plus the pieces of other haplotypes that hang off it, as `vg paths -u` computes
it under the sample name `gref_<ref>`. The cover's short paths start and end inside regions, which
would give those regions fragments.

- `--drop-sample` (repeatable) names the samples to remove instead, and with `--keep-all-samples`
  `prepare` removes none. A graph whose reference sample is not CHM13 needs `--drop-sample` set to
  its own cover, `gref_<ref>`, when it has one. vg warns about a sample the graph lacks.
- `prepare` removes the samples with `vg gbwt -R`, and then sets the GBZ's reference samples to the
  original reference samples that were not dropped. It refuses to drop the reference sample, and
  stops if dropping would leave no reference sample.

The realigned graph has no cover either ([Known gaps](#known-gaps) gives the commands that add one).

#### A graph already cut

`--graph`, in place of `--gbz`, names a graph that is already the contig's: a GFA, a GBZ, or a graph
in another of vg's formats, such as the `<contig>.vg` that Minigraph-Cactus writes for each reference
chromosome. `prepare` then writes only `<contig>.snarls.pb`, `<contig>.ref.tsv.gz` and the manifest,
and the later stages read the graph itself, by its own name. It cuts nothing, so the graph may hold
the reference sample's paths on other contigs, and it drops no sample: it takes neither
`--drop-sample` nor `--keep-all-samples`, so a cover the graph holds must be removed beforehand.

### regions

`regions` places each target in the contig's chains and merges targets whose regions overlap. It
reads the targets in `--targets` that lie on `--contig` (a catalogue or BED,
[format](#target-beds-and-the-catalogue)), and `--snarls` and `--ref`, the files `prepare` wrote. It
writes `regions.tsv` to `-o`, gzipped when the name ends in `.gz` ([format](#regionstsv)).

It reads the snarls and the reference index only when the contig has targets; with none, it writes
the header alone. When the contig has targets but no chain has two boundary nodes that the reference
visits once, forward, `regions` stops with an error.

#### Anchoring and merging

`regions` pads each target by `--pad` bp on each side and takes its anchors from the chain, top-level
or nested, that encloses the padded target in the fewest reference bases. Taking the regions in
order of their left anchors, it merges one into the region before when its left anchor starts before
that region's right anchor starts, and anchors the merged targets again as one; regions that only
share a boundary node stay separate ([method.md](method.md#anchors) and
[method.md](method.md#merging) give the rules in full).
Raise `--pad` when realignment should be able to move variation that the graph places just outside a
target; a larger pad also merges more nearby targets into one region.

#### Size

A region whose anchors cannot lie within `--max-span` bp of each other, measured from L's start to
R's end, gets the status `too_long`, and later stages leave it as it is. A target that no chain
encloses within that distance keeps its own interval, has no anchors and takes no part in merging,
so it can overlap an `ok` region; merged targets that cannot be anchored as one keep the outermost
anchors of the regions merged. With a larger `--max-span`, `regions` admits longer regions, whose
alignment costs more memory and time in `realign`.

The summary holds the number of the contig's `targets`, of `chains` found (null with no targets), of
`regions` and of each `status`, the number of targets in `ok` regions, including targets alone in
their region (`merged_targets`), and the reference span of the `ok` regions in bp (`span_bp_ok`).

### extract

`extract` reads every walk of a contig's graph once, and records every run of every walk through
every region. It reads:

- `--graph`: the `<contig>.gbz` from `prepare`, or a GFA. `extract` looks up node sequences as it
  reads the W lines, so a GFA must have every S line before its first W line, as `vg convert -f`
  writes them, and unique walk keys.
- `--regions`: the `regions.tsv` from `regions`. Only regions that are `ok` there get a package.

In the directory `-o` it writes `packages.jsonl.gz`, one package per region
([format](#packagesjsonlgz)), `occurrences.tsv.gz`, one row per run ([format](#occurrencestsvgz)),
`interiors.tsv.gz`, each region's interior nodes ([format](#interiorstsvgz)), `summary.json`
([format](#summaryjson)), and the manifest `extract.manifest.json`.

`extract` holds every node's sequence in memory while it reads, so its memory grows with the
contig's total sequence. `--ref-sample` names the sample whose spanning run gives each package's
`ref` allele.

#### Statuses

`extract` gives each package the first of these statuses that applies, in this order. Only `ok`
packages go on to `realign` and `replace`.

- `complex`: a run is complex, runs visit an interior node both ways once turned to the region's
  forward orientation, as at an inversion ([method.md](method.md#orientation)), an anchor-free walk
  ([Regions and runs](#regions-and-runs)) uses the interiors of more than one region or runs both
  ways through one, or an anchor-free walk longer than `--max-internal-bp` uses the interior.
- `no_runs`: no walk has a run in the region.
- `invariant`: there are no fragments, and all spanning runs spell the same sequence (one allele).
- `ok`: every other region.

`extract` keeps an anchor-free walk whole, as a possible internal run, only when it spells at most
`--max-internal-bp` bases, and never records a longer one as a run
([method.md](method.md#walks-that-visit-no-anchor)). A larger `--max-internal-bp` leaves fewer
regions `complex` for this reason (the summary's `internal_too_long`), at the cost of memory.

### realign

`realign` builds the star alignment of each region with abPOA. It reads `--packages`, the
`packages.jsonl.gz` from `extract` (or from `union`), and writes `msas.jsonl.gz` to `-o`
([format](#msasjsonlgz)). Each star is built around the region's reference allele, where it has one
of positive weight and no more than `--max-excess` (`realign.MAX_EXCESS`, 0.1) of the haplotypes'
bases lie beyond its length, and around the medoid of its alleles otherwise; `--max-excess 0` takes
the reference allele whatever, and `--centre medoid` takes the medoid always
([method.md](method.md#the-centre)).

It aligns the packages whose status is `ok` and that `--select` selects; its summary lists the rest
of its shard ([below](#parallel-runs-and-resuming)) as `skipped`, by status or as `not_selected`.
With `--select any`, the default, it realigns every region whose alleles are not all identical.
With `--select length`, it realigns only regions whose alleles differ in length, which saves the
time of regions whose alleles all have one length, such as those that differ only by substitutions.
It always realigns a region that has fragments.

#### Limits

`replace` leaves a region that `realign` cannot align as it is, and the region's record names the
reason in its status ([statuses](#msasjsonlgz)). Whether a region is aligned is decided from its
package alone, before anything is aligned, so it does not depend on the machine, on the memory
budget or on how long anything takes: the same packages give the same records anywhere. Each limit
on a region's size below leaves a region that exceeds it `too_big`, with a message that says which.

- `--max-cells` bounds the cells of one unbanded pairwise alignment, and so its memory: abPOA
  touches 20 bytes per cell, so the default allows about 8 GB. An alignment to the centre whose
  matrix, the product of the two sequences' lengths, would exceed the cap is made by
  `star.long_pair` instead ([method.md](method.md#pairs-too-large-to-align-unbanded)): with the
  widest adaptive band that keeps within the cap (`star.band_width`), and seeded, the one that
  scores higher kept, and a kept seeded one refined in windows of at most the cap, each aligned
  unbanded. Every other alignment to the centre is unbanded. The record's `banded_pairs` counts
  the pairs over the cap, `seeded_pairs` those whose seeded alignment was kept, and
  `refined_windows` the windows refinement improved.
- Such a pair's memory is about that of an unbanded alignment of the cap's cells, and more where the
  two lengths differ by more than a band within the cap allows, since its band then spans their
  difference: abPOA reserves 20 bytes of address space per cell of the whole matrix, banded or not,
  and a band touches at least one page of each of its five score matrices per base of the centre,
  and more where the two lengths differ ([Memory](#memory)).
- A region is `too_big` when even abPOA's own band would exceed the cap aligning its longest allele
  or fragment to one of the same length: when `realign.banded_cells(n)`, min(n², 2n(10 + ⌊n / 100⌋))
  for a sequence of n bases, is larger than `--max-cells`, so the cap sets the longest sequence a
  region can hold. Raise it, with the memory to match, to align such regions.
- A region whose alleles hold more than `realign.MAX_TOTAL_BP` (40,000,000) bases in all is also
  `too_big`, because choosing the centre holds every allele's k-mer counts in memory.
- `--max-work` bounds the region's work, about the cells of dynamic-programming matrices that its
  alignments fill (`realign.work`): for each allele and fragment aligned to the centre, the product
  of the two lengths, or, where that exceeds `--max-cells`, `realign.long_pair_cells(n, m)`: the
  centre's length times the band, which spans the two lengths' difference and `star.band_width`
  cells either side (`realign.banded_cells(n, m, w)`), about `--max-cells` in all where the lengths
  allow a band wider than abPOA's own. Seeding, which takes a small fraction of the band's time, and
  refinement, which runs only where the seeded alignment is kept, are not counted. And for the
  insertion slots, the bases by which each member is longer than the centre, its excess, taken to
  fall in one slot, each member's excess aligned, banded, to the longest. Around a medoid, which is
  not known before the region runs, the longest allele stands in for the centre, and the slots are
  not counted. Run time follows the work
  ([method.md](method.md#limits-and-regions-that-keep-their-old-subgraph)). The default is
  `realign.MAX_WORK`, 2 × 10^11 cells, and 0 sets no limit.
- `--max-region-memory` bounds the region's memory estimate, in GB ([Memory](#memory)). The default
  is `realign.MAX_REGION_MEMORY`, 64, and 0 sets no limit.
- `--max-address-space` bounds the address space, in GB, that abPOA reserves for one of the
  region's alignments (`realign.reservation`): the alignments to the centre, and the slots' as the
  work takes them, of which the second longest excess aligned to the longest reserves the most;
  around a medoid, the longest allele stands in for the centre, and the longest sequence against
  one of the same length bounds the slots ([Memory](#memory)). Linux refuses a reservation larger
  than the machine's memory and swap, so without this limit a region that needs one would be
  aligned on a large machine and fail on a small one. The default is `realign.MAX_ADDRESS_SPACE`,
  128, which admits alignments of two sequences of up to 82,896 bp, or of 62 kb to 110 kb, and 0
  sets no limit.
- `--timeout` sets a time limit, in seconds, for one region; a region stopped by it is `timeout`.
  There is none unless one is given, since a time limit makes a region's outcome depend on the
  machine and on what else runs there. `realign` checks the limit while choosing a medoid, once per
  allele, and at each abPOA call, which it does not start past the limit and stops when the limit is
  reached. The work before the first check (counting k-mers) and after the last abPOA call is not
  limited, so a region can run past the limit and still be `ok`.
- `--min-share` refuses a region, as `not_homologous`, once its centre is known and before anything
  is aligned, when an allele of positive weight shares less than that fraction of its distinct
  15-mers with the centre, over those of whichever of the two has fewer. An allele shorter than
  `star.HOMOLOGY_BP` (100) bases is not judged, nor is any allele when the centre is that short, nor
  any fragment. Such an allele is unrelated to the centre or reversed against it, and aligning it
  forward would thread it through the centre with chance matches
  ([method.md](method.md#homology)). The default is 0.5, and 0 aligns every region.

`--scores` sets abPOA's scores for every alignment ([Scores](#scores-and-scores)).

`--plan` aligns nothing and writes nothing, but prints what `realign` would do with the packages
under the same options: the regions it would align (`admitted`), those it would leave `too_big`,
counted by the limit each exceeds (`too_big`), the rest (`skipped`), the admitted regions' `work` in
all, and `need_gb`, the largest admitted region's memory estimate in GB, which is the least
`--max-memory` that aligns them, with that region's id (`largest`); and `reservation_gb`, the largest
reservation abPOA would make for them, against `granted_gb`, the largest this machine grants (null for
no limit), which must be at least as large ([Memory](#memory)).

#### Memory

Before it starts, `realign` estimates from each region's package the most memory the region can
take. It starts regions in package order while the estimates of the regions running, the next
one's included, add up to no more than `--max-memory` GB, and while fewer than `-j` run. The
default is half the machine's physical memory. The budget decides only how many regions run at
once, never which are aligned: a region whose estimate exceeds `--max-region-memory` is `too_big`
under any budget, and if one within it needs more than `--max-memory` on its own, `realign` stops
before it aligns anything, with an error that names the region and its estimate. Give `realign` that
much (`--plan`'s `need_gb`), or lower `--max-region-memory` below it, which leaves such regions
`too_big`. Since the budget bounds the total, `-j` can be the number of cores.

abPOA aligns a sequence of m bases to one of n bases given before it (the centre, when `realign`
aligns an allele to the centre) in five matrices of 4-byte scores, with one row per base of the
first sequence. It reserves them whole, and only the pages it touches take memory. Unbanded, it
touches every cell: 20nm bytes. Banded, with a margin of w cells either side, by default abPOA's own,
10 + ⌊m / 100⌋, each row of each matrix touches a run of about |n − m| + 2w cells that starts
anywhere in a page, so the alignment touches about 5n(P + 4(|n − m| + 2w)) bytes, P being the page
size, and never more than unbanded (`memory.abpoa_bytes`). Measured as macOS's peak memory
footprint, which counts memory compressed or swapped out as well as resident (abPOA 1.5.7, 16 KB
pages, the default `--scores`), against this formula, in GB of 2^30 bytes:

| alignment | footprint | formula |
|---|---|---|
| unbanded, two sequences of 7.5, 15, 30 kb | 1.05, 4.20, 16.8 | 1.05, 4.19, 16.8 |
| banded, two near-equal sequences of 7.5, 15, 30, 58.7 kb | 0.60, 1.24, 2.65, 5.85 | 0.60, 1.23, 2.64, 5.80 |
| banded, a 58.7 kb centre with a 40.1 kb or a 103.2 kb allele | 19.2, 47.1 | 25.7, 55.4 |

The last two are only 10.1 and 11.4 GB resident: an alignment's resident size understates what it
needs. With 4 KB pages, as on most Linux systems, the per-page part is a quarter as large. There,
the band `star.long_pair` gives a pair over the default `--max-cells` peaked at 7.6, 6.9 and 8.1 GB
resident for alleles of 27.8, 37.2 and 34.4 kb aligned to centres of 27.7, 27.7 and 32.3 kb (Linux,
abPOA 1.5.7), against 8.9, 8.9 and 9.0 GB by the formula, and abPOA's own band, the margin it gives
where the lengths differ too much for a wider one, 10.2 GB for a 48.8 kb allele to a 27.7 kb
centre, against 11.9 GB.

abPOA reserves its five matrices whole, 20 bytes for each of (n + 1)(m + 1) cells, rounded up to a
power of two, banded or not (`memory.abpoa_reservation`): two sequences of 62,146 bp reserve 2^37
bytes, 128 GB, though banded they touch about 6 GB. macOS grants such a reservation; Linux grants
one only up to the machine's memory and swap, and abPOA then exits with `posix_memalign fail`.
`--max-address-space` leaves such regions `too_big` on every machine alike ([Limits](#limits)). A
machine must also grant the reservations of the regions within it: `realign` stops before it aligns
anything, naming the largest, if one would need more than this machine's memory and swap (its commit
limit under Linux's strict overcommit policy, and nothing under its 'always' policy or elsewhere). A
region the machine fails all the same, by refusing abPOA's reservation, killing abPOA (as the OOM killer
does) or running the worker out of memory, is not recorded: `realign` aligns the others, then stops
with an error that names it, and run again where it can be aligned, or with lower limits, it resumes
from the regions it finished ([below](#parallel-runs-and-resuming)).

A region's estimate, `realign.peak_bytes`, is the sum of:
- 100 MB for the worker and abPOA's own tables (`memory.BASE`);
- the k-mer counts that choose the centre, 80 bytes per base of the alleles, and 50 more per base of
  the longest allele for comparing two alleles or sketching one (`realign.KMER_BYTES` and
  `realign.COMPARE_BYTES`);
- abPOA's matrices for the worst alignment the region can need, the largest value of the formula:
  around the reference allele, over its alignments to the reference (`memory.pair_bytes`), and over
  the slots' as the work takes them, each member's excess aligned, banded, to the longest excess;
  around a medoid, which is not known before the region runs, over every two of its alleles and
  fragments, in either order, and the longest sequence banded against one of the same length, the
  bound taken for the slots. A pair whose product of lengths exceeds `--max-cells` is taken at the
  larger of its band's, with the margin `star.band_width`, and refinement's windows, 20 bytes for
  each of `--max-cells` cells. Seeding's alignment is not counted: it took at most 1.3 GB resident
  on the pairs of up to 71 kb measured (Linux, abPOA 1.5.7), aligning window by window under abPOA's
  own band, against about 8 GB for the band's.

The page size P in the formula is `realign.ESTIMATE_PAGE`, 16 KB, whatever the machine's, so that the
estimate, and with it which regions are `too_big`, is the same on every machine; with Linux's 4 KB
pages it overstates a banded alignment. A region that the checks on cells and on total bases make
`too_big`, or that has no allele, is estimated at 100 MB, since it is refused before any alignment. The estimate is an upper bound in
practice, not a guarantee: a band can reach further than |n − m| from the diagonal. It leaves out
the parent process, which holds every package to align and every finished record.
[whole-genome.md](whole-genome.md#memory) gives the estimates of two contigs' regions.

#### Parallel runs and resuming

With `-j N`, `realign` aligns up to N regions at once, each in its own process, as many as
`--max-memory` allows ([Memory](#memory)). With `--shard i/n`, it aligns only the packages at
positions i, i + n, i + 2n, ... of the packages file, counting from 0, one shard of n. Run the n
shards as separate jobs, each with its own `-o` and, on one machine, its own share of the memory in
`--max-memory`. `replace` takes all their outputs after one `--msas`; `qc` and `project` each read
one file, so concatenate the outputs with `cat` for them, since concatenated gzip files read as one:

```sh
pgrealign realign --packages work/extract/packages.jsonl.gz --shard 0/2 -o work/msas.0.jsonl.gz
pgrealign realign --packages work/extract/packages.jsonl.gz --shard 1/2 -o work/msas.1.jsonl.gz
pgrealign replace --graph work/chr1.gbz --extract work/extract \
    --msas work/msas.0.jsonl.gz work/msas.1.jsonl.gz -o work/chr1.realigned.gfa
cat work/msas.0.jsonl.gz work/msas.1.jsonl.gz > work/msas.jsonl.gz
pgrealign qc --extract work/extract --msas work/msas.jsonl.gz -o work/qc.tsv
```

As each region finishes, `realign` appends its record to `<out>.partial`
([format](#msasjsonlgzpartial)). A run that finds that file keeps the records in it that read back
whole and does not align those regions again, whatever their status or the options that made them
([Known gaps](#known-gaps)). Once every region is done, it writes `<out>` in package order and
deletes the partial file. To retry failed regions with larger limits, delete `<out>.partial` if one
is left and run `realign` again.

The summary holds the number of regions `realigned`, by status, and `skipped`, by reason, the total
`aligner_seconds`, the number of regions `resumed` from a partial file, and `peak_estimate_gb`, the
largest sum of the estimates of the regions running at once ([Memory](#memory)). The limits on a
region's size apply afresh to a resumed run: a region the partial file holds is `too_big` if the
limits given now make it so, whatever its record, and aligned if its record is `too_big` but they no
longer make it so, so that the run ends as a fresh one would.

### replace

`replace` writes a contig's graph with the interior of each realigned region replaced by the
region's induced graph, the graph of its MSA: one node per distinct base in each column, with nodes
merged where nothing branches, and one walk per row. It reads:

- `--graph`: the graph `extract` read, since the step indices in `occurrences.tsv.gz` refer to its
  walks. A GFA must have every S and L line before its first W line, as `vg convert -f` writes them.
- `--extract`: `extract`'s output directory, with its `packages.jsonl.gz`, `occurrences.tsv.gz` and
  `interiors.tsv.gz`.
- `--msas`: one or more `msas.jsonl.gz` files, given after one `--msas`, such as one per shard. A
  region may appear in only one.

It writes a GFA to `-o`, under a temporary name that it renames when the GFA is complete. With
`--gbz`, it also builds that GBZ with `vg gbwt -G` and then deletes the GFA unless `--keep-gfa` is
given; the manifest goes next to the GBZ in that case, and next to the GFA otherwise. It replaces a
region only when both its package and its MSA are `ok`. The summary is described
[below](#replaces-summary).

#### What changes in a replaced region

`replace` removes the interior's nodes and their edges, adds the induced graph
([induce](#induceinduce-and-graph)) joined to the anchors, and rewrites every run recorded in
`occurrences.tsv.gz` to its member's walk in the induced graph
([method.md](method.md#what-changes-and-how-it-is-checked)). A node between the anchors that no
recorded run visits stays, without its edges to removed nodes. H lines, the header lines that hold
tags about the whole graph, lose vg's graph-name tags (NM, SG and TL), so that `vg gbwt -G` does not
record the output as a subgraph of the input. P lines are renumbered under the `dense` numbering,
and lines of any other type copied unchanged.

#### Node ids

`--id-mode` sets how nodes are numbered:

- `dense`, the default, renumbers every node. The nodes that stay keep their order, each region's new
  nodes come right after its left anchor, and ids run consecutively from `--id-start`, or else from
  the contig's smallest id. Every S, L, W and P line is renumbered, inside the regions or not. New
  nodes stay next to their neighbours, which suits tools that read the graph in windows of node ids,
  such as vg call.
- `reuse` keeps the old nodes' ids. A region's new nodes take its freed ids, then ids from
  `--id-start`, which this mode requires and which must start a range that is free in the whole
  genome. New nodes can then land far from their regions. Choose it when files made against the input
  graph, such as read alignments, must keep pointing at the same nodes outside the regions.

`vg stats -r` on the whole-genome GBZ prints its smallest and largest node ids, so an `--id-start`
above the largest starts a range that is free in the whole genome. Under `reuse`, the summary's
`ids_above_start` says how many ids the contig took from that range, and `merge` refuses a contig
that took any if other contigs' ids lie between its old ids and the new ones
([Known gaps](#known-gaps)). Under `dense`, the summary's `id_range` gives the first and last id of
the output.

#### Checks

`replace` stops, and leaves no output, at any of these:

- `--id-mode reuse` without `--id-start`, or a region in two `--msas` files;
- an anchor of one replaced region in the interior of another;
- an MSA whose induced graph does not spell its rows, or a member with no row in its MSA;
- under `dense`, a left anchor that is not a node of the graph; under `reuse`, two new nodes given
  one id, or a new node given the id of a node of the graph;
- an S or L line after a W line;
- a walk whose recorded runs overlap, or that visits a removed node outside a recorded run;
- a run whose old steps, or whose member's new walk, do not spell the member's sequence in the
  package, or that visits a node no replaced interior holds or whose S line comes after the walk;
- a P line that visits a removed node ([Known gaps](#known-gaps)).

Most of these mean that the `extract` or `realign` outputs do not belong to this graph, or are
wrong.

### verify

`verify` checks that every path of the realigned graph, in vg's sense, spells what it spelled
before. Under the `dense` numbering node ids change throughout the contig, so it compares the graphs
by what their paths spell rather than step by step. It reads `--before`, the graph given to
`replace --graph`, and `--after`, the realigned GBZ; either can be a GBZ or a GFA. It requires the
path metadata from `vg paths -M` (name, sense, sample, haplotype, locus, phase block and subrange)
to be identical, and every path's sequence, from `vg paths -F` and upper-cased, to have the same
length and MD5 in both.

It prints its report ([format](#verifys-report)), with the node count, edge count and sequence
length of both graphs, and writes it to `-o` as well when given, with a manifest. It exits with
status 1 when the check fails. It reads every base of every path once per graph.

### qc

`qc` measures how far each region's MSA falls from optimal pairwise alignments
([method.md](method.md#qc)). An MSA implies an alignment of every pair of its rows: the columns
where either row has a base. `qc` compares the implied alignments of a sample of allele pairs with
each pair's optimal alignment, which it computes with abPOA as a global, unbanded alignment under
`--scores`, left-normalised as `realign` does ([method.md](method.md#left-normalisation)). Give it
the scores `realign` used.

It reads `--extract`, `extract`'s directory, and one `--msas` file, and writes a TSV to `-o` and its
summary to `<out>.json` ([format](#qctsv-and-qctsvjson)). It measures each region whose MSA is `ok`
with more than one row and whose package has at least two alleles.

- `--graph` names the graph `extract` read, a GBZ or a GFA. `qc` then also measures that graph's own
  alignment of the same pairs, against the same optimal alignments
  ([method.md](method.md#the-graphs-own-alignment)), which takes one more pass over the graph, by
  `vg convert` for a GBZ (`--vg`, `-t`). `run` gives it the contig's graph.

- `--pairs` sets how many random pairs of alleles, the reference allele excluded, it measures per
  region, and `--refs` how many of the heaviest other alleles it compares with the reference allele.
  More pairs measure a region more closely, at the cost of more abPOA runs.
- A pair whose optimal alignment would exceed `--max-cells` cells, the product of the two alleles'
  lengths, is compared with abPOA's banded alignment instead; the `banded` column counts such pairs.
- `-j` and `--max-memory` work as in `realign` ([Memory](#memory)): `qc` starts regions in order
  while the estimates of those running add up to no more than `--max-memory` GB, and while fewer
  than `-j` run. A pair's estimate (`qc.pair_bytes`) is 100 MB and the [formula](#memory) for its
  alignment, unbanded or banded by `--max-cells`, plus, with `--graph`, 120 bytes per base of the
  region's alleles for the graph positions `qc` holds while it aligns (`qc.GRAPH_BYTES`); a
  region's (`qc.peak_bytes`) is that of its largest compared pair, with 16 KB pages whatever the
  machine's. A pair whose estimate exceeds `--max-pair-memory` GB (`qc.MAX_PAIR_MEMORY`, 64; 0 sets
  no limit) is not compared, on any machine, and counts in the `skipped` column and in no other. The
  budget decides only how many regions run at once: if a region needs more than `--max-memory` on its
  own, `qc` stops before it measures anything and names it.
- With `--sample N`, `qc` draws N of the `ok` MSAs with more than one row, and measures those whose
  package has at least two alleles, so it can measure fewer than N.
- `--ref-sample`: for a package without a `ref` field, `qc` reads the reference allele from this
  sample's runs in `occurrences.tsv.gz`. A package whose `ref` is null has the field, so `qc` makes
  no reference comparisons for it.

`qc` draws each region's pairs with a seed computed from the region's id, and draws `--sample` with
a fixed seed, so a rerun measures the same pairs.

To try other scores, run both `realign` and `qc` with them and compare the summaries;
[Reading the results](#reading-the-results) says what the values mean.

### run

`run` runs every contig stage in order for one contig of a whole-genome GBZ, skipping the stages
whose output exists. It reads `--gbz`, `--contig` and `--targets`, and writes under `-o`, in
`OUT/<contig>/`:

| stage | output | manifest | skipped when |
|---|---|---|---|
| prepare | `C.gbz`, `C.snarls.pb`, `C.ref.tsv.gz` | `C.gbz.manifest.json` | each output exists |
| regions | `regions.tsv.gz` | `regions.tsv.gz.manifest.json` | it exists |
| extract | `extract/` | `extract/extract.manifest.json` | `extract/summary.json` exists |
| realign | `msas.jsonl.gz` | `msas.jsonl.gz.manifest.json` | it exists |
| replace | `C.realigned.gbz`, `replace.json` | `C.realigned.gbz.manifest.json` | the GBZ exists |
| verify | `verify.json` | `verify.json.manifest.json` | it exists |
| qc | `qc.tsv`, `qc.tsv.json` | `qc.tsv.manifest.json` | `qc.tsv` exists, or `--no-qc` is given |

C stands for the contig's name. An interrupted `realign` resumes from `msas.jsonl.gz.partial`, and
`replace` deletes its GFA, `C.realigned.gfa`, once the GBZ is built. `replace.json` holds
`replace`'s summary. `run` decides what to skip only by whether files exist
([Known gaps](#known-gaps)). Each file the last column names is written after the stage's other
files, under a temporary name renamed into place once complete, so a `run` killed partway through
a stage redoes that stage; the manifest comes last. A `verify.json` that records a failure stops
`run` again on every later run, until it and the realigned GBZ are deleted.

`run` passes `-j`, `--max-memory`, `--scores` and `--max-cells` to `realign` and `qc`, which run one
after the other, each within the whole of `--max-memory`; `-t` to vg; and `--vg`, `--abpoa`,
`--ref-sample`, `--drop-sample`, `--keep-all-samples`, `--timeout`, `--max-internal-bp` and
`--id-start` to the stages that take them. It uses the defaults of the stage options it lacks:
`regions`' `--pad` and `--max-span`, `realign`'s `--select` with one shard, `replace`'s `dense`
numbering, and `qc`'s `--pairs` and `--refs`, without `--sample`, so `qc` measures every region. It
gives `qc` the prepared `C.gbz` as `--graph`. To change those, run the stages one by one.

`run` stops with an error when `verify` fails, after writing `verify.json`. It prints progress lines
on standard error, and at the end a JSON object with `id_range`, the smallest and largest node id of
the realigned GBZ (`vg stats -r`), and an entry for each stage it ran: `prepare`'s output paths,
always; `ok`, `paths`, `bases` and `metadata_identical` from `verify`; and every other stage's
summary.

### merge

`merge` joins per-contig GBZs into one whole-genome GBZ without renumbering any node. Its arguments,
before or after `-o`, are GBZs whose node ids do not overlap: each realigned contig, and every part
of the genome that was not realigned, such as the contigs without targets and the components with no
reference contig. [whole-genome.md](whole-genome.md) gives the commands that cut out those parts,
with the same dropped samples as the realigned contigs. `merge` writes the merged GBZ to `-o`.

It stops, naming two inputs, when their node-id ranges (`vg stats -r`) overlap ([Known
gaps](#known-gaps) says how to move a contig). It sets the reference samples to every sample with a
reference path in any input, and checks that the output's node count, edge count and length are the
sums of the inputs', and its path names exactly the union of theirs ([method.md](method.md#merge)).
An output that fails the check is renamed to `<out>.failed`, and `merge` stops with an error.

The summary holds the number of `inputs`, their `id_ranges`, the output's `nodes`, `edges`, `length`
and `paths`, and its `reference_samples`.

### union and project

`union` and `project` realign regions on the haplotypes of one graph, the panel, and apply the
result to another graph, the recipient
([method.md](method.md#realigning-on-a-panel-union-and-project)). The recipient's walks must use the
panel graph's nodes under the same ids, as they do in the panel graph with some samples removed.
(`union`'s and `project`'s help, the fields `target_map` and `target_status`, and `union`'s summary
keys call the recipient the target.)
For one contig, run `prepare` on both whole-genome GBZs, which keeps their node ids, and `regions`
on the panel's output alone; `extract` both contig graphs with that one `regions.tsv`; `union` the
two packages files, `realign` the union, and `project` its MSAs onto the recipient; and `replace`
the recipient's contig graph with its own `extract` directory and the projected MSAs.

#### union

`union` reads `--base`, the panel's `packages.jsonl.gz`, and `--add`, the recipient's, and writes the
union's packages to `-o`, in the panel's order. It keeps only the regions in both files. A region's
two packages must have the same anchors, or `union` stops with an error. Both extracts read one
`regions.tsv`, so this check catches only packages made from different regions files, not a
recipient whose node ids differ from the panel's.

- A recipient allele with the same sequence as a panel allele is that allele, and the panel's weight
  stands. A recipient allele the panel lacks joins with weight 0, so the centre is chosen from the
  panel's alleles alone. Recipient fragments join in the same way, matched by side and sequence. A
  joining member's id is `t` followed by its id in the recipient.
- The union's status is the panel's, except that a region `complex` in the recipient is `complex`,
  because `replace` cannot rewrite the recipient's runs there.
- Each union package adds `target_map`, which maps each of the recipient's member ids to the
  union's, and `target_status`, the recipient's own status.

The summary holds the number of `regions` written, of regions in only one file (`base_only`,
`add_only`), and of members only the recipient has (`target_only_alleles`,
`target_only_fragments`).

#### project

`project` reads `--union`, the union's `packages.jsonl.gz`, and `--msas`, `realign`'s output for
it, and writes `msas.jsonl.gz` for the recipient's `replace` to `-o`. For each `ok` MSA of a region
in the union, it keeps the rows of the recipient's members, under the recipient's member ids, and
drops the columns that are gaps in all of them. The recipient's `replace` then rewrites the regions
that are `ok` in the recipient's own `extract` and have a projected MSA. The summary holds the
number of MSAs `projected` and `skipped`.

## File formats

The TSV files that pgrealign writes are tab-separated, and each starts with a line of column names
that begins with `#`, except `qc.tsv`. A file whose name ends in `.gz` is gzip-compressed.

### Target BEDs and the catalogue

`regions` and `run` read their targets with `catalog.read_catalog`, which accepts two kinds of file.
A catalogue is a file in which a line before the first data line is exactly this header:

```
#chrom	start	end	id	class	period	motif	n_members	member_classes	sources
```

Each line is one catalogue entry, made by merging repeat records
([catalog](#records-trfind-and-merging)). In the column names, a member is one of those records.

| column | meaning |
|---|---|
| `chrom` | contig name |
| `start`, `end` | the entry, 0-based and half-open |
| `id` | `TR_<chrom>_<start>`, unique in the file |
| `class` | `STR`, `VNTR`, `SAT` or `LC` ([method.md](method.md#the-catalogue)) |
| `period` | the motif's length; 0 when the motif is not a `(MOTIF)n` repeat, as for most satellites |
| `motif` | the motif; for a RepeatMasker record whose name is not `(MOTIF)n`, its name |
| `n_members` | the number of records merged into the entry |
| `member_classes` | the merged records by class, such as `STR:2,VNTR:1` |
| `sources` | the merged records by source, `rm` (RepeatMasker), `trfind`, `trf` or `strf`, such as `rm:2,trfind:1` |

`catalog` writes the contigs `chr1`, `chr2`, ... in numeric order, then other contigs in sorted
order, and each contig's entries by start. A repeated id is an error.

Any other file is read as BED: `chrom`, `start`, `end` and an optional name in column 4, with
further columns (score, strand, BED12 blocks) ignored, and blank, `#`, `track` and `browser` lines
skipped. A BED target has class `.`, period 0 and motif `.`. Every line must have
0 <= start <= end; a line that does not parse stops the reader, which names the line. A target's
id is:

- its name, when no other line has that name;
- otherwise `<name>_<chrom>_<start>`, or `R_<chrom>_<start>` for a line without a name (none, empty
  or `.`);
- with `_2`, `_3`, ... appended, in file order, where that would repeat a name kept by the first
  rule or an id already given, so `chr1 100 200 A` becomes `A_chr1_100_2` when another line also
  has the name `A` and a third holds the name `A_chr1_100`.

### `<contig>.ref.tsv.gz`

The reference path's nodes, one row per step, in the path's order:

| column | meaning |
|---|---|
| `node` | node id |
| `orient` | `+` when the reference visits the node forward, `-` in reverse |
| `start`, `end` | the node's interval on the contig, 0-based and half-open |

Positions start at the reference walk's start offset, or at 0 for a P line. [prepare](#prepare)
says which path it takes.

### `regions.tsv`

One row per region, in order of `ref_start`, then `ref_end`:

| column | meaning |
|---|---|
| `id` | the id of the region's first target |
| `contig` | the contig |
| `L`, `R` | the anchors' node ids, both forward on the reference; `.` for a `too_long` target that no chain encloses |
| `ref_start`, `ref_end` | the interior on the reference, from L's end to R's start; for `.` anchors, the target's own interval |
| `span` | the reference length from L's start to R's end; for `.` anchors, the target's length |
| `targets` | the merged targets' ids, comma-separated, by start |
| `classes` | their distinct classes, sorted and comma-separated (`.` for BED targets) |
| `status` | `ok` or `too_long` |
| `note` | for `ok`, `merged N targets` when N > 1 and otherwise empty; for `too_long`, why |

### `packages.jsonl.gz`

One JSON object per line, one per region that is `ok` in `regions.tsv`, in the same order. An
example, with spaces and line breaks added for reading, and its sequences shortened:

```json
{"id": "TR_chr1_11635", "contig": "chr1", "L": 1041, "R": 1046,
 "ref_start": 11435, "ref_end": 11973, "targets": "TR_chr1_11635", "classes": "VNTR",
 "status": "ok", "n_runs": {"spanning": 3, "suffix": 1}, "complex": [], "ref": "a1",
 "alleles": [{"id": "a1", "weight": 2, "seq": "ACGTTGCA"}, {"id": "a2", "weight": 1, "seq": ""}],
 "fragments": [{"id": "f1", "weight": 1, "side": "left", "seq": "ACG"}]}
```

| field | type | meaning |
|---|---|---|
| `id`, `contig`, `L`, `R`, `ref_start`, `ref_end`, `targets`, `classes` | as in `regions.tsv` | the region's row (`L` and `R` are integers) |
| `status` | string | `ok`, `invariant`, `no_runs` or `complex` ([statuses](#statuses)) |
| `n_runs` | object | the number of runs of each type that has any: `spanning`, `prefix`, `suffix`, `internal` |
| `complex` | array of strings | why the region is complex, each as `<walk key>: <reason>`, at most 20 |
| `ref` | string or null | the id of the allele a spanning run of the reference sample spells; if there are several such runs, the last in the order `extract` reads the walks; null if there is none |
| `alleles` | array | one object per allele: `id`, `weight` (its number of runs) and `seq` |
| `fragments` | array | one object per fragment: `id`, `weight` (its number of runs), `side` (`left`, `right` or `internal`) and `seq` |

Sequences are spelled from the graph's node sequences, in the region's forward orientation apart
from internal fragments ([Known gaps](#known-gaps)); an allele's may be empty. `union` adds two
fields, `target_map` (an object mapping each recipient member id to a union member id) and
`target_status` (a string). A program writing packages for the later stages needs only the fields
they read:

| stage | fields |
|---|---|
| `realign` | `id`, `status`, `alleles` (`id`, `seq`, `weight`), `fragments` (`id`, `seq`, `side`) |
| `replace` | `id`, `status`, `L`, `R`, and the members' `id`, `seq` and fragments' `side`; the `--extract` directory must also hold `occurrences.tsv.gz` and `interiors.tsv.gz` |
| `qc` | `id`, `alleles` (`id`, `seq`, `weight`), and `ref` if known |

### `occurrences.tsv.gz`

One row per run, apart from complex runs, which are not recorded:

| column | meaning |
|---|---|
| `walk` | the walk key, `sample#haplotype#contig#start` |
| `region` | the region's id |
| `type` | `spanning`, `prefix`, `suffix` or `internal` |
| `orient` | for a spanning, prefix or suffix run, `+` when the walk runs in the region's forward orientation and `-` otherwise; `+` for every internal run |
| `first`, `end` | the run's steps in the W line's walk, as 0-based indices [first, end); an internal run's are 0 and the walk's step count |
| `member` | the id of the allele or fragment the run spells |

Rows come in the order `extract` reads the walks, with internal runs last. `replace` rewrites the
steps that each row names.

### `interiors.tsv.gz`

One row per package, in package order, with two columns: `region`, and `nodes`, the ids of the nodes
the region's runs visit, sorted and comma-separated (empty when there are none). These are the nodes
`replace` removes.

### `summary.json`

`extract`'s summary, also printed and stored in its manifest:

| key | meaning |
|---|---|
| `walks` | the number of W lines read |
| `regions` | the number of regions given a package |
| `status` | the number of packages of each status |
| `internal_candidates` | the number of anchor-free walks held whole, as possible internal fragments |
| `internal_too_long` | the number of regions made `complex` by an anchor-free walk longer than `--max-internal-bp` |
| `both_ways` | the number of regions made `complex` because runs visit an interior node both ways |
| `internal_flipped` | the number of anchor-free walks that run against their region, and so were reversed and flipped |

### `msas.jsonl.gz`

One JSON object per region that `realign` attempted, failures included, in package order; the
regions it skipped have none. `realign` writes each record as its own gzip member, and
`realign.write_msas` and `project` write one gzip stream. Every reader of these files takes any gzip
file, including gzip files concatenated with `cat`.

| field | present | meaning |
|---|---|---|
| `id` | always | the region's id |
| `method` | always | `medoid_star` |
| `status` | always | see below |
| `seconds` | always | the time spent on the region |
| `centre` | `ok` | the centre allele's id |
| `rows` | `ok` | `{member id: row}`, one row per allele and fragment, all the same length |
| `columns` | `ok` | the rows' length |
| `slot_alignments` | `ok` | the number of insertion slots that needed an alignment: an insertion slot is a place before, between or after the centre's bases where other rows insert bases, and the distinct insertions of one slot, when there are two or more, are aligned to each other |
| `fragments_exact` | `ok` | the number of fragments placed in an allele's row |
| `fragments_aligned` | `ok` | the number of unplaced fragments, aligned to the centre |
| `banded_pairs` | `ok` | the number of alignments to the centre that would exceed `--max-cells` unbanded, aligned instead with a band and seeded ([limits](#limits)) |
| `seeded_pairs` | `ok` | the number of those for which the seeded alignment scored higher and was kept |
| `refined_windows` | `ok` | the number of windows in which refinement improved a kept seeded alignment |
| `message` | not `ok` | why the region was not aligned |
| `projected_from` | after `project` | the number of rows in the union's MSA |

| status | meaning |
|---|---|
| `ok` | the MSA was built |
| `too_big` | even a banded alignment of the longest allele or fragment would exceed `--max-cells`, or the alleles hold more than `realign.MAX_TOTAL_BP` bases in all, or the region's work, memory estimate or abPOA's reservation exceeds `--max-work`, `--max-region-memory` or `--max-address-space`; the message says which ([limits](#limits)) |
| `timeout` | over `--timeout`, which is unset unless given |
| `aligner_failed` | abPOA exited with an error, or returned rows that do not match its input. A failure of the machine (abPOA's reservation refused, abPOA killed by a signal, the worker out of memory) is not recorded at all ([Memory](#memory)) |
| `no_alleles` | the package has fragments but no allele |
| `not_homologous` | an allele shares less than `--min-share` of its 15-mers with the centre; the message names the allele that shares least, its weight, and what it shares reverse-complemented ([limits](#limits)) |
| `error` | any other failure, as the exception's type and message; it means a bug |

`replace` reads only `id`, `status` and `rows`, so a program can feed it MSAs made another way
([realign_package](#realignrealign_package)). After `project`, `rows` are keyed by the recipient's
member ids, `columns` is the projected width, and every other field, `centre` included, is the
union's.

### `msas.jsonl.gz.partial`

The records `realign` has finished, in the same form, each its own gzip member, appended in the
order the regions finish. `realign` reads it back up to the first record that does not decompress or
parse, which an interruption can leave cut off, and rewrites it without that record before carrying
on. A finished run deletes it.

### replace's summary

`replace` prints its summary, stores it in its manifest, and under `run` writes it to
`replace.json`:

| key | meaning |
|---|---|
| `regions` | the number of regions replaced |
| `nodes_removed`, `nodes_added` | the number of interior nodes removed, and of induced-graph nodes added |
| `edges_removed` | the number of input edges dropped because they touched a removed node |
| `runs_rewritten` | the number of runs rewritten |
| `walks_rewritten` | the number of walks with at least one rewritten run |
| `id_mode` | `dense` or `reuse` |
| `id_range` | `dense`: the first and last node id of the output |
| `next_id` | `reuse`: the next unused id from `--id-start` |
| `ids_above_start` | `reuse`: the number of ids taken from `--id-start` on |
| `gbz` | the GBZ built, when `--gbz` is given |

A count that is zero is left out.

### verify's report

| key | meaning |
|---|---|
| `before`, `after` | `nodes`, `edges` and `length` of each graph, from `vg stats -lz` |
| `metadata_identical` | whether the path metadata matches |
| `paths`, `bases` | the number of paths of `before`, and their total length |
| `paths_missing`, `paths_extra`, `paths_differ` | the number of paths missing from `after`, extra in it, and in both with different sequences |
| `paths_missing_names`, `paths_extra_names`, `paths_differ_names` | up to `verify.SHOW` (20) names of each, the first in sorted order |
| `ok` | true when the metadata matches and all three counts are 0 |

### `qc.tsv` and `qc.tsv.json`

`qc.tsv` has a plain header line, without `#`, and one row per measured region, in package order. An
aligned base pair is a pair of bases, one from each sequence, in one column of an alignment,
matching or not. The `pair_*`, `ref_*` and `forced*` columns are summed over the region's measured
pairs, a banded pair included and a skipped pair left out. The `graph_*` columns are present only
with `--graph`.

| column | meaning |
|---|---|
| `id` | the region's id |
| `alleles` | the number of alleles |
| `longest` | the longest allele's length |
| `sequence` | the sum over columns of the distinct bases in each, the bases of the region's induced graph; `sequence / longest` is the alignment's growth |
| `pair_tp`, `pair_opt`, `pair_msa` | over the random pairs: the aligned base pairs that both alignments have, that the optimal one has, and that the MSA has |
| `ref_tp`, `ref_opt`, `ref_msa` | the same, for the reference allele against the heaviest other alleles |
| `ref_sv`, `ref_sv_opt` | runs of 50 or more gap columns opposite bases, in the implied and optimal reference alignments |
| `forced`, `forced_opt` | columns of the implied and optimal alignments (random and reference pairs) in a forced-homology window |
| `banded` | the number of random and reference pairs compared with a banded alignment because their optimal alignment would exceed `--max-cells` |
| `skipped` | the number of random and reference pairs not compared, because the memory estimate of their alignment exceeds `--max-pair-memory` ([qc](#qc)) |
| `graph_pair_tp`, `graph_pair_msa` | over the random pairs: the aligned base pairs that both the graph and the optimal alignment have, and that the graph has |
| `graph_ref_tp`, `graph_ref_msa` | the same, for the reference pairs |
| `graph_sequence` | the bases of the nodes the region's runs visit in the graph |

A forced-homology window is 40 consecutive columns, all-gap columns aside, of which at least 20 have
a base in both rows, and in which more than 30% of those columns are mismatches.

`qc.tsv.json` holds `regions`, the number of rows, the sum of each numeric column, and these ratios,
rounded to four places and null when their denominator is 0:

| key | value |
|---|---|
| `pair_precision`, `ref_precision` | tp / msa |
| `pair_recall`, `ref_recall` | tp / opt |
| `pair_f1`, `ref_f1` | 2 tp / (opt + msa) |
| `growth` | sum of `sequence` / sum of `longest` |
| `graph_pair_precision`, ..., `graph_ref_f1`, `graph_growth` | with `--graph`, the same for the graph: its tp and msa columns over the shared `pair_opt` and `ref_opt`, and `graph_sequence` over `longest` |

It also holds `peak_estimate_gb`, the largest sum of the memory estimates of the regions measured at
once ([qc](#qc)).

#### Reading the results

- An F1 of 1 means the MSA aligns every measured pair as the pair's optimal alignment does. A pair
  that includes the centre is aligned in the MSA by its own unbanded alignment unless `realign`
  banded it, so under `realign`'s scores `ref_*` says little about a region whose centre is the
  reference allele.
- `forced` above `forced_opt` means the MSA aligns unrelated sequence base against base more than
  optimal alignments do; `ref_sv` above `ref_sv_opt` means indels the reference comparison does not
  need. `growth` is at least 1, and rises with the sequence the graphs hold beside the longest
  alleles.
- The `graph_` ratios are the region's alignment before `replace`, and the others after it, both
  against one set of optimal alignments. A graph that folds a repeat's copies into a cycle aligns
  more pairs than an alignment can, so its precision falls; its `graph_growth` can be below 1. Read
  the per-region rows for the regions where the graph does better: the totals are dominated by the
  many regions where both agree.

### `<output>.manifest.json`

A stage that writes a file names its manifest `<output>.manifest.json`; `prepare` and `extract`,
which write directories, put theirs inside, and `verify` writes one only when given `-o`. Its keys,
sorted:

| key | meaning |
|---|---|
| `stage` | the stage's name |
| `pgrealign`, `python` | their versions |
| `command` | the command line, as a list |
| `started`, `finished` | local times in ISO 8601 with the UTC offset; `started` is null if unknown |
| `inputs` | `{name: {path, bytes, mtime}}`, with the absolute path; `catalog` adds `md5` |
| `params` | the stage's parameters |
| `tools` | `{name: {path, version, error}}` for each external tool the stage can run (vg or abPOA) |
| `summary` | the stage's summary, as printed (`prepare` has `outputs`, its output paths, instead; `verify` has neither) |

A tool's `path` is its absolute path, or null when it cannot be found. Its `version` is the first
line that `vg version`, or `<tool> --version` for any other tool, prints; it is null when the tool
cannot be run, exits with an error or prints nothing, and `error` then names the cause. `catalog`
records TRF when `--trf` runs it, and no other tool; its summary names the builds of trfind and strf
instead.

| stage | manifest | `inputs` | `params` | `tools` |
|---|---|---|---|---|
| `catalog` | `<out>.manifest.json` | `repeatmasker`, `reference_0`, ... | `gap`, `contigs`, `min_period`, `max_period` | TRF with `--trf`, else none |
| `prepare` | `<out>/<contig>.manifest.json` | `gbz` | `contig`, `ref_sample`, `drop_samples` | vg |
| `regions` | `<out>.manifest.json` | `targets`, `snarls`, `ref` | `contig`, `pad`, `max_span` | vg |
| `extract` | `<out>/extract.manifest.json` | `graph`, `regions` | `ref_sample`, `max_internal_bp` | vg |
| `realign` | `<out>.manifest.json` | `packages` | `method`, `select`, `shard`, `timeout`, `max_cells`, `max_memory`, `scores`, `min_share`, `centre`, `max_work`, `max_region_memory`, `max_address_space`, `max_excess` | abpoa |
| `replace` | `<gbz>.manifest.json` with `--gbz`, else `<out>.manifest.json` | `graph`, and `msas` or `msas_0`, ... | `extract`, `id_mode`, `id_start` | vg |
| `verify` | `<out>.manifest.json`, only with `-o` | `before`, `after` | none | vg |
| `qc` | `<out>.manifest.json` | `msas`, and `graph` with `--graph` | `extract`, `pairs`, `refs`, `sample`, `scores`, `max_cells`, `max_memory`, `max_pair_memory` | abpoa, and vg with `--graph` |
| `merge` | `<out>.manifest.json` | `gbz_0`, ... | none | vg |
| `union` | `<out>.manifest.json` | `base`, `add` | none | none |
| `project` | `<out>.manifest.json` | `union`, `msas` | none | none |

`max_memory` is the budget, in GB, that the stage ran within: half the machine's physical memory
unless `--max-memory` was given.

`run` writes the manifests at the paths in [its table](#run). They are these stages' manifests,
except that `prepare`'s is named `C.gbz.manifest.json` rather than `C.manifest.json`, that their
`command` is the `run` command line, and that they have fewer `params`: `regions` records only
`contig`, `realign` leaves out `select` and `shard`, and `qc` leaves out `pairs`, `refs` and
`sample`.

## The Python API

The alignment method lives in `pgrealign/star.py`. Its functions take sequences and return rows;
they read no graph and no pipeline file. A graph builder such as Minigraph-Cactus can call it once
per region: collect the region's distinct sequences with their weights, as alleles, and any partial
sequences, as fragments; call `star.star_msa` with an aligner, an `Abpoa` entered with `with`
([example](#example-abpoa)) or any callable that keeps the [aligner contract](#the-aligner-contract);
and use the returned rows, or turn them into a graph with `induce.induce`.

The other functions below run whole stages or parts of them. In the signatures, an argument written
`...` defaults to the default of the matching command-line option (`pgrealign <stage> --help`). A
`vg` of None means vg on `PATH`, and an `abpoa` of `'abpoa'` abPOA on `PATH`.

### `star.star_msa`

```python
star_msa(alleles, fragments=(), align=None, centre=None, k=K, deadline=None, max_cells=None,
         min_share=None, scores=SCORES)
```

| argument | meaning |
|---|---|
| `alleles` | a sequence of `(id, sequence, weight)`; at least one, with distinct ids |
| `fragments` | a sequence of `(id, sequence, side)`, side being `left`, `right` or `internal` |
| `align` | the aligner (required) |
| `centre` | an allele id to use as the centre instead of the medoid |
| `k` | the k-mer size for choosing the centre; `star.K` is 15 |
| `deadline` | a `time.time()` value; choosing the centre past it raises `TimedOut` |
| `max_cells` | an alignment to the centre whose matrix, the product of the two sequences' lengths, would exceed this is made by `star.long_pair`, with a band and seeded, and every other one unbanded; None, the default, aligns every one unbanded |
| `min_share` | once the centre is chosen, raise `NotHomologous` if an allele of positive weight and at least `star.HOMOLOGY_BP` bases shares less than this fraction of its k-mers with the centre (`star.share`), the centre having that many bases too; None, the default, or 0 checks nothing |

It returns a [`StarMSA`](#starmsa). Its rows spell the inputs upper-cased. The aligner sees only A,
C, G, T and N: `star_msa` gives it any other character as N and restores the character in the rows.

`star_msa` builds the star that [method.md](method.md#realign-the-star) describes. The
centre's distances, one minus the Jaccard index of two alleles' k-mer multisets, are exact unless
comparing every pair would take more than `star.MEDOID_WORK` (10^9) steps; they are then estimated
from `star.sketch`, each multiset's `star.SKETCH` (1,000) smallest hashes. `star_msa` drops the
columns that are gaps in every row, and checks that every row spells its input and that all rows
have the same length.

### `StarMSA`

| attribute | meaning |
|---|---|
| `rows` | `{id: row}` for every allele and fragment, all the same length, `-` for a gap |
| `centre` | the centre allele's id |
| `slot_alignments` | the number of insertion slots that needed a call to the aligner |
| `fragments_exact` | the number of fragments placed in an allele's row |
| `fragments_aligned` | the number of unplaced fragments, aligned to the centre |
| `banded_pairs` | the number of alignments to the centre made by `star.long_pair` because of `max_cells` |
| `seeded_pairs` | the number of those whose seeded alignment scored higher and was kept |
| `refined_windows` | the number of windows in which refinement improved a kept seeded alignment |

### The aligner contract

`star_msa` calls its aligner as `align(seqs)`, `align(seqs, banded=...)` or `align(seqs, banded=True,
seeded=True)`, and so does `qc.region_qc`. The aligner returns a global alignment of `seqs`.

- `seqs` is a list of two or more distinct, non-empty, upper-case strings. From `star_msa` they hold
  only A, C, G, T and N; `qc.region_qc` passes alleles as they are, upper-cased.
- From `star_msa`, a call is either an alignment to the centre, of the centre and another allele or
  an unplaced fragment, the centre first, or the insertions of one slot, two or more, longest
  first. `star_msa` aligns an allele equal to the centre, or an empty allele or centre, without
  calling the aligner.
- `banded=True` marks a slot's insertions. An alignment to the centre whose matrix would exceed
  `max_cells` is asked for twice (`star.long_pair`): with `banded` a number, the band's margin
  either side (`star.band_width`), and with `banded=True, seeded=True`; where the seeded one scores
  higher under `scores`, the windows that refine it are asked for unbanded. Every other alignment to
  the centre is unbanded.
- From `qc.region_qc`, a call is the two alleles of a measured pair, asked for in the same ways when
  their lengths multiply to more than `max_cells`, and unbanded otherwise: the reference allele first
  when the pair holds it, and otherwise the allele that comes first in the package.
- The aligner returns a list of strings, one row per sequence in order, all the same length, with
  `-` for a gap; each row without its gaps is its sequence. Columns that are gaps in every row are
  allowed.

`banded` and `seeded` are hints that the aligner may ignore; `Abpoa` uses them to choose abPOA's
band and seeding. The aligner owns its scoring, since `star_msa` passes it none, but `star_msa`
scores the two alignments of a long pair, and a refined window, with its `scores`, which should be
the aligner's; an aligner that ignores the hints returns the same alignment twice, and the band's,
the first, is kept. `deadline` bounds only the choice of the centre, so an aligner that must stop in
time checks the clock itself and raises `TimedOut`, as `Abpoa` does.

A pairwise alignment whose rows do not spell their inputs makes `star_msa` raise `AlignerFailed`;
one that is not two rows raises `ValueError`, and two rows of unequal length `IndexError` or
`AlignerFailed`. A slot alignment is checked only by the final check, which raises `ValueError`,
and one with too few rows raises `KeyError`. `realign` records every exception other than a
`RealignError` as the status `error`.

### `Abpoa`

```python
Abpoa(exe='abpoa', scores=SCORES, deadline=None, workdir=None)
```

`Abpoa` is the aligner `realign` and `qc` use: it runs the abPOA program once per call. Use it as a
context manager, which owns a temporary directory, created in `workdir` or else in the system's
temporary directory. Each call runs `abpoa -m 0 -r 1` (global alignment, MSA output) followed by
`flags(seqs, banded, seeded)`, and returns abPOA's rows upper-cased.

- `flags(seqs, banded, seeded=False)` returns every flag after `-m 0 -r 1`: `scores.flags()`, plus
  `-b -1` (no band) unless `banded` is true. abPOA's adaptive band then applies: its own, 10 cells
  plus 1% of the second sequence's length either side, if `banded` is True, and `banded` cells
  either side, `-b <banded> -f 0`, if it is a number. With `seeded` and `banded` true, the flag is
  `-S` instead, abPOA's own band within each window. A subclass can override it to choose flags per
  call; its flags replace these, so it must include the score flags and `-b -1` where it wants
  them.
- Past `deadline`, a `time.time()` value, a call raises `TimedOut`, before starting abPOA or by
  stopping it. When abPOA exits with an error, or returns no row for some input, a call raises
  `AlignerFailed`.
- Each call writes `in.fa` and `out.fa` in the object's directory, so one `Abpoa` must not be called
  from two threads at once. Give each thread or process its own.

### `Scores` and `SCORES`

```python
Scores(match=2, mismatch=4, gap_open=(4, 24), gap_ext=(2, 1))
```

`Scores` holds abPOA's scores: a match score, a mismatch penalty, and two-piece gap penalties, open
(o1, o2) and extension (e1, e2). A gap of length n costs min(o1 + n e1, o2 + n e2). `Scores()` holds
abPOA's own defaults. `flags()` returns the scores as abPOA's flags `-M`, `-X`, `-O` and `-E`, and
`Scores.parse(text)` reads such flags, as `--scores` takes them; a flag left out keeps abPOA's
default, and an unknown flag or a missing value raises `ValueError`. Two `Scores` are equal when
their flags are.

`star.SCORES`, the scores `realign` and `qc` use unless `--scores` is given, is
`Scores(match=8, mismatch=16, gap_open=(16, 96), gap_ext=(8, 1))`: abPOA's defaults times four,
except the long-gap extension, which stays 1 ([method.md](method.md#scores) gives the reason).

### Exceptions

`star.RealignError` is the base class of the exceptions for a region that cannot be aligned. Its
`status` attribute is the status that `realign` writes to `msas.jsonl.gz` for such a region.

| class | `status` | raised when |
|---|---|---|
| `RealignError` | `error` | `star_msa` is given no allele |
| `TooBig` | `too_big` | `realign.realign_package`'s size checks fail ([limits](#limits)); `star_msa` does not raise it |
| `TimedOut` | `timeout` | choosing the medoid runs past `deadline`, or `Abpoa` runs past its own |
| `AlignerFailed` | `aligner_failed` | abPOA fails, or a pairwise alignment does not spell its inputs |
| `ReservationRefused` | `reservation_refused` | abPOA fails because the system refused its reservation (`posix_memalign fail`); a subclass of `AlignerFailed`, which `realign.run` does not record but stops on ([Memory](#memory)) |
| `AlignerKilled` | `aligner_killed` | abPOA is killed by a signal, as by the OOM killer; a subclass of `AlignerFailed`, which `realign.run` does not record but stops on |
| `NotHomologous` | `not_homologous` | `star_msa` is given `min_share`, and an allele shares less than that with the centre (`star.check_homology`) |

`realign` never passes `star_msa` a package without alleles: it records `no_alleles` for one itself.
`star_msa` raises `ValueError` for the other arguments it cannot use (no aligner, repeated ids, an
unknown side, a centre that is not an allele) and when its final check fails.

### Example: a custom aligner

This example is runnable with pgrealign alone. Its aligner aligns two sequences with Python's
`difflib`, and pads three or more on the right, which is a valid if crude global alignment. It
ignores the `banded` hint, as the contract allows, so a slot with exactly two insertions goes through
`difflib` too.

```python
import difflib

from pgrealign.induce import induce
from pgrealign.star import star_msa


def align(seqs, banded=False):
    """A global alignment of `seqs`: one row per sequence, in order, all the same length."""
    if len(seqs) > 2:
        # the distinct insertions that share one slot, longest first: pad them on the right
        return [s.ljust(len(seqs[0]), '-') for s in seqs]
    a, b = seqs
    rows_a, rows_b = [], []
    for _, i1, i2, j1, j2 in difflib.SequenceMatcher(None, a, b, autojunk=False).get_opcodes():
        width = max(i2 - i1, j2 - j1)
        rows_a.append(a[i1:i2].ljust(width, '-'))
        rows_b.append(b[j1:j2].ljust(width, '-'))
    return [''.join(rows_a), ''.join(rows_b)]


unit = 'ACGTTGCAAGT'
alleles = [('a1', 'GATTACA' + unit * 4 + 'CCGG', 5),     # (id, sequence, weight)
           ('a2', 'GATTACA' + unit * 6 + 'CCGG', 1),
           ('a3', 'GATTACA' + unit * 2 + 'CCGG', 1),
           ('a4', '', 1)]                                 # deletes the whole region
fragments = [('f1', 'GATTACA' + unit * 5, 'left'),        # (id, sequence, side)
             ('f2', unit + 'CCGG', 'right')]

msa = star_msa(alleles, fragments, align=align)
print('centre', msa.centre)
for member, row in msa.rows.items():
    print('%-3s %s' % (member, row))

graph = induce(msa.rows)
print('%d nodes; a3 walks %s' % (len(graph.seq), graph.walks['a3']))
```

It prints:

```
centre a1
a1  GATTACA----------------------ACGTTGCAAGTACGTTGCAAGTACGTTGCAAGTACGTTGCAAGTCCGG
a2  GATTACAACGTTGCAAGTACGTTGCAAGTACGTTGCAAGTACGTTGCAAGTACGTTGCAAGTACGTTGCAAGTCCGG
a3  GATTACA--------------------------------------------ACGTTGCAAGTACGTTGCAAGTCCGG
a4  -----------------------------------------------------------------------------
f1  GATTACAACGTTGCAAGTACGTTGCAAGTACGTTGCAAGTACGTTGCAAGTACGTTGCAAGT---------------
f2  --------------------------------------------------------------ACGTTGCAAGTCCGG
5 nodes; a3 walks [1, 4, 5]
```

Left normalisation ([method.md](method.md#left-normalisation)) put each allele's extra or missing
copies of the repeat unit at the start of the repeat. `f1` was found in `a2` and `f2` in `a1`, so
neither fragment needed aligning.

### Example: abPOA

This example needs `abpoa` on `PATH`. It gives the region a time limit, asks for any alignment to
the centre larger than 10,000,000 cells to be banded, and treats a `RealignError` as `realign`
does: the region is left as it is. The deadline goes to both `Abpoa`, which applies it to each
abPOA run, and `star_msa`, which applies it to the choice of the centre.

```python
import time

from pgrealign.star import SCORES, Abpoa, RealignError, star_msa

unit = 'ACGTTGCAAGT'
alleles = [('a1', 'GATTACA' + unit * 4 + 'CCGG', 5),
           ('a2', 'GATTACA' + unit * 6 + 'CCGG', 1),
           ('a3', 'GATTACA' + unit * 3 + 'GCGG', 2)]
fragments = [('f1', 'GATTACA' + unit * 5, 'left')]

deadline = time.time() + 60                     # give up on this region after a minute
try:
    with Abpoa('abpoa', scores=SCORES, deadline=deadline) as align:
        msa = star_msa(alleles, fragments, align=align, deadline=deadline, max_cells=10_000_000)
except RealignError as e:                       # TimedOut or AlignerFailed
    print('left as it is:', e.status, e)
else:
    print('centre %s, %d slot alignments, %d banded' % (msa.centre, msa.slot_alignments,
                                                         msa.banded_pairs))
    for member, row in msa.rows.items():
        print('%-3s %s' % (member, row))
```

### `realign.realign_package`

```python
realign_package(pkg, abpoa='abpoa', timeout=None, max_cells=..., workdir=None, scores=SCORES,
                min_share=..., centre='ref', max_work=..., max_region_memory=..., max_address_space=...)
excess(pkg, aid)
centre_of(pkg, centre='ref', max_excess=...)
work(pkg, centre='ref', max_cells=...)
peak_bytes(pkg, max_cells=..., page=ESTIMATE_PAGE, centre='ref')
reservation(pkg, centre='ref')
refusal(pkg, max_cells=..., centre='ref', max_work=..., max_region_memory=..., max_address_space=...)
plan(packages_path, select='any', shard=(0, 1), max_cells=..., centre='ref', max_work=...,
     max_region_memory=..., max_address_space=...)
run(packages_path, out_path, select='any', shard=(0, 1), jobs=1, max_memory=None, **kw)
```

`realign_package` aligns one package with `star_msa` and `Abpoa`, and returns its record in the
[`msas.jsonl.gz` format](#msasjsonlgz). The package is a dict in the [package
format](#packagesjsonlgz), of which it reads `id`, `ref`, `alleles` and `fragments`. It applies
`realign`'s [limits](#limits), passing `max_cells` and `min_share` to `star_msa` and the centre that
`centre_of` names, and records a failure to align in the record's status rather than raising it:
`no_alleles` for a package without alleles, `too_big` when `refusal` gives a reason, the status of
any `RealignError` that `star_msa` or `Abpoa` raises, `reservation_refused` among them, and `error`
for any other exception. `timeout` None sets no time limit, and `max_work`, `max_region_memory` and
`max_address_space` 0 or None set no such limit. `realign.write_msas(path, records)` writes records
to a `msas.jsonl.gz`.

`excess` is the share of the haplotypes' bases that lie beyond allele `aid`'s length, each allele
weighted by its haplotypes. `centre_of` is the id of the allele the star is built around, the reference
allele when `centre` is 'ref', the package has one of positive weight and its excess is no more than
`max_excess` (None for no limit), or None for the medoid; the estimates, `refusal`, `plan`,
`realign_package` and `run` take `max_excess` too. `work`, `peak_bytes` and
`reservation` are a package's [work](#limits), [memory estimate](#memory) in bytes, `page` being the
page size, and abPOA's largest reservation in bytes, and `refusal` is the message of the first
limit the package exceeds, or None. `plan` is `realign --plan`. `run` is the `realign` stage: it
passes `kw` (`abpoa`, `timeout`, `max_cells`, `workdir`, `scores`, `min_share`, `centre` and the
three limits) to `realign_package`, records the packages `refusal` gives a reason for as `too_big`
itself, runs the others within `max_memory` GB, None meaning half the machine's physical memory,
raises `RuntimeError` before aligning anything when one of them alone needs more, and returns the
summary. It also raises `RuntimeError` before aligning anything when a region it admits would need a
larger reservation than `memory.granted_bytes()`; a region whose status is one of
`realign.MACHINE_FAILURES` is left out of the partial file, and `run` raises `RuntimeError` once the
others are done. `refusal` and the estimates count only the fragments that occur in no allele, since
`star_msa` places the others without aligning them.

`realign_package` always uses abPOA. To put another aligner into the pipeline, call `star_msa` with
it for each `ok` package and write records with `realign.write_msas`, each with `id`, `status` set
to `ok`, and `rows`: one row per member, keyed by member id, all the same length, each spelling its
member's sequence upper-cased.

### `induce.induce` and `Graph`

```python
induce(rows, max_node=MAX_NODE)
```

`induce` builds the induced graph of an MSA, `{member: row}`
([method.md](method.md#the-graph-of-an-alignment)), with nodes cut into pieces of at most
`max_node` bases. `induce.MAX_NODE` is 1,024, vg's default, so the graph needs no id translation
when it becomes a GBZ. It returns a `Graph`, a forward-only acyclic graph with ids 1 to n in
topological order:

| attribute | meaning |
|---|---|
| `seq` | `{node id: sequence}` |
| `edges` | a set of `(from, to)` pairs, each from the end of one node to the start of the next |
| `walks` | `{member: [node id, ...]}`; an all-gap row has an empty walk |

`induce` checks that every walk spells its row's characters, and raises `ValueError` if one does
not. It treats characters as they are, so an upper-case and a lower-case base in one column become
two nodes.

### `qc.region_qc`

```python
region_qc(pkg, rows, align, pairs=..., refs=..., max_cells=None, graph=None, max_bytes=None)
```

`region_qc` measures one region's MSA, `rows`, as [`qc`](#qc) does, and returns a dict keyed by
`qc.FIELDS`, the columns of [`qc.tsv`](#qctsv-and-qctsvjson). It reads the package's `id`, `alleles`
and `ref`, if present, and needs at least one allele. `align` computes the optimal pairwise
alignments under the [aligner contract](#the-aligner-contract); pass an `Abpoa` with the scores the
MSA was built with. A pair whose two lengths multiply to more than `max_cells` is aligned with a
band and counted in `banded`; None, the default, bands none. `graph`, if given, is the graph's own
alignment of the region, as `(walks, lengths, interior_bp)`: `walks` maps each allele id to its walk
as a list of `('>' or '<', node id)` in the region's forward orientation, `lengths` maps each node
id to its length, and `interior_bp` is `graph_sequence`; the dict then has `qc.GRAPH_FIELDS` too. A
pair whose estimate, `qc.pair_bytes(n, m, max_cells)` for alleles of n and m bases plus
`qc.GRAPH_BYTES` per base of the alleles when `graph` is given, exceeds `max_bytes` is not compared
and counts in `skipped`; None, the default, skips none. `qc.run` takes it as `max_pair_memory`, in
GB. `qc.compared_pairs(pkg, rows, pairs, refs)`
lists the pairs `region_qc` compares, in order, as `(kind, a, b)`, and
`qc.peak_bytes(pkg, rows, pairs, refs, max_cells, max_bytes, graph)` is the largest estimate among
those it does not skip ([qc](#qc)), `graph` being true when `region_qc` is given one.
`qc.graph_walks(extract_dir, graph, region_ids, vg=None, threads=4)` reads the graph's walks for a
set of regions from `extract`'s output and its graph. `qc.summarise(list_of_dicts)` returns the sums and ratios of `qc.tsv.json`, and
`qc.run(extract_dir, msas_path, out_path, ...)` is the `qc` stage, with keyword arguments named
after its options, `max_cells`, `max_memory` (in GB, None meaning half the machine's physical
memory), `graph`, `vg` and `threads` included, and `workdir` for its temporary files.

### `pipeline.run_contig` and `merge.merge`

```python
run_contig(gbz, contig, targets, out_dir, jobs=..., threads=..., vg=None, abpoa='abpoa',
           ref_sample='CHM13', drop_samples=None, scores=SCORES, timeout=None, max_cells=...,
           max_internal_bp=..., id_start=None, do_qc=True, max_memory=None, log=print, min_share=...,
           centre='ref', max_work=..., max_region_memory=..., max_address_space=...)
merge(gbzs, out_gbz, vg=None, workdir=None)
```

`run_contig` is [`run`](#run): it runs every stage for `contig` of the whole-genome `gbz` under
`out_dir/<contig>/`, skipping stages whose outputs exist. `drop_samples` None removes `prepare`'s
default samples, and an empty sequence removes none. `max_memory` None means half the machine's
physical memory, and a value not more than 0 raises `ValueError` before any stage runs. `log`
receives one line of progress for each stage it runs. It returns the dict that `run` prints, and
raises `RuntimeError` when `verify` fails.

`merge` is [`merge`](#merge): it joins the GBZs in `gbzs` into `out_gbz` and returns the summary.
`workdir` holds the temporary GFA, in the system's temporary directory when None. It raises
`ValueError` when two inputs' node ids overlap, before building anything, and when the output fails
its check, after renaming it to `<out_gbz>.failed`.

## Known gaps

- `replace` rewrites walks (W lines) only, because `extract` records runs of walks only. A P line
  that visits a removed node stops `replace`. A reference path given as a P line does so in every
  replaced region whose interior nodes it shares with a haplotype's walk.
- A walk that visits no anchor and spells more than `--max-internal-bp` bases is never recorded as a
  run, and `replace` cannot rewrite a walk through an interior outside a recorded run. So `extract`
  makes a region whose interior such a walk uses `complex`, and the region is left as it is.
- `star_msa` upper-cases every row, so lower-case (soft-masked) bases in a replaced region come out
  upper-case. `verify` upper-cases both graphs' paths before comparing them, so it does not report
  this.
- With `prepare`'s default, the realigned graph has no reference path cover. To add one, compute
  it on a mutable copy of each realigned contig, and of each contig cut with `prepare`, before
  `merge`, and rebuild the GBZ; pgrealign does not run these commands. `vg paths -u` writes the
  cover as the sample `gref_<ref>`, which the rebuilt GBZ lists as a reference sample:

  ```sh
  vg convert -p chr1.realigned.gbz > g.pg
  vg paths -x g.pg -u -Q <ref> > g.gref.pg
  vg convert -f g.gref.pg > g.gref.gfa
  vg gbwt -G g.gref.gfa --gbz-format -g chr1.realigned.gref.gbz
  ```

- abPOA reserves address space for the whole matrix of each alignment, banded or not, rounded up
  to a power of two: 2^37 bytes, 128 GB, for two sequences of 62,146 bp. macOS grants it and only
  the pages touched take memory; Linux refuses a reservation larger than its memory and swap.
  `--max-address-space` keeps the outcome the same on every machine by refusing, on all of them, the
  regions whose reservation would exceed it, but a machine with less than that limit in memory and
  swap stops `realign` at the first such region it admits ([Memory](#memory)).
- The work, memory and address-space estimates take every member's excess over the reference to
  fall in one slot and be aligned to the longest excess. A member shorter than the reference can
  insert bases too, and abPOA aligns each of a slot's insertions to the graph it has grown from the
  ones before, which can be longer than the longest alone, so a slot can take more than the
  estimates say: they are judgements, not guarantees.
- `realign`'s memory estimate around a medoid takes the worst centre, since the medoid is chosen
  only once the region runs. For a region whose alleles differ in length by tens of kb it can be
  twice what the region's own centre needs, as for the 58.7 kb centre and 103.2 kb allele of
  [Memory](#memory), so such a region can be left `too_big` though it would fit.
- Each `realign` or `qc` process keeps to its own `--max-memory`. Several on one machine at once,
  such as the contigs or shards of a genome, each need `--max-memory` set to their share.
- On the command line, `realign` and `qc` align with abPOA only. Another aligner needs the Python
  API and an `msas.jsonl.gz` written with `realign.write_msas`.
- `qc`'s pairwise alignments have no time limit: `realign`'s `--timeout` does not apply to them.
  `qc` measures a pair over `--max-cells` against a banded alignment, which may not be optimal, so
  in a region of long alleles its yardstick is weaker; the `banded` column shows how many pairs that
  applies to.
- `union` keeps the panel's status, so a region invariant in the panel is not realigned even when
  the recipient has alleles the panel lacks. The recipient's `replace` changes only regions that are
  `ok` in the recipient's own `extract`.
- `replace`'s `dense` numbering starts at the contig's smallest input id, so a contig that grows can
  run into the next contig's ids, and how many ids a contig needs is known only once `replace` has
  run. `merge` stops when this happens. Run `replace` on the contig `merge` names first again, with
  `--id-start` above every other input's largest id (`vg stats -r` prints a GBZ's range), and then
  `merge` again; under `run`, delete `C.realigned.gbz` and `verify.json` first.
  [whole-genome.md](whole-genome.md) gives the commands.
- `merge` compares node-id ranges, not ids. A contig replaced under `reuse` that took ids from an
  `--id-start` above other contigs' ids has a range that spans theirs, so `merge` refuses it.
- `run` decides what to skip only by whether files exist, so it does not notice changed inputs or
  options; delete the outputs of the stage that changed and of every stage after it.
- `realign` reuses the records in `<out>.partial` without checking the options that made them, except
  the limits on a region's size, which it applies afresh, so a run with changed `--scores`,
  `--centre`, `--max-cells`, `--min-share` or `--timeout` mixes old and new results. Delete
  `<out>.partial` after changing those.
