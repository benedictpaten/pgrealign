# Building a whole-genome realigned graph

pgrealign realigns chosen intervals of a pangenome graph one contig at a time. This page describes
how to build a realigned graph of a whole genome: realign each contig that has intervals to realign,
cut out the parts of the graph that are not realigned, merge them all into one graph and check the
result. pgrealign does this work in stages, each a subcommand of the program. [The
reference](reference.md) describes every stage, option and file, and [the method](method.md) how
each stage works; `pgrealign <stage> --help` lists every option, with the default of most.

## Overview

### Paths

A GBZ is vg's compressed file format for a pangenome graph together with its paths, which it holds
in a GBWT, vg's compressed index of paths. A node of the graph holds a sequence and is named by an
integer, its node id. A path is a named route through the nodes: a reference path, of a sample that
the GBZ marks as a reference sample, such as CHM13; a haplotype path, which covers a stretch of one
haplotype; or a generic path, which has only a name. `vg paths -M -x <gbz>` prints a table of every
path, whose first column is its name and fifth its contig.

### Contigs

pgrealign's reference sample (`--ref-sample`) is the reference sample whose paths give pgrealign its
coordinates. A contig is one of that sample's sequences, such as chr1, with the path
`<sample>#0#<contig>`, and its component is the connected component of the graph that holds it. vg
divides a component into snarls: a snarl is a part of the graph cut off from the rest by two
boundary nodes, so that a path can enter or leave it only through them
([snarls and chains](method.md#snarls-and-chains)).

### Regions

A target is an interval of a contig to realign. pgrealign places each target in a region, the part
of the graph between two anchors, snarl boundary nodes that every path crossing the target passes
through. A region's alleles are the distinct sequences that the paths crossing it spell, and its
fragments the partial ones spelled by paths that start or end inside it. pgrealign aligns the other
alleles, and the fragments found in no allele, to one allele, the centre, with abPOA, a sequence
aligner, merges those alignments into one, and replaces the region's nodes with a graph built from
it. Every path spells the same sequence afterwards as before, apart from letter case. [The
reference's overview](reference.md#overview) defines these terms at more length.

### Stages

For each contig, prepare cuts the contig's component out of the whole-genome GBZ; regions chooses
the regions; extract collects each region's alleles and fragments; realign aligns them; replace
writes the realigned graph; verify checks it against the prepared one; and qc measures the
alignments. `pgrealign run` runs these seven stages for one contig, and merge joins the contigs into
one graph. regions, extract and realign each give a region a status, and only a region that is `ok`
in all three is replaced.

## Inputs

`pgrealign run` takes the whole-genome GBZ (`--gbz`), such as the HPRC v2.1 Minigraph-Cactus GBZ,
the targets (`--targets`), the reference sample and the samples to drop. The commands on this page
use CHM13 as the reference sample; with another, change it in each.

### The reference sample's paths

The reference sample must have exactly one path on each contig, not one for each piece of a contig
assembled in pieces: prepare stops when it finds none or several. It cuts a contig's component with
`vg chunk --contig <contig>`, which selects every component that holds a path whose contig is that
name, of any sample, and stops unless that gives exactly one. It does not check that the component
holds the path of only one of the reference sample's contigs; see
[why two ranges overlap](#why-two-ranges-overlap).

### Targets

The targets are a BED file of intervals, such as the catalogue `catalogs/chm13v2.0_tr.bed.gz` of the
tandem repeats of CHM13 v2.0 on chr1-22, chrX and chrY; `pgrealign catalog` builds one for another
reference ([catalog](reference.md#catalog)). The BED's first column holds the contig field of the
reference sample's path names, such as `chr1`, not the whole name, `CHM13#0#chr1`.

### Dropped samples

A reference sample's cover, which `vg paths -u` computes from its paths, is a sample named
`gref_<sample>`: a copy of each reference path, such as `gref_CHM13#0#chr1`, plus short paths that
branch off it, such as `gref_CHM13#0#chr1_1_alt`. The short paths start and end inside regions and
would add fragments to align, so by default prepare removes CHM13's cover, gref_CHM13, from each
contig's graph right after cutting it out ([dropped samples](reference.md#dropped-samples)).
`--drop-sample` replaces the default rather than adding to it, so to drop another sample as well as
the cover, name both; with another reference sample, name its cover. `--keep-all-samples` removes
none.

## Running one contig

```
pgrealign run --gbz hprc-v2.1.gbz --contig chr1 --targets catalogs/chm13v2.0_tr.bed.gz \
    -o out -j 16 -t 8
```

`run` runs every stage for the contig in order, and writes under `out/chr1/`:
- prepare's `chr1.gbz`, the contig's component; `chr1.snarls.pb`, its snarls; and `chr1.ref.tsv.gz`,
  its reference index, which gives the position on the contig of each node of the reference path;
- regions' `regions.tsv.gz`;
- extract's directory `extract/`, whose `packages.jsonl.gz` holds a package for each region that
  regions made `ok`, with its anchors, alleles, fragments and status;
- realign's alignments, `msas.jsonl.gz`;
- replace's `chr1.realigned.gbz` and `replace.json`;
- verify's report, `verify.json`;
- qc's `qc.tsv` and `qc.tsv.json`.

Each stage also writes a manifest, a JSON record of its inputs, options and tools;
[run](reference.md#run) lists every file. At the end `run` prints a JSON summary to stdout, with
`id_range`, the smallest and largest node id of `chr1.realigned.gbz`.

### Options that run does not pass

`run` lacks some stage options, such as regions' `--max-span` and realign's `--select`, and those
stages then keep their defaults; [run](reference.md#run) lists them. To change one, run that stage
and those before it with the subcommands, writing the files that `run` writes and giving each the
values `run` would, and then let `pgrealign run` do the rest. prepare writes into the directory that
its `-o` names, where `run` adds the contig's name. The prepare subcommand names its manifest
`<contig>.manifest.json`; `run` looks for `<contig>.gbz.manifest.json` instead, and adds it beside
the other, recording its own command line. To realign only the regions whose alleles differ in
length, for example:

```
d=out/chr1
pgrealign prepare --gbz hprc-v2.1.gbz --contig chr1 -o $d -t 8
pgrealign regions --targets catalogs/chm13v2.0_tr.bed.gz --contig chr1 \
    --snarls $d/chr1.snarls.pb --ref $d/chr1.ref.tsv.gz -o $d/regions.tsv.gz
pgrealign extract --graph $d/chr1.gbz --regions $d/regions.tsv.gz -o $d/extract -t 8
pgrealign realign --packages $d/extract/packages.jsonl.gz -o $d/msas.jsonl.gz -j 16 --select length
pgrealign run --gbz hprc-v2.1.gbz --contig chr1 --targets catalogs/chm13v2.0_tr.bed.gz \
    -o out -j 16 -t 8
```

qc's options are simpler to change: run `pgrealign run --no-qc`, then
`pgrealign qc --extract $d/extract --msas $d/msas.jsonl.gz -o $d/qc.tsv` with the options wanted.

### Resuming

`run` skips a stage whose output exists, judging only by existence, not by the options or inputs
that made it. So running the same command again after an interruption continues from the first
unfinished stage, and an interrupted realign resumes from `msas.jsonl.gz.partial`, which holds each
region's result as it finishes. The file by which `run` judges each stage is renamed into place once
complete, after the stage's other files, so a run killed partway through a stage redoes it. Each
stage writes its manifest last, so a run killed in between leaves the stage done without one.

To redo a stage, for example with other options, delete in `out/<contig>/` its outputs and those of
every stage that `run` runs after it, with their manifests and realign's partial file, and run
again. qc reads neither replace's output nor verify's, so redoing replace or
verify keeps `qc.tsv`. To redo prepare, delete the whole directory.

## Running a genome in parallel

These commands list the reference sample's contigs, those that have targets in `contigs.txt` and
those that have none in `notargets.txt`. The last prints the contigs that have targets but are not
in the graph, such as chrY in a graph without it; their targets cannot be realigned:

```
vg paths -M -S CHM13 -x hprc-v2.1.gbz | grep -v '^#' | cut -f5 | sort -u > ref-contigs.txt
gzip -dc catalogs/chm13v2.0_tr.bed.gz | grep -v '^#' | cut -f1 | sort -u > target-contigs.txt
comm -12 ref-contigs.txt target-contigs.txt > contigs.txt
comm -23 ref-contigs.txt target-contigs.txt > notargets.txt
comm -13 ref-contigs.txt target-contigs.txt
```

Contigs are independent until merge, so each contig's `run` can be a process of its own, on one
machine or on several. Each writes only its own directory and, for abPOA's files, the system
temporary directory; two runs of the same contig must not run at once. Within a contig, `-j` sets
how many regions realign and qc work on at once, each in its own worker process, which runs abPOA
once per alignment, and `-t` sets the threads given to `vg snarls` and to `vg convert`, which
streams the graph in prepare, extract and replace. `vg chunk`, `vg gbwt`, `vg paths`, `vg stats` and
`vg view` use vg's own default thread counts. Each contig's prepare loads the whole-genome GBZ, so
memory as well as cores limits how many contigs one machine runs at once. This runs every contig in
`contigs.txt`, four at a time:

```
mkdir -p out
xargs -P 4 -I{} sh -c 'pgrealign run --gbz hprc-v2.1.gbz --contig {} \
    --targets catalogs/chm13v2.0_tr.bed.gz -o out -j 8 -t 4 >out/{}.json 2>out/{}.log' <contigs.txt
```

`out/<contig>.json` is empty unless the contig's run finished. This lists each contig whose run did
not finish or whose `verify.json` does not record a pass; `out/<contig>.log` says why, and
[the reference's known gaps](reference.md#known-gaps) list the cases, such as a generic path
through a region to be replaced, that stop replace. Run each again once the cause is fixed:

```
while read c; do [ -s "out/$c.json" ] && grep -q '"ok": true' "out/$c/verify.json" 2>/dev/null ||
    echo "$c"; done < contigs.txt
```

### Sharding a contig

To spread a slow contig over several machines, split realign into shards:
`pgrealign realign --shard i/n` aligns every n-th package, starting from the i-th, counting from 0.
The shards' outputs, joined end to end, make one valid `msas.jsonl.gz`, and `run` then skips
realign and does the rest. Every machine needs `extract/packages.jsonl.gz`, and each shard's
statuses are in its own manifest. Run prepare, regions and extract as above, and then:

```
# on each of four machines, with i = 0, 1, 2 and 3:
pgrealign realign --packages $d/extract/packages.jsonl.gz -o $d/msas.shard$i.jsonl.gz \
    --shard $i/4 -j 16

# once every shard has finished:
cat $d/msas.shard[0-3].jsonl.gz > $d/msas.jsonl.gz.tmp && mv $d/msas.jsonl.gz.tmp $d/msas.jsonl.gz
pgrealign run --gbz hprc-v2.1.gbz --contig chr1 --targets catalogs/chm13v2.0_tr.bed.gz \
    -o out -j 16 -t 8
```

## Alignment size

### Cells and bands

abPOA aligns two sequences by filling a matrix of cells, one per pair of their bases. `--max-cells`
bounds, roughly, the cells of each alignment with the centre. An alignment whose full matrix, the
product of the two sequences' lengths, would exceed that many cells is banded: abPOA's adaptive band
fills only the cells near the alignment's best path. Every other alignment with the centre is
unbanded, filling every cell. The distinct insertions that fall in one slot, between two centre
bases or at either end, are aligned to one another with a banded alignment, which `--max-cells` does
not bound.

### Regions too large to align

A region whose longest allele or fragment would need more than `--max-cells` cells even banded is
left as it is, with status `too_big`. So is a region whose alleles hold more than 40 million bases
in all, which bounds the k-mer counts that choose the centre (`MAX_TOTAL_BP` in `realign.py`; no
option changes it).

## Resources per contig

pgrealign streams each graph from `vg convert` rather than loading it whole, but several stages hold
a table per node, per region or per path, and some start a vg process that loads a graph:

| stage | holds in memory | writes to disk |
|---|---|---|
| prepare | `vg chunk`: the whole-genome GBZ; `vg snarls`: the contig's graph; then the length of every node of the contig | the contig's GBZ, snarls and reference index; a temporary directory in `out/<contig>/` for the cut component and, when samples are dropped, two GBWT files, removed at the end |
| regions | the snarls' boundary nodes and the reference index | `regions.tsv.gz` |
| extract | the sequence of every node; each region's alleles and fragments; the paths that visit no anchor, until it knows which region each lies in | `extract/`, mostly `occurrences.tsv.gz`, about one row per path per region it crosses |
| realign | per worker, abPOA's alignment with the centre (about 15 bytes per cell it fills) or the k-mer counts that choose the centre (about 80 bytes per base of the region's alleles, each allele counted once, plus about 50 per base of the longest allele when comparing every pair of alleles exactly would take too long); in the parent, the package of every region to align and every finished result | `msas.jsonl.gz.partial`, then `msas.jsonl.gz` |
| replace | every replaced region's alleles, fragments, alignment rows and new nodes; a record of every path's passage through those regions; every node id and its new id; then `vg gbwt -G`, building the realigned GBZ | `<contig>.realigned.gfa`, uncompressed, with every path as text; the GBZ built from it, after which the GFA is deleted |
| verify | `vg paths -F` on each graph; a length and an MD5 digest per path | `verify.json` |
| qc | the alignment rows of every realigned region; per worker, unbanded pairwise alignments by abPOA | `qc.tsv`, `qc.tsv.json` |

prepare's `vg chunk` needs the same memory for every contig, and replace holds the tables that grow
most with the contig and with the number of regions it replaces. realign runs up to `-j` workers,
each using the larger of its two figures in the table, plus the parent.

### qc's limits

qc skips a pair whose unbanded alignment would exceed `--max-cells` cells, which `run` passes to qc
as well as to realign, so the same limit bounds a qc worker's memory. But qc has no time limit, so a
region can take longer in qc than realign's `--timeout` allows it; `--no-qc` skips qc, and
`pgrealign qc --sample` measures a random subset of the regions.

### Disk

replace needs the most disk. Its GFA spells out every path's steps as text, so it is many times the
size of the contig's GBZ, and `vg gbwt -G` reads it while it writes the GBZ beside it. Keep every
output until the merged graph has passed its checks: a contig whose node ids must move is replaced
again, which needs its prepared `<contig>.gbz`, `extract/` and `msas.jsonl.gz`. abPOA's files, in
realign and qc, and merge's large working files go in the system temporary directory, `$TMPDIR`
when it is set, so set `TMPDIR` to a disk with room for them. prepare puts its temporary files in
the contig's directory.

## Parts that are not realigned

merge builds the genome only from the GBZs it is given, so it must also be given the contigs without
targets, listed in `notargets.txt`, such as chrM, and the components that hold no path of the
reference sample.

Cut each contig without targets with `pgrealign prepare`, given the same `--ref-sample` and
`--drop-sample` or `--keep-all-samples` options as the runs, so that it drops the same samples.
prepare writes into the directory that `-o` names, and only its `<contig>.gbz` goes into the merge:

```
while read c; do pgrealign prepare --gbz hprc-v2.1.gbz --contig "$c" -o "parts/$c" ||
    echo "prepare failed: $c" >&2; done < notargets.txt
```

`vg chunk` with `--gbz` and no `--contig` writes every component of the graph as a GBZ of its own,
`parts/part_<i>_<name>.gbz`, named after the contig of a reference or generic path in it, or
`component_<i>` when it has neither. Every component that holds a path of the reference sample is a
contig, already realigned or cut above, so this loop removes those; `vg paths -R -L` lists the
names of a GBZ's reference paths, and the loop stops if it cannot read a file:

```
mkdir -p parts
vg chunk -x hprc-v2.1.gbz --gbz -b parts/part
for f in parts/part_*.gbz; do
    names=$(vg paths -R -L -x "$f") || { echo "cannot read $f" >&2; break; }
    if printf '%s\n' "$names" | grep -q '^CHM13#'; then rm "$f"; fi
done
```

This `vg chunk` holds the whole-genome GBZ and every component cut from it in memory at once, and
the parts need about as much disk as the graph until the loop has removed the contigs.

`vg chunk` keeps every sample in every component. With the default `--drop-sample`, nothing more is
needed, because a component without a path of the reference sample holds no cover path. If other
samples are dropped, remove them from each part as prepare does, with
`pgrealign.prepare.drop_samples_from(part, samples, out, tmp, vg)`, which runs `vg gbwt -R` with
its working files in the directory `tmp`, and writes the part without them to `out`. `vg` is the vg
binary's name or path, such as `'vg'`, and the call raises an error when the part has reference
samples and dropping `samples` would leave none.

## Merging

```
pgrealign merge out/*/*.realigned.gbz parts/*/*.gbz parts/part_*.gbz -o hprc-v2.1.realigned.gbz
```

Give merge each contig's `.realigned.gbz`, not the prepared `<contig>.gbz` beside it, which holds
the same contig under its old node ids. merge renumbers no node, stops if two inputs' node-id ranges
overlap, and checks its output's sizes and path names only against its inputs
([merge](reference.md#merge)). Its working files hold every node sequence of the genome and one
GBWT per input, so, like prepare's `vg chunk`, merge works at the scale of the whole genome.

merge writes the output GBZ before its checks, and its manifest after them. A merge that fails them
leaves a GBZ that must not be used: delete it. Before merging again under the same name, delete any
manifest an earlier merge left, so that a manifest beside the output marks a merge that passed.
Indexes and annotations made for the whole-genome GBZ, such as its snarls, its distance and
minimizer indexes, or anything else keyed by node id, do not apply to the merged GBZ, because
replace renumbers every node of each realigned contig. Build its own.

### The node-id rule

merge needs each input's ids to lie in a range, from its smallest id to its largest, that overlaps
no other input's range; `vg stats -r` prints a GBZ's range. The parts cut with `vg chunk` and
prepare keep the ids of the whole-genome GBZ. replace renumbers each realigned contig densely: its
nodes take consecutive ids in their old order, starting at the contig's smallest id (or at
`--id-start`), and each region's new nodes follow the region's left anchor. Without `--id-start`,
the contig's range therefore starts where it did and ends where its node count takes it, and it
runs into the next input's range when its node count is more than the next input's smallest id
minus the contig's smallest id. merge then stops with a Python traceback whose last line is:

```
ValueError: node ids overlap: out/chr7/chr7.realigned.gbz ends at X, out/chr8/chr8.realigned.gbz starts at Y; replace one of them again with --id-start above every other input's ids
```

### Why two ranges overlap

The range of the first input named runs into the second's. Their original ranges, those of the
prepared `<contig>.gbz` of a realigned contig and of any part, tell the causes apart:
- The original ranges are disjoint. The first input is a realigned contig that grew: move it, as
  below.
- The original ranges are the same. If the two inputs hold the same contig, it was given twice, for
  example as its prepared GBZ beside its realigned one, or as a `vg chunk` part beside its realigned
  GBZ: give merge each input once. If they hold different contigs, the two contigs share a
  component. Each contig's prepared GBZ then holds the paths of both, and
  `vg paths -M -S CHM13 -x out/<contig>/<contig>.gbz | grep -v '^#' | cut -f5 | sort -u` prints
  more than one contig. Realign only one of the two, and leave the other out of the merge.
- The original ranges interleave without sharing an id; merge compares ranges, not ids. If one of
  the two inputs is a realigned contig, move it above the other; if neither is, no pgrealign
  command can move either.

### Moving a contig

To move a contig, choose an id N above the largest id of every input (`id_range` in each
`out/<contig>.json` gives a realigned contig's range, and `vg stats -r` any input's). Delete the
contig's realigned GBZ and its verify report, with their manifests, and run it again with
`--id-start`, writing its summary to `out/<contig>.json` as before:

```
rm out/chr7/chr7.realigned.gbz* out/chr7/verify.json*
pgrealign run --gbz hprc-v2.1.gbz --contig chr7 --targets catalogs/chm13v2.0_tr.bed.gz \
    -o out -j 16 -t 8 --id-start N >out/chr7.json
```

Only replace and verify do their work again, and qc's outputs stay. If the contig's cover was
rebuilt ([the gref_CHM13 cover](#the-gref_chm13-cover)), rebuild it again. Then run merge again.
When several contigs need moving, move them one at a time, each above every id of the inputs,
including the ids of those already moved.

## Checking the result

### verify

verify compares each contig's prepared GBZ with its realigned GBZ: the path metadata that
`vg paths -M` prints must be identical, and every path must spell the same sequence, ignoring case
([verify](reference.md#verify)). `run` stops with an error when a contig fails, and leaves
`verify.json`, which records the failure, so a later run stops again; the loop under
[running a genome in parallel](#running-a-genome-in-parallel) lists the contig. replace checks every path
it rewrites, so a failed verify points to a fault in pgrealign or in the vg it ran. `verify.json`
names some of the paths that are missing, extra or spell something else; keep the contig's directory
for a bug report. To finish the genome without realigning that contig, move its realigned GBZ out of
`out/<contig>/`, so that the merge command's `out/*/*.realigned.gbz` does not find it, and add its
prepared `out/<contig>/<contig>.gbz` to the merge command, which leaves the contig as it was.

### What was realigned

regions, extract and realign each count the regions of each status they give, in regions' manifest,
extract's summary and realign's manifest, or each shard's. The method defines the statuses of
[regions](method.md#regions), [extract](method.md#statuses) and
[realign](method.md#limits-and-regions-that-keep-their-old-subgraph). This totals them over every
contig:

```
python3 - <<'EOF'
import collections, glob, json
for pattern, keys in [('out/*/regions.tsv.gz.manifest.json', ['summary', 'status']),
                      ('out/*/extract/summary.json', ['status']),
                      ('out/*/msas*.jsonl.gz.manifest.json', ['summary', 'realigned'])]:
    total = collections.Counter()
    for f in glob.glob(pattern):
        d = json.load(open(f))
        for k in keys:
            d = d[k]
        total.update(d)
    print(pattern, dict(total))
EOF
```

realign's manifest also counts, under `summary.skipped`, the packages it did not align: those that
extract did not make `ok`, and those that `--select` passed over (`not_selected`). `replace.json`
gives, as `regions`, the number of regions replaced, and each `ok` record in `msas.jsonl.gz`
counts, in `banded_pairs`, its alignments with the centre that were banded.

### Statuses to act on

These statuses can be changed by an option, or point to a fault:

| status | what it means | what to do |
|---|---|---|
| `too_long` | no pair of anchors within `--max-span` encloses the target | raise regions' `--max-span`, at the cost of larger regions |
| `complex` | a path passes through the region in a way replace cannot rewrite | if it visits no anchor and spells more than `--max-internal-bp` bases, raise that |
| `too_big` | the longest allele or fragment would need more than `--max-cells` cells even banded, or the alleles hold more than 40 million bases | for the first, raise `--max-cells`, if memory allows; no option changes the second |
| `timeout` | the region took longer than `--timeout` seconds | raise `--timeout` |
| `aligner_failed`, `error` | abPOA failed, or realign hit a bug | read the region's `message` in `msas.jsonl.gz`, and report an `error` with it and the region's package |

### qc

A multiple alignment implies an alignment of every pair of its alleles, the pair's projected
alignment. qc compares a sample of projected alignments with the pairs' optimal alignments, which
abPOA computes without a band under the scores the multiple alignment was built with (`--scores`,
which `run` passes to both realign and qc). It reports precision, recall and F1 of the projected
alignments' aligned base pairs, for random pairs of alleles other than the reference allele
(`pair_*`) and for the reference allele, the one the reference path spells, against the alleles
carried by the most paths (`ref_*`); [the method](method.md#qc) defines every measure. qc skips a
pair whose optimal alignment would exceed `--max-cells` cells, the same test by which realign bands
an alignment with the centre, so it does not measure pairs as long as those that realign bands; the
column `skipped` counts them.

`qc.tsv.json` holds the totals over the contig, with precision, recall and F1, and `growth`, the
bases of the realigned regions' graphs over the bases of their longest alleles. `qc.tsv` has the
same counts per region, but not F1, to find the regions that pull a total down. A region's reference
F1 is 2 `ref_tp` / (`ref_opt` + `ref_msa`), and the regions where it is lowest are the first to look
at.

### Completeness

merge's checks compare the output only with its inputs, so they cannot detect a part that was left
out. The merged graph should hold every path of the whole-genome GBZ except those of the dropped
samples, given in `drop`, separated by `|`; `diff` prints nothing when the merged graph is complete:

```
drop='gref_CHM13'
vg paths -M -x hprc-v2.1.gbz | grep -v '^#' | cut -f1 | grep -Ev "^($drop)#" | sort > want.txt
vg paths -M -x hprc-v2.1.realigned.gbz | grep -v '^#' | cut -f1 | sort > got.txt
diff want.txt got.txt
```

A line marked `<` is a path the merged graph lacks: give merge the part that holds it. A line marked
`>` is a path that should not be there, such as a dropped sample's path in a `vg chunk` part. When
the cover has been rebuilt, its short paths need not match the original's in name or number, so
compare the rest: `grep -v '^gref_CHM13#' got.txt | diff want.txt -`.

## The gref_CHM13 cover

With the default `--drop-sample`, the merged graph has no gref_CHM13 cover. To give it one, rebuild
the cover in each realigned contig and each contig cut with prepare, before merge, on a mutable copy
of the graph, as the docstring of `prepare.py` describes. For a realigned contig:

```
c=chr1; g=out/$c/$c.realigned
vg convert -p $g.gbz > $g.pg
vg paths -x $g.pg -u -Q CHM13 > $g.gref.pg
vg convert -f $g.gref.pg > $g.gref.gfa
vg gbwt -G $g.gref.gfa --gbz-format -g $g.gref.gbz
mv $g.gref.gbz $g.gbz
rm $g.pg $g.gref.pg $g.gref.gfa
```

`vg paths -u` computes the cover from the paths that `-Q` selects, those whose names begin with
CHM13. The node ids do not change, so the rebuilt GBZ takes the realigned GBZ's place in the merge.
The contigs cut with prepare, such as chrM, need the same, with `g=parts/$c/$c`. Rebuild a contig's
cover only after its verify has passed, since verify compares with the prepared GBZ, which has no
cover, and would report the cover's paths as extra; when replace runs again for the contig, rebuild
the cover again. pgrealign does not run these commands, and checks nothing that they build.

## Checklist

1. Install pgrealign ([the README](../README.md#install)), and put vg and abPOA on `PATH`.
2. Choose the whole-genome GBZ, the reference sample, the targets and the samples to drop.
3. Set `TMPDIR` to a disk with room for merge's working files.
4. Split the reference sample's contigs into `contigs.txt` and `notargets.txt`.
5. Run `pgrealign run` for every contig in `contigs.txt`. List those whose run did not finish or
   whose verify failed, and run them again once the cause is fixed.
6. Read the realign statuses and qc totals.
7. Cut each contig in `notargets.txt` with `pgrealign prepare`, with the same sample options.
8. Cut the components with no reference-sample path with `vg chunk --gbz`, and remove from them any
   sample other than gref_CHM13 that the runs dropped.
9. If the cover is needed, rebuild it in each realigned contig and each prepared part.
10. Merge the realigned GBZs and the parts. If merge reports overlapping ids, find the cause
    ([the node-id rule](#the-node-id-rule)), move a contig that grew with `--id-start`, rebuild its
    cover if step 9 did, and merge again.
11. Compare the merged graph's path names with the whole-genome GBZ's.
12. Build the indexes that later tools need on the merged graph.
