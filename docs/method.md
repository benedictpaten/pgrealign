# The pgrealign method

pgrealign realigns chosen regions of a pangenome graph. For each region it gathers the sequence that
every haplotype spells there, builds one multiple alignment of those sequences, turns the alignment
into a small graph, and puts that graph in place of the region's old nodes, rewriting every
haplotype's walk through it.

This document describes each stage as the code implements it, with the reasoning behind each step.
It names the option that controls a step; `pgrealign <stage> --help` gives the default of most
options, and [reference.md](reference.md) gives every option, file format and summary field.
pgrealign runs two external programs, vg, the pangenome toolkit, and abPOA, a sequence aligner, and
`catalog` compiles a repeat finder of its own with a C compiler
([requirements](reference.md#installation-and-tools)).

## Graphs and walks

pgrealign works on a pangenome graph together with the haplotypes in it.
- A node is a stretch of sequence with an integer id. An edge joins two node ends.
- A walk is one haplotype's route through the graph: a list of node visits, called steps. A step
  visits its node forward, written `>12`, reading its sequence, or in reverse, written `<12`, reading
  its reverse complement. A haplotype assembled in pieces has one walk per piece.
- A generic path is a named route that is not a sample's haplotype. A path, in vg's sense, is either
  a walk or a generic path.

## Graph files

The graph is stored as a GBZ, vg's compressed file format for a graph together with its paths. vg
also writes a graph as GFA, a text format with one record per line: an S line holds a node, an L line
an edge, a W line a walk, and a P line a generic path or, in some graphs, the reference path (see The
reference); H lines hold the header. A W line names its walk by its sample, haplotype, contig and
start, as `sample#haplotype#contig#start`, and pgrealign takes that name, unchecked, as unique, as it
is in a GBZ.

## The reference

The reference sample, named by `--ref-sample`, gives the graph its coordinates. A contig is one of
its sequences, such as chr20, named as in the contig field of its walk. The contig's reference path
is the single W line of the reference sample on that contig, or a P line named `<sample>#0#<contig>`
or `<sample>#<contig>`; `prepare` stops unless there is exactly one. Reference coordinates are
0-based positions on the contig: a node's position on the reference path is the W line's start field
(0 for a P line) plus the length of the steps before it. A target, an interval of the reference to
realign, is given in the same coordinates, and under the same contig names.

## Snarls and chains

vg divides a graph into snarls. A snarl is a part of the graph cut off from the rest by two boundary
nodes: a walk can enter or leave it only through its boundary nodes. A trivial snarl is two boundary
nodes joined by an edge, with nothing between them. A chain is a series of snarls in which each
shares a boundary node with the next. A snarl can hold chains of its own: a chain inside a snarl is
nested, and a chain inside no snarl is top-level. Snarls and chains, nested by containment, form the
snarl tree.

## Regions and statuses

- A region is the part of the graph between two anchor nodes, L and R, chosen around one or more
  targets. Realigning a region replaces the nodes between its anchors that walks use, its interior
  (Extract defines it exactly). The anchors stay.
- A region gets a status at each of the three steps before it is replaced: when it is anchored,
  when its sequences are collected, and when they are aligned. A region that one of them cannot
  handle keeps its old subgraph, and its status says why. Only a region that is `ok` after all three
  is replaced.
- Every walk spells the same sequence after replacement as before, apart from case: realigned bases
  are upper case (see Known gaps). This is checked as each walk is rewritten, and again end to end
  for every path.

## The stages

Each stage is a command that reads and writes files, so it can be run on its own
([files](reference.md#the-stages-and-their-files)):

| command | run once per | does |
|---|---|---|
| `catalog` | genome | lists the reference's tandem repeats: the targets |
| `prepare` | contig | cuts the contig from the whole-genome GBZ, and indexes its snarls and reference path |
| `regions` | contig | anchors the targets, from the catalogue or any BED, in the contig's snarls |
| `extract` | contig | collects each region's sequences, and where every walk passes through it |
| `realign` | contig, or shard | builds one multiple alignment per region |
| `replace` | contig | puts each alignment's graph in place of its region, as GFA, or as GBZ with `--gbz` |
| `verify` | contig | checks that every path of the realigned GBZ spells what it did before |
| `qc` | contig | measures how good each alignment is |
| `merge` | genome | joins the per-contig GBZs into one |
| `union` | contig | merges the extracts of a panel graph and a recipient graph ([Realigning on a panel](#realigning-on-a-panel-union-and-project)), for `realign` |
| `project` | contig | cuts the union's alignments down to the recipient graph's members, for its `replace` |

`regions`, `extract` and `realign` are the three stages that give a region its status, and `replace`
turns each alignment into a graph as part of its own work. `pgrealign run` runs every stage from
`prepare` to `qc` for one contig, skipping any whose output exists
([reference.md](reference.md#run)), and a genome is built with `catalog`, a `run` for each contig
that has targets, the parts of the graph that are not realigned, cut out as they are, and `merge`
([whole-genome.md](whole-genome.md)). Only `catalog` knows that the targets are tandem repeats; the
later stages take any intervals, and the sequences in them, as they come.

## The catalogue

`catalog` lists the tandem repeats of a reference, from its RepeatMasker annotation
(`--repeatmasker`), its FASTA (`--reference`), or both, on the contigs that `--contigs` names. With
`--trf`, it also runs TRF on the FASTA, in place of a RepeatMasker annotation the reference lacks, and
with `--strf` it runs strf, pgrealign's reimplementation of that part of TRF, instead. The
chrom names of both sources must be the graph's contig names, because `regions` keeps only the
targets on its own contig. Each entry of the catalogue is a target, of one of four classes:

| class | meaning |
|---|---|
| STR | short tandem repeat: a motif of at most 6 bp, repeated |
| VNTR | variable-number tandem repeat: a motif of 7 bp or more, repeated |
| SAT | a satellite, as RepeatMasker annotates it |
| LC | a low-complexity stretch, as RepeatMasker annotates it |

### RepeatMasker records

`catalog` keeps the RepeatMasker records of three classes. A Simple_repeat record becomes VNTR when
its name has the form `(MOTIF)n` with a motif longer than 6 bp, and STR otherwise; Satellite records
become SAT, and Low_complexity records LC. A name of the form `(MOTIF)n` gives the record's motif and
period, the motif's length; any other name gives period 0, unknown, with the name as the motif.
RepeatMasker reports no Simple_repeat motif longer than 30 bp, so on its own it misses the VNTRs of
longer periods. The second source, trfind, finds them.

### trfind

trfind is a seed-based tandem-repeat finder in C (`pgrealign/c/trfind.c`). `catalog` compiles it on
first use and runs it on each FASTA file, on up to `-j` files at once. trfind scans each record of a
file on its own, upper-cased, so soft-masked (lower-case) bases count like any others.

1. Seeds. At each exact 10-mer of A, C, G and T, trfind takes the distance d to the nearest earlier
   copy of the same 10-mer as a candidate period, when d lies between `--min-period` and
   `--max-period`. A seed is such a 10-mer with its d. Only the nearest copy counts, so a 10-mer
   whose nearest copy is closer than `--min-period` gives no seed.
2. Chains. Seeds of about the same period form a chain. Each distance falls in a bin: d itself
   below 40 bp, and 40 + ⌊log(d / 40) / log 1.03⌋ from 40 bp on, so each bin above 40 bp is 3% wide.
   A bin holds at most one open chain. For each seed, trfind tries its own bin, then the bin below,
   then the bin above, and the seed joins the first of them whose open chain's last seed starts at
   most max(40, min(⌊d / 2⌋, 300)) bp before it, d being the seed's own distance. If none
   qualifies, the seed opens a new chain in its own bin, which closes the chain that bin held. Every
   open chain closes at the end of the record.
3. Qualifying. When a chain closes, its period p is its mean seed distance, rounded, and its span
   is the length from the start of its first seed to the end of its last. It is reported when it
   has at least 3 seeds, its span is at least 0.75 p, it has at least 0.15 seeds per bp of span, and
   its interval, from p bp before its first seed to the end of its last seed, clipped to the record,
   is at least 40 bp. Unclipped, the interval is the span plus a period, so it is then at least
   1.75 p long: about two copies.
4. Refining. A chain can lock onto a multiple of the true period. The reported period is the
   smallest divisor q of p whose self-identity over the interval at offset q is no more than 0.03
   below the self-identity at offset p. The self-identity at offset q is the fraction of the
   interval's positions i, among those for which i + q also lies in the interval, at which the base
   equals the base at i + q; N matches nothing.

A trfind repeat is VNTR when q is 7 or more and STR otherwise. Its motif is the first q bases of its
interval in the reference.

### TRF

RepeatMasker 4.1 finds its Simple_repeat records with TRF (Tandem Repeats Finder), and `--trf`, naming
a TRF binary, runs TRF the same way, for a reference that has no RepeatMasker annotation.

1. Batches. Each record of each FASTA file is cut into RepeatMasker's batches: windows of 60 kb, each
   starting 58 kb after the one before, so neighbours overlap by 2 kb, the last reaching the record's
   end. Besides following RepeatMasker, this bounds the work TRF does on any one repeat: run on the
   whole of CHM13's chr20, TRF had not finished after 48 minutes, and on its batches it took 10.
2. Stages. Every batch goes through RepeatMasker's two TRF stages, each with its own parameters
   (match, mismatch, indel, PM, PI, minimum score and maximum period) and its own copy number that a
   repeat must exceed to be kept: one for young repeats, `2 7 7 80 10 50 10` and 4 copies, and one
   for old repeats, `2 3 5 75 20 33 7` and 5 copies. `-j` runs that many TRF processes at once, each
   on its share of the batches.
3. Records. Each repeat kept is placed on its record by its batch's offset, and becomes a
   Simple_repeat record named `(MOTIF)n` by TRF's consensus motif, so STR or VNTR by the rule above,
   from the source `trf`. A repeat in the overlap of two batches can be found in both, and the two
   records merge.

RepeatMasker then rescores each repeat's alignment with its own matrix and drops the low-scoring
ones, and it masks interspersed repeats, such as Alu elements with their poly-A tails, before looking
for simple ones; `catalog` does neither. So TRF finds about every Simple_repeat that RepeatMasker would,
and more besides. On CHM13 chr20, TRF and trfind together find 98% of the entries of the shipped
catalogue, built from RepeatMasker's annotation and trfind (98% of its STRs, 99% of its VNTRs, within
200 bp), and write 1.9 times as many entries.

### strf

strf (`pgrealign/c/strf.c`) is pgrealign's own reimplementation of the part of TRF that `--trf` uses,
so that a reference without a RepeatMasker annotation needs no TRF, whose licence is the AGPL; strf,
like the rest of pgrealign, is MIT. `--strf` runs it in TRF's place, on the same batches and through
the same stages, and its records come from the source `strf`. It takes TRF's command line,
`strf File Match Mismatch Delta PM PI Minscore MaxPeriod`, and writes TRF's `-ngs` output, so the
steps above apply to it unchanged. It covers only what the two stages need:

- a PM of 80 or 75, the two values for which TRF has detection thresholds;
- periods up to 20, for which TRF allows no indels in detection, so that PI has no effect; strf
  ignores PI;
- the `-ngs` output.

It follows Benson (1999), "Tandem repeats finder: a program to analyze DNA sequences" (Nucleic Acids
Research 27:573-580), where the paper is specific, and TRF's observed behaviour where it is not. The
header of c/strf.c states the rules in full.

1. Detection. For each candidate period d up to MaxPeriod + 3, and at most 20, a position that ends
   a run of at least k bases each equal to the base d before it is a hit (k = 4 for PM 80 and 3 for
   PM 75, as in the paper's Table 1) when at least 5 of the last max(d, 20) bases lie in such runs,
   and the first such run in that window ends soon enough after its start: within a waiting time,
   computed from PM, of 18 bases for PM 80 and 15 for PM 75.
2. Alignment. The d bases ending at the hit are a pattern, aligned to the sequence around the hit by
   wraparound dynamic programming with the stage's match, mismatch and indel scores. A backward pass
   finds the left end, and a forward pass from there the right end and the alignment.
3. Tests. The alignment must hold at least 1.9 copies, and d must be among the three strongest peaks
   of the region's histogram of distances between identical 2-mers, once a linear trend is removed,
   so that a repeat of period p is analysed at p, 2p and 3p only. At d = 1, one base must instead
   make up at least 80% of the region.
4. Consensus. Each column of the pattern takes a majority vote, ties going to A, C, G and T in that
   order and a gap winning a tie with a base, which deletes the column; a base is inserted where at
   least half of the copies have an insertion. The region is aligned again to the consensus, and the
   repeat is kept when its score reaches the stage's minimum and it still holds 1.9 copies.
5. Period and redundancy. The reported period is the most common distance between matching bases of
   adjacent copies, which can differ from the consensus size. A repeat is dropped when its period
   exceeds MaxPeriod, and when at least 90% of it lies within another repeat of the same period and
   no lower score, or of a period that divides its own and a score at least 1/1.1 of its own.

strf was written without TRF's source code: from the paper, TRF's documentation, a description of
TRF's behaviour checked against the binary, and runs of the TRF binary, which settled each rule the
paper leaves open. It reads FASTA as TRF does, and treats N and every other letter but A, C, G and T as
TRF does; the header of c/strf.c says how. On the batches of CHM13 chr20, GRCh38 chr21 and yeast
S288C, it reproduces the start, end and period of 99.9% to 100% of TRF's records exactly in the first
stage and 99.1% to 99.3% in the second; the records it misses or adds are marginal ones, scoring near
the stage's minimum. With trfind, on chr20, it finds 98.0% of the shipped catalogue's STRs, 98.9% of
its VNTRs and 97.8% of all its entries, where TRF finds 98.2%, 99.1% and 98.0%. It runs both stages
on chr20's batches in 52 s of CPU, against 73 s for the fastest TRF tested, in 7 MB of memory.

### Merging records into entries

The records of each contig, from every source, are taken in order of start, then end. A record joins
the current entry when its start lies at most `--gap` bp past the entry's end so far, so records
that overlap or lie within `--gap` bp of each other become one entry.
- The entry's class is the class whose records cover the most bp of it. Ties go to VNTR, then SAT,
  then STR, then LC.
- Its period and motif are those of the longest record of that class, the first of equally long
  ones.
- Its id is `TR_<chrom>_<start>`.

The catalogue is a BED with a header line, one entry per line, that also records how many records of
each class and source merged into each entry ([format](reference.md#target-beds-and-the-catalogue)).
The later stages recognise a catalogue by its header line, and stop if its ids repeat. They read any
other BED as plain intervals, with class `.`. A line's id is its name when no other line has that
name, and otherwise `<name>_<chrom>_<start>`, or `R_<chrom>_<start>` without a name, with `_2`, `_3`,
... appended, in file order, where that would repeat an id.

## Preparing a contig

`prepare` cuts the connected component that holds the contig `--contig` names from the whole-genome
GBZ, with `vg chunk --contig`, and stops unless vg writes exactly one component. Node ids stay
genome-wide. When several contigs share a component, the contig's GBZ holds them all (see Known
gaps).

It then removes the samples that `--drop-sample` names, by default `gref_CHM13`, the gref cover of
CHM13: a copy of the reference path together with fragments covering the parts of the graph that the
reference path misses. With another reference sample, its cover, `gref_<sample>`, is removed only
when `--drop-sample` names it. The cover's fragments start and end inside regions rather than
crossing them, and each such fragment that ends inside a region adds a row to its alignment, and
often a pairwise alignment. The reference sample cannot be dropped, and the dropped samples are
absent from every later output.

Last, it writes the two inputs that `regions` needs: the contig's snarls, from
`vg snarls -T -P <sample>#0#<contig>`, trivial snarls included and the snarl tree oriented along the
paths whose names begin `<sample>#0#<contig>` (see Known gaps); and the reference index, each node
of the reference path in order, with its orientation and its start and end in reference
coordinates. `prepare` skips each output that exists.

## Regions

`regions` takes the targets whose chrom is the contig that `--contig` names, turns each into a region
by choosing its two anchor nodes, and merges regions that would interfere. It writes `regions.tsv`,
one row per region in reference order, with its anchors, its interior on the reference, from L's end
to R's start, its targets and its status ([format](reference.md#regionstsv)).

### Anchors

The candidate anchors are the snarl boundary nodes that the reference path visits exactly once, and
forward. Each belongs to one chain, and has one position on the reference, so every region's forward
orientation is the reference's. `regions` builds the chains from every snarl in the snarls file,
top-level and nested: snarls that share a boundary node belong to one chain. A chain with fewer than
two candidates is not used.

For each target:
1. The target is padded by `--pad` bp on each side, within the contig.
2. The padded interval is clamped to lie between the end of the first candidate on the reference,
   of any chain, and the start of the last, so a target that reaches into either is anchored on it.
   An interval that ends at or before the first candidate's end, or starts at or after the last
   one's start, has no candidates on both sides, and its target is `too_long`.
3. In each chain, L is the last candidate that ends at or before the padded start, and R is the
   first that starts at or after the padded end.
4. Of the chains that give both, the one whose pair spans the least reference, from L's start to R's
   end, wins, provided that span is at most `--max-span`. Between chains of equal span, which one
   wins is not defined.

Both anchors come from one chain because two boundary nodes of one chain separate the graph: the
stretch between them is a series of whole snarls, each entered and left only through its boundary
nodes. So a walk that crosses the region visits both anchors, and a walk that visits only one either
starts or ends inside the region, or turns back inside it, which `extract` marks complex. Boundary
nodes of two different chains carry no such guarantee. With padding, each anchor lies, where the
reference allows, at least `--pad` bp outside the target, so the region holds some flank on each
side.

### Merging

Regions are taken in order of their left anchor's start. A region merges into the one before it when
its left anchor starts before that region's right anchor starts. This catches both interiors that
overlap and an anchor that lies inside another region's interior, where replacing one region would
change the other. The merged targets are anchored again as one, from the smallest of their padded
starts to the largest of their padded ends, and the pass repeats until nothing merges.

Two regions that only share a boundary node, one's R being the next one's L, stay separate. Each
replaces only its own interior, and the shared node stays in the graph.

A region whose anchors cannot lie within `--max-span` of each other gets the status `too_long`, so
that later stages leave it as it is; every other region is `ok`. A target that no chain anchors
within `--max-span` keeps its own unpadded interval, with `.` for both anchors, and takes no part in
merging, so another region can overlap it. Merged targets that cannot be anchored again as one keep
the outermost anchors of the regions that merged.

## Extract

`extract` reads the contig's graph as a stream of GFA lines, from `vg convert -f` or from a GFA file,
plain or gzipped, that lists every S and L line before its first W line, as vg writes it. It records,
for every region that `regions` made `ok`, how each walk passes through it. One pass over the graph
serves every region. It reads walks only; generic paths take no part.

### Runs

A run is a stretch of consecutive steps of one walk inside one region. The order of a walk's visits
to the region's two anchors decides its runs:

| run | anchor visits along the walk | the run's steps |
|---|---|---|
| spanning, forward | `>L`, then `>R` | those between the two visits |
| spanning, reverse | `<R`, then `<L` | those between the two visits |
| prefix | the walk starts inside, and its first anchor visit leaves (`>R` or `<L`) | from the walk's start to that visit |
| suffix | the walk enters (`>L` or `<R`) and ends inside | from after the entry to the walk's end |
| internal | the walk visits no anchor of any region (see the next section) | the whole walk |
| complex | any other order, such as a U-turn or a second entry before leaving | none |

A walk can pass through one region more than once, and then has one run for each pass. A spanning
run is empty when L and R are adjacent in the walk; a prefix or a suffix never is, since a walk that
starts or ends on an anchor has no run at that end. A run's `orient` is `+` when its walk runs in the
region's forward orientation and `-` when it runs the other way, and its sequence is spelled in the
forward orientation, a reverse run's steps reversed and each flipped. An internal run has no anchor
visit to orient it, and takes its `orient` from the region's other runs (see Orientation).

### Walks that visit no anchor

A walk that visits no anchor of any region, an anchor-free walk, may lie wholly inside one region. A
region's interior nodes, which `replace` removes, are the nodes that its runs visit, internal runs
included. Which region an anchor-free walk lies in is known only once the pass has found every other
run, so `extract` holds a short one, which spells at most `--max-internal-bp`, whole until the end,
and a long one only as the nodes it visits, which bounds the memory the pass needs. At the end:
- a short walk that visits nodes of the spanning, prefix or suffix runs of exactly one region
  becomes an internal run of that region, and its nodes join the interior;
- a short walk that visits such nodes of several regions makes each of them `complex`;
- a long walk that visits a region's interior nodes makes that region `complex`. Such a walk is
  never recorded as a run, and `replace` refuses to rewrite a walk through a replaced interior
  outside a recorded run;
- an anchor-free walk that visits no interior node is left alone.

### Orientation

Turned to the region's forward orientation, as each run's sequence is, every run should visit each
interior node the same way. A node can be stored reversed, as when an insertion was first added to
the graph from a contig on the reverse strand: every run then visits it `<`, and spells the right
sequence through it. But a node that runs visit both ways marks an inversion between haplotypes: the
reference crosses its nodes forward and an inverted haplotype crosses them backward, or one walk
turns back on itself, as in a hairpin. `realign` would spell the inverted allele as a plain forward
sequence and align it to the centre forward, losing the reversal, so such a region is made `complex`
and keeps its subgraph, the reversed walks included. Each run is checked against the ways earlier
runs, and its own earlier steps, visited its nodes.

An internal run is oriented by the nodes it shares with the region's spanning, prefix and suffix
runs. When it visits each of them the way those runs do, it is `+`. When it visits each the other
way, it runs against the region: it is `-`, and spelled reversed and flipped, as a reverse run is,
so that its fragment can be found in the alleles. When it agrees on some nodes and not on others, it
runs both ways through the region, and the region is made `complex`.

### Alleles and fragments

The runs' sequences are collected into the region's members: its alleles and its fragments.
- An allele is a distinct sequence that spanning runs spell. Its weight is the number of spanning
  runs that spell it. Alleles are numbered `a1`, `a2`, ... in the order the pass first meets them.
- A fragment is a distinct sequence of a prefix, suffix or internal run, together with its side: the
  end of the region it touches, in the forward orientation. `left` touches L (a forward suffix, or a
  reverse prefix), `right` touches R (a forward prefix, or a reverse suffix), and `internal` (an
  internal run) touches neither. Fragments are numbered `f1`, `f2`, ..., and carry weights in the
  same way.

Runs that spell the same member share it, so `realign` aligns each distinct sequence once. Fragments
exist because some walks start or end inside a region, for example where a haplotype's assembly
breaks; they become rows of the alignment too, so that their walks can be rewritten like any other.
The reference allele, which only `qc` uses, is the allele of the last spanning run of a
reference-sample walk; a region with no such run has none.

### Statuses

Each region gets the first status that applies:

| status | when |
|---|---|
| `complex` | a walk has a complex run, runs visit an interior node both ways, or an anchor-free walk makes the region complex as above |
| `no_runs` | no walk has a run in the region |
| `invariant` | there is one allele and no fragment, so there is nothing to realign |
| `ok` | any other region: two or more alleles, or at least one fragment |

A reference path that is a W line crosses every region `regions` made `ok`, so `no_runs` means that
the graph given to `extract` has a P-line reference path or is not the graph the regions were
anchored in. `extract` writes a package per region, holding its status, members and reference
allele; every run that is not complex, which `replace` rewrites; and each region's interior nodes,
which `replace` removes ([packages](reference.md#packagesjsonlgz),
[runs](reference.md#occurrencestsvgz), [interiors](reference.md#interiorstsvgz)).

## Realign: the star

`realign` builds one multiple alignment for each `ok` region that it selects (see Which regions are
realigned). A multiple alignment here is a set of rows of equal length, one per member, each row
being the member's sequence with gap characters (`-`) inserted. A column is one position across all
the rows.

The method is a star alignment around one allele, called the centre: the region's reference allele,
where it has one, or else a medoid, the allele with the least weighted distance to the others:
1. Centre. The reference allele, or the medoid of the alleles (see The centre).
2. Pairwise alignment. Every other allele is aligned to the centre alone.
3. Left normalisation. Each pairwise alignment's indels are moved as far left as they can go without
   lowering its score.
4. Merge. The pairwise alignments are merged into one, with the centre's bases as columns.
5. Fragments. Each fragment either takes the row of an allele that contains it, or is aligned to the
   centre like an allele.

The method lives in `star.py`, which knows nothing of graphs or files and takes the aligner as a
function ([the aligner contract](reference.md#the-aligner-contract)). `realign.py` applies it to
extract's packages.

### The aligner

Every alignment is computed by abPOA, a partial-order aligner, under the scores below, run as
`abpoa -m 0 -r 1`: globally, from end to end of every sequence, writing the alignment's rows. To
align two sequences abPOA fills a dynamic-programming matrix of one cell per pair of positions, one
from each sequence, at 20 bytes per cell. A band limits it to the cells near the alignment's
path. abPOA's own band is adaptive, reaching about 10 cells plus 1% of the sequence's length either
side of the path, and `-b -1` turns it off. A band saves time, but less memory than its cells
suggest: abPOA reserves the whole matrix either way, rounded up to a power of two, and a banded
alignment still touches at least one page per base of the centre in each of its five score matrices
([limits](reference.md#limits)). The aligner is told, for each alignment, whether a band
is wanted: alignments to the centre are unbanded unless too large, and alignments of inserted
sequences to one another are banded, as the sections below say.

abPOA is given A, C, G, T and N only. Every sequence is upper-cased first, and any other character is
given to abPOA as N and restored in the returned rows. Every row is checked to spell its member's
upper-cased sequence, and all rows to have one length. The centre is chosen from the same sequences.

### The centre

The centre is the region's reference allele, the package's `ref`, when it has a positive weight and no
more than `--max-excess`, 0.1, of the haplotypes' bases lie beyond its length (see below); and the
medoid of the alleles otherwise, as when the reference does not span the region (`--centre medoid`
takes the medoid always, and `--max-excess 0` the reference allele whatever). Every other allele is aligned to the centre on its own, so around
the reference, the alignment that a genotyper reads between the reference and any allele, the pair's
projected alignment (see Merging into one alignment), is that allele's own pairwise alignment to the
reference: optimal under the scores, unless banded, and left-normalised against the reference, as a
variant caller writes its records. Around any other centre it is the composition of two alignments
to the centre, the reference's and the allele's, which is neither: an indel that each places as far
left as it can go on the centre can land anywhere on the reference, so the same difference from the
reference, carried by several alleles, can come out as several records at different places.

Measured with `qc` on regions of the HPRC graph of 30 haplotypes realigned on CHM13 (a sample of
1,468 regions of chr2, chr6, chr14 and chr16 with three or more alleles, one of 100 bp or more,
aligned around each centre in turn), the reference as centre gives:

| centre | ref F1 | pair F1 | growth | gap runs of 50+ bp against the reference | seconds |
|---|---|---|---|---|---|
| medoid | 0.9934 | 0.9949 | 1.030 | 830 | 1,169 |
| reference allele | 1.0000 | 0.9942 | 1.027 | 625, as many as the optimal alignments | 981 |

The reference's own alignments become optimal (ref F1 1, by construction), at a small cost to pairs
of other alleles, which are now composed through the reference rather than through the medoid
(pair F1), and the graph does not grow. The 726 of these regions whose medoid is not the reference
allele carry the whole difference. ([QC](#qc) defines the metrics.)

That cost is larger where many alleles are much longer than the reference, as at a tandem repeat
expanded in most haplotypes. Each allele's extra copies are inserted into the reference where its own
alignment puts them, and two alleles' insertions are aligned to each other only where they fall in the
same slot (see Merging into one alignment). At chr2:232.31M (30 alleles of up to 48.8 kb around a
27.7 kb reference), pairs of other alleles keep 213,917 of the 499,974 aligned pairs of bases of their
optimal alignments, against 416,142 around the medoid, while the region's variant sites, runs of
columns whose rows differ, fall from 1,672 to 1,068. And an allele's optimal alignment to the
reference is not always the one with the fewest records: at chr2:90.32M, where alleles hold diverged
copies of the repeat, it aligns them base against base, and the region's 29 other alleles differ from
the reference in 20,391 runs of columns (16,986 of mismatches only), against 7,676 (6,588) for the
alignments composed through the medoid, which leave more of the copies unaligned, as long indels, in
an alignment of 64,487 columns against 39,315.

Where most haplotypes are much longer than the reference allele, as where they carry an expansion it
lacks, the cost becomes the graph's size: their extra sequence is all inserted around the reference
and aligned only slot by slot, and much of it stays unaligned. Joined both ways, the 30-way graph's
realigned regions held 2.8 Mb more sequence around the reference allele than around the medoid, and
228 regions whose graph the reference made at least half as large again carried 5.5 Mb, all with a
reference allele far shorter than most haplotypes': at chr19:8.79M, a 23 kb reference allele against
75 kb at the median, 717,811 bases of graph against 129,260. `excess`, the share of the haplotypes'
bases that lie beyond the reference allele's length, each allele counted once per haplotype, picks
these regions out. Taking the medoid where it exceeds `--max-excess`:

| `--max-excess` | regions on the medoid | bases of graph against the medoid everywhere | blow-ups caught |
|---|---|---|---|
| none | 0 | +2.81 Mb | 0 of 228 |
| 0.25 | 980, 0.11% | -2.59 Mb | 150 |
| 0.15 | 1,683, 0.19% | -2.86 Mb | 201 |
| 0.1 | 2,508, 0.28% | -2.89 Mb | 220 |
| 0.05 | 4,592, 0.51% | -2.74 Mb | 227 |

At 0.1, the default, the realigned graph is smallest, and the eight blow-ups left carry 32 kb. The
regions at chr2:232.31M and chr2:90.32M above keep the reference allele, their excess 0.095 and 0.018.

### The medoid

The distance between two alleles comes from their 15-mers. Each allele's k-mers, with k = 15, form a
multiset: a k-mer that occurs three times counts three times. k-mers that contain a character other
than A, C, G or T are skipped. The distance between alleles i and j is one minus the Jaccard index
of their multisets:

    J(i, j) = shared / (|i| + |j| - shared)

where `shared` sums, over every k-mer, the smaller of its two copy numbers, and |i| is the number of
k-mers in i. When neither allele has a k-mer, as when both are shorter than 15 bp, J is 1, so the two
are at distance 0; an allele without k-mers is at distance 1 from one that has them. The medoid is
the allele i of positive weight with the least weighted sum of distances to the others, the sum over
every other allele j of w_j (1 - J(i, j)), where w_j is j's weight. Ties go to the earlier allele.
If no allele has a positive weight, every weighted sum is 0 and the medoid is the first allele.

The multiset, rather than the set, carries copy number: tandem-repeat alleles differ mostly in how
many copies of the motif they hold, and two that differ only in copy number often contain the same
distinct k-mers. Because the distances are weighted, the medoid lies near the sequences that most
runs spell. An allele can have weight 0 when `union` adds it (see Realigning on a panel).

### Estimating distances from sketches

Comparing one pair exactly costs about as many steps as the smaller allele has distinct k-mers. When
the sum over all pairs exceeds `star.MEDOID_WORK`, 10^9 steps, every distance is estimated from
sketches instead, which makes each pair cost about the same however long the alleles are. An
allele's sketch is the `star.SKETCH`, 1,000, smallest hash values of its k-mer multiset, each copy
of a k-mer hashed as an element of its own, so that copy number still counts. J(i, j) is estimated
as the fraction of the 1,000 smallest values of the two sketches' union that lie in both sketches;
two empty sketches are at distance 1 (see Known gaps). The hash is a fixed function, so the
sketches, and the medoid, are the same in every run.

### Homology

The star assumes that every allele is a variant of the centre. An allele that is unrelated to it, as
at a non-homologous structural variant, or reversed, as at an inversion, has no good forward
alignment to it, but the aligner returns one anyway: the allele threaded through the centre by
hundreds of short chance matches, each of which becomes a variant of its own in the graph, and a VCF
record for every haplotype that carries the allele. So once the centre is known, and before anything
is aligned, `realign` compares each allele of positive weight with it, using their 15-mers (those
that choose a medoid), as sets:

    share = shared distinct k-mers / the distinct k-mers of whichever of the two has fewer

so an allele that is the centre plus an insertion, or minus a deletion, shares them all. If any
allele shares less than `--min-share`, 0.5, the region is `not_homologous`, and keeps its subgraph.
An allele shorter than 100 bp is not judged, nor is any allele when the centre is that short: one
variant can remove most of a short allele's k-mers, and a short allele cannot be threaded into many
variants. Fragments are not judged. The status's message names the allele that shares least, its
weight, and what it shares reverse-complemented, which is high for an inversion.

Measured on the HPRC graph realigned on CHM13 (1,122,539 realigned regions, 915.7 Mb of reference),
against the change in each region's non-reference VCF sample-records that realignment made, on chrX
and chr15:

| `--min-share` | regions refused | span | chrX: added records taken back / removed records given up | chr15: the same |
|---|---|---|---|---|
| 0.3 | 1,041 | 3.5 Mb | 402,166 / 30 | 55,355 / 2,242 |
| 0.4 | 1,761 | 4.3 Mb | 409,282 / 1,232 | 90,409 / 3,044 |
| 0.5 | 2,843 | 5.6 Mb | 410,601 / 2,904 | 159,993 / 6,396 |
| 0.6 | 4,267 | 7.4 Mb | 425,028 / 5,135 | 164,057 / 9,812 |
| 0.7 | 6,762 | 10.0 Mb | 429,317 / 11,936 | 172,202 / 12,367 |

These count every allele; judging only alleles and centres of 100 bp or more refuses 2,794 regions at
0.5 with the same chrX figures. At 0.5 the refused regions include the eight regions at four chrX
loci of inversion-polymorphic segmental duplications (MAGED4B, MAGEA, H2BW and HSFX), the five of
chr15's acrocentric p-arm cluster at 14.2-14.4 Mb and the one at GOLGA6L1, where realignment had
added most of the records, and few others: the other regions refused on chrX took back 55,227 records
and gave up 2,904. A lower threshold misses most of the chr15 regions, and a higher one gives up more than it
takes back. Counting events, the variants an allele has against the centre once aligned, was tried
as well and not adopted: among the regions the share passes, every event threshold tried gave up
about as many records on chrX as it took back.

These were measured with the medoid as centre. Around the reference allele the guard compares each
allele with the reference instead, and refuses nearly the same regions: on nine chromosomes of the
same graph (chr1, 2, 4, 5, 7, 11, 15, 22 and X), 1,401 regions against 1,385, 1,364 of them the same,
among them every region of the chrX and chr15 loci above. On chrX the regions it refuses take back
410,605 records and give up 3,470 (410,601 and 2,904 around the medoid), and on chr15 159,764 and
7,466 (159,142 and 6,396), so 0.5 stays.

### Pairwise alignment to the centre

Every allele is aligned to the centre with abPOA, so the number of alignments grows linearly with
the number of alleles. The alignment is unbanded, and so optimal under the scores, unless its
matrix, the centre's length times the allele's, would exceed `--max-cells` cells. Such an alignment
is banded instead, which bounds its memory and time but finds only the best alignment within the
band; the output's `banded_pairs` counts these, fragments aligned to the centre included. Two cases
need no aligner: a sequence equal to the centre, once both are upper-cased and other characters made
N, aligns to it base for base, and an empty allele, or an empty centre, aligns as one run of gaps
opposite the other sequence. The rows abPOA returns must spell their two inputs, or the region fails
as `aligner_failed`.

### Left normalisation

An indel inside a repeat can sit at any of several positions with the same score. The merge places
each insertion by its position on the centre, so two alleles that carry the same insertion share
columns only if both alignments put it in the same place. Every pairwise alignment is therefore
left-normalised. A run of gaps in one row moves one column left when the column before it holds two
bases and the other row's base in that column equals the other row's base in the run's last column.
The column the run leaves then holds the pair of bases that the column before it held, so the score
is unchanged, unless the run comes to abut another gap in its row, which joins the two and can only
raise the score; that happens only in a banded alignment, which need not be optimal. Moves repeat,
in both rows, until none applies. For example, the alignment on the left becomes the one on the
right:

    GCACACAT    GCACACAT
    GCACA--T    G--CACAT

### Merging into one alignment

Each pairwise alignment is read against the centre, whose bases are numbered 0 to n - 1:
- for each centre base, the character opposite it: a base, or a gap;
- for each insertion slot, the bases inserted there. Slot q lies just before centre base q, and slot
  n lies after the last centre base, so there are n + 1 slots.

In each slot, the distinct inserted sequences of every allele and every fragment aligned to the
centre are collected. A single sequence is used as it is. Several are aligned to each other by
abPOA, banded, given longest first, and those of equal length in alphabetical order; the output's
`slot_alignments` counts these alignments. Each row is then slot 0, centre column 0, slot 1, and so
on to slot n. A member with nothing inserted in a slot has gaps across that slot's width, and
columns that are gaps in every row are dropped. Members that insert the same sequence in one slot
take the same line of that slot's alignment, so in the graph they share its nodes. A slot
alignment's lines are not checked as they come back; the final check on every row catches a wrong
one.

A multiple alignment implies an alignment of every pair of its rows, the pair's projected
alignment: the columns in which either row has a base. A pair of the centre and an allele projects
to exactly the allele's pairwise alignment. Two alleles neither of which is the centre were never
aligned to each other: their projected alignment follows from their alignments to the centre, and
`qc` measures how far it falls from their own optimal alignment.

### Fragments

Each fragment is first looked for in the alleles, in their order: a `left` fragment as a prefix of
the allele, a `right` fragment as a suffix, and an `internal` fragment anywhere, at its leftmost
occurrence. A fragment found in an allele takes that allele's merged row, with every base outside
the matching stretch replaced by a gap, so in the graph it adds no sequence: its walk follows the
allele's, though nodes may be cut where it starts and ends (`fragments_exact`). A fragment found in
no allele is aligned to the centre like an allele, banded under the same rule, and joins the merge
(`fragments_aligned`).

### Scores

abPOA scores an alignment with a match score M, a mismatch penalty X, and two-piece gap penalties: a
gap of length n costs min(O1 + n E1, O2 + n E2). The first piece is the cheaper for short gaps, the
second for long ones. pgrealign's scores are abPOA's defaults times four, except the long-gap
extension E2, which stays 1:

| score | abPOA's default | pgrealign |
|---|---|---|
| match M | 2 | 8 |
| mismatch X | 4 | 16 |
| gap open O1, O2 | 4, 24 | 16, 96 |
| gap extension E1, E2 | 2, 1 | 8, 1 |

Under these scores a long gap costs 1 per base against a match's 8, a quarter of its relative cost
under the defaults, and the long piece is the cheaper one from 12 bases on, against 21. Under the
defaults, two unrelated stretches are cheaper to align base against base, threading short gaps
through chance matches, than to delete one and insert the other, so one replacement becomes many
small bubbles in the graph; under these, the deletion and insertion score better. abPOA takes
whole-number scores, and its long-gap extension is already 1, so it can become cheaper relative to
the rest only by scaling the rest up, by one factor so that their ratios to each other stay as
abPOA's defaults set them.

`--scores` takes other scores as abPOA's own flags; a flag left out takes abPOA's default, not the
value above ([syntax](reference.md#the-command-line)). `qc` should be given the scores the alignment
was built with; `run` passes its `--scores` to both.

### Limits, and regions that keep their old subgraph

A region the star cannot be built for gets one of these statuses, and `replace` leaves its subgraph
as it is:

| status | when |
|---|---|
| `too_big` | the region exceeds one of the limits on its size below |
| `timeout` | only with `--timeout`: the region ran longer than that many seconds |
| `aligner_failed` | abPOA exited with an error, left out a row, or returned pairwise rows that do not spell their inputs |
| `no_alleles` | the region has fragments but no allele: no walk spans it |
| `not_homologous` | an allele shares less than `--min-share` of its 15-mers with the centre (see Homology) |
| `error` | any other failure, recorded with its message: a bug, a slot alignment whose rows the final check rejects, or abPOA could not be started |

Whether a region is aligned is decided from its package alone, before anything is aligned, so the
same packages give the same alignments, and the same graph, on any machine. Nothing that depends on
the machine decides it: not the time a region takes, which depends on the machine's speed and on
what else runs there, nor the memory `realign` is given, nor the machine's page size. The limits on
a region's size, each of which makes it `too_big`, are all computed from the lengths of its alleles
and fragments ([reference.md](reference.md#limits)):
- even a banded alignment of the region's longest allele or fragment would exceed `--max-cells`
  cells. A banded alignment of a sequence of length n to one of the same length is taken to fill
  the smaller of n² and 2n (10 + ⌊n / 100⌋) cells, abPOA's adaptive band. For long sequences this
  is about n² / 50, so the limit admits sequences about seven times as long as an unbanded
  alignment of two equal lengths would;
- its alleles hold more than 40,000,000 bases in all, since choosing the centre holds every allele's
  k-mer counts at once;
- its work exceeds `--max-work`: the cells of dynamic-programming matrices its alignments to the
  centre fill, a banded one's band spanning the two lengths' difference, and its slots', taking each
  member's excess over the centre to fall in one slot and be aligned to the longest;
- its memory estimate exceeds `--max-region-memory` GB ([reference.md](reference.md#memory));
- abPOA would reserve more than `--max-address-space` GB for one of its alignments, which abPOA does
  for the whole matrix, rounded up to a power of two, banded or not. Linux refuses a reservation
  larger than the machine's memory and swap, so without this limit such a region would be aligned
  on a large machine and fail on a small one.

Run time follows the work. On the 7,288 regions of chrX in the HPRC graph realigned on CHM13 that
took a second or more, the logarithms of the two correlate at 0.90, and the 160 that took a minute
or more ran at 3, 6 and 9 × 10^7 cells a second (the 10th, 50th and 90th percentiles), so the
default `--max-work`, 2 × 10^11 cells, is about an hour of one core.

On that graph's 1,123,110 regions, the defaults leave 346 `too_big`: 155 for `--max-cells` or their
total bases, as before, 40 for their work, 88 for their memory estimate and 63 for abPOA's
reservation. They admit 298 of the 344 regions that the earlier time limit, 900 s, left as they were,
and leave out 106 that it aligned, on machines of hundreds of GB. For such machines,
`--max-region-memory 128 --max-address-space 256` leaves out 24 of those and admits 320 of the 344.
The largest region the defaults admit on any chromosome needs 62.5 GB, by its estimate.

`--max-memory` decides only how many regions run at once: `realign` starts them in order while the
estimates of those running fit within it. If a region within `--max-region-memory` needs more than
`--max-memory` on its own, `realign` stops before it aligns anything, with an error that names the
region and how much it needs, rather than leave it out; `--plan` gives that figure, and the regions
each limit leaves `too_big`, without aligning anything. Nor does the machine decide: `realign` stops
before it aligns anything if a region it admits would need a larger reservation than the machine
grants, and a region the machine fails all the same (abPOA's reservation refused, abPOA killed, the
worker out of memory) is not recorded: `realign` aligns the others and then stops, and run again where
the region can be aligned, it resumes with it. A time limit, `--timeout`, applies
only when given: it is checked while choosing a medoid and before each call to abPOA, and abPOA
itself is stopped at the deadline. Every other failure is recorded in the region's own result, so
one region cannot stop a contig.

### Which regions are realigned

`realign` aligns the `ok` packages that `--select` selects: `any` selects every `ok` package, and
`length` only those with fragments or with alleles of different lengths. A package that is not
selected gets no result, and keeps its subgraph. `--shard` splits the packages between separate
jobs, and an interrupted run resumes from the results it finished
([reference.md](reference.md#parallel-runs-and-resuming)). The output, `msas.jsonl.gz`, holds one
result per aligned region, in package order: its status, and for an `ok` region the centre, the rows
and the counts named above ([format](reference.md#msasjsonlgz)).

## Replace

`replace` writes the contig's graph with each realigned region's interior replaced. It reads the
graph as `extract` does and rewrites it in one pass; under dense numbering, an earlier pass reads
every node id (see Numbering the nodes). It replaces a region when the region's package is `ok` and
its alignment is `ok`. Several `msas.jsonl.gz` files, one per shard, can be given, and a region
found in two of them is an error.

### The graph of an alignment

For each region it replaces, `replace` first turns the alignment into a graph in which each row is a
walk, the induced graph (`induce.py`).

1. Columns. Each column gets one node for each distinct base in it. A row with a base in the column
   steps to that base's node, and an edge joins the row's previous node to it. A gap in a row adds
   nothing to its walk.
2. Unchopping. A node u merges into a node v when v is u's only successor, u is v's only
   predecessor, no walk ends at u and none starts at v. Every walk that visits either then visits
   both, one after the other, so merging them loses nothing.
3. Ordering. Every edge runs from one column to a later one, so the graph has no cycle. The nodes are
   numbered 1 to n in topological order: at each iteration, of the nodes whose predecessors are all
   numbered, the next is the one that begins with the earliest node of step 1, which creates nodes
   column by column, and a column's nodes in the order in which the rows, taken in member-id order
   as strings (`a10` before `a2`), first hold their bases.
4. Chopping. A node longer than 1,024 bases, vg's default maximum node length, is cut into pieces of
   at most 1,024 bases, numbered in turn, so the ids stay topological. `vg gbwt` then has no node to
   cut when it builds a GBZ, so the GBZ keeps the GFA's node ids.

Every edge joins the end of one node to the start of another. The step checks that every row's walk
spells the row without its gaps. An empty row, an allele that deletes the whole interior, has an
empty walk.

### Numbering the nodes

`--id-mode` chooses one of two numberings ([reference.md](reference.md#node-ids)).
- Dense numbering, the default, renumbers every node, from the contig's smallest id, or from
  `--id-start`. The nodes that stay keep their order, each region's new nodes come right after its
  left anchor, in topological order, and the ids run without gaps. New nodes then sit next to their
  neighbours in id order, so tools that read a graph in windows of node ids, such as vg call, find
  them in the same window.
- Reuse numbering keeps every old node's id. A region's new nodes take its freed ids, in increasing
  order, and then ids counted up from `--id-start`, which must be free in the whole genome. New
  nodes can land far from their regions in id order.

`merge` needs each input's id range, from its smallest id to its largest, to overlap no other
input's; [whole-genome.md](whole-genome.md#the-node-id-rule) says how to renumber a contig whose ids
grew into the next one's.

### What changes, and how it is checked

The anchors L and R stay. Inside:
- Nodes. The interior nodes are removed, and the induced graph's nodes are added.
- Edges. Edges that touch a removed node are dropped. The induced graph's edges are added, and so
  are the edges that join it to the anchors: L to the first node of each allele's walk and of each
  left fragment's, the last node of each allele's walk and of each right fragment's to R, and, for
  an empty allele, L to R unless the input already has that edge.
- Walks. Every run that `extract` recorded is replaced by its member's new walk, forward when the
  run's orient is `+`, and reversed and flipped when it is `-`. Steps outside every run are copied,
  renumbered under dense numbering.
- Generic paths. P lines are copied in the same way; one that visits a removed node is an error (see
  Known gaps).

H lines keep the reference-samples tag but lose vg's graph-name tags, which would make `vg gbwt -G`
record the output as a subgraph of the input. The output is a GFA; with `--gbz` it is turned into a
GBZ with `vg gbwt -G`, and the GFA is removed unless `--keep-gfa` is given.

`replace` checks, as it writes, that each rewritten run's old steps, spelled from the input's nodes,
and its new walk, spelled from the induced graph, both equal its member's upper-cased sequence. It
also checks that no walk visits a removed node outside a recorded run, that a walk's runs do not
overlap, that every node a run visits is in a replaced interior, that no anchor of one region lies
in another region's interior, that every S and L line comes before the first W line, so that the old
runs can be spelled, and that every new node's id is free in the contig's graph; ids elsewhere in
the genome are checked only by `merge`, which compares ranges. A failed check stops the stage and
leaves no output. It means that the inputs are wrong, not that one region is hard, so `replace` does
not skip the region and carry on.

## Verify

`verify` compares the contig's graph before `replace`, from `prepare`, with the realigned GBZ. Runs
through a replaced region have new steps, and dense numbering can shift ids, so it compares what
the paths spell and how they are named rather than their steps. The path metadata (`vg paths -M`)
must be identical, and every path's sequence (`vg paths -F`), upper-cased, must have the same length
and MD5 digest in both graphs. The report counts the paths missing from the realigned graph, extra
in it, and in both with different sequences ([format](reference.md#verifys-report)). It is `ok` when
the metadata is identical and all three counts are 0; otherwise `pgrealign verify` exits with
status 1, and `pgrealign run` stops before `qc`. Because it reads the GBZ built from `replace`'s
output, it also covers the conversion to GBZ.

## QC

`qc` measures, region by region, how far the projected alignments of pairs of alleles (see Merging
into one alignment) fall from the pairs' own optimal alignments, how much sequence the alignment
puts in the graph, and how much unrelated sequence it forces into shared columns. Given the graph
`extract` read, with `--graph`, it measures that graph's own alignment of the same pairs against the
same optimal alignments, so each region has a before and an after; `run` always gives it.

### What is compared

`qc` reads extract's output directory and one realign output. It takes the regions that `realign`
made `ok` with more than one row, or with `--sample` a fixed random sample of them, and measures
those whose package has more than one allele. Within a region it compares two sets of allele pairs:
- random pairs: up to `--pairs` pairs, drawn without replacement from all pairs of alleles other
  than the reference allele, the centre included, with a seed fixed by the region's id, so a rerun
  compares the same pairs;
- reference pairs: the reference allele against each of up to `--refs` other alleles, the heaviest
  first, ties broken by id in string order, so `a10` comes before `a2`. A genotyper's records,
  written against the reference, follow these alignments. A region with no reference allele has no
  reference pairs, and its `ref_*` fields are 0.

Each pair's projected alignment is compared with its optimal alignment: the pair's global,
unbanded, left-normalised abPOA alignment under the scores given with `--scores`. A pair whose
optimal alignment would need more than `--max-cells` cells, the product of the two lengths, is
compared with abPOA's banded alignment instead, which can fall short of the optimum, and is counted
in `banded`; `run` gives `qc` the same `--max-cells` as `realign`, so a pair through the centre that
`realign` banded is compared with the very alignment `realign` used. A pair whose alignment's
memory estimate alone exceeds `--max-memory` is not compared at all, and is counted in `skipped`
([reference.md](reference.md#qc)). An aligned pair is a base of one sequence and a base of the other
in one column, whether they match or not.

### The graph's own alignment

A graph aligns two bases when both alleles' walks read them from one node, at one offset, in one
direction. Each allele takes the walk of its first spanning run in `occurrences.tsv.gz`, turned to
the region's forward orientation; another run can spell the same allele along other nodes. A walk
that visits a node more than once, as when the graph has folded copies of a repeat into a cycle,
aligns each of its bases there with each base the other walk reads there. Such a graph can align
more pairs than any alignment can, and those pairs cost it precision. Two walks that read one node
in opposite directions align nothing there, since no column of a global alignment can hold an
inversion.

### Metrics

Each region's row of the output TSV ([format](reference.md#qctsv-and-qctsvjson)) gives its number of
alleles, `longest`, the length of its longest allele, and `sequence`, the sum over the alignment's
columns of the number of distinct bases in each, fragments included: the number of bases in the
induced graph. Its other fields are counts summed over the region's compared pairs:

| field | definition |
|---|---|
| `pair_opt`, `pair_msa`, `pair_tp` | aligned pairs in the random pairs' optimal alignments, in their projected alignments, and in both |
| `ref_opt`, `ref_msa`, `ref_tp` | the same, for the reference pairs |
| `ref_sv`, `ref_sv_opt` | in the reference pairs' projected and optimal alignments, runs of 50 or more consecutive columns in which one row has a gap opposite a base of the other, whichever row holds the gap |
| `forced`, `forced_opt` | columns of the projected and optimal alignments, of random and reference pairs alike, that lie in a forced-homology window |
| `banded` | random and reference pairs compared with a banded alignment because of `--max-cells` |
| `skipped` | random and reference pairs not compared because of `--max-memory`; every other field leaves them out |
| `graph_pair_msa`, `graph_pair_tp` | with `--graph`: aligned pairs in the graph's alignment of the random pairs, and in both it and the optimal alignment |
| `graph_ref_msa`, `graph_ref_tp` | the same, for the reference pairs |
| `graph_sequence` | with `--graph`: the bases of the nodes the region's runs visit, which `replace` removes |

A forced-homology window is any 40 consecutive columns of a two-row alignment in which at least 20
columns hold two bases and more than 30% of those are mismatches. Homologous tandem-repeat alleles
differ mostly by gaps, so such a window shows unrelated sequence aligned base against base. An
excess of `ref_sv` over `ref_sv_opt` shows long indels that the reference comparison does not need.

The summary sums every field over the regions, and adds, for the random and for the reference pairs,
precision, tp / msa, recall, tp / opt, and F1, 2 tp / (opt + msa). `growth` is the sum of `sequence`
over the sum of `longest`. It is at least 1, because each base of a region's longest allele needs a
column of its own, and its excess over 1 is the bases the graphs hold beyond their longest alleles,
as a fraction of those alleles' total length. With `--graph`, the summary has the same ratios for
the graph, prefixed `graph_`, over the same optimal alignments, and `graph_growth`, the sum of
`graph_sequence` over the sum of `longest`. A graph that folds a repeat into a cycle can have a
`graph_growth` below 1.

## Realigning on a panel: union and project

`union` and `project` realign the regions of one graph, the recipient, with the alleles of a panel:
a larger graph in which the recipient's haplotypes are walks too, as when the recipient was built
from part of the panel. `extract` is run on both graphs with one regions file, so that each region
has the same anchors in both. Then:
1. `union` merges the two extracts' packages, in the panel's order, keeping the regions present in
   both and stopping if a region's anchors differ. A recipient allele with the sequence of a panel
   allele is that allele, with the panel's weight; one that the panel lacks joins with weight 0, so
   it cannot be the centre while any allele has a positive weight. Fragments join in the same way,
   matched by side and sequence. A union package's status is the panel's, except that a region
   complex in the recipient is complex, so a region invariant in the panel is not realigned even
   when the recipient has alleles the panel lacks.
2. `realign` aligns the union's packages, so each centre is the panel's reference allele, or chosen
   from the panel's alleles.
3. `project` keeps, from each `ok` alignment of the union, the rows of the recipient's members,
   under the recipient's ids, and drops the columns that are gaps in all of them.
4. The recipient's `replace` takes the recipient's extract and project's output, and rewrites the
   regions that are `ok` in the recipient's own extract and have an alignment from `project`.

## Merge

`merge` joins per-contig GBZs into one whole-genome GBZ, without renumbering any node. Its inputs
must cover the genome: each realigned contig, and every part of the graph that was not realigned, cut
out as it is ([whole-genome.md](whole-genome.md)). It stops if two inputs' node-id ranges overlap,
comparing ranges rather than ids. It concatenates the inputs' nodes and edges into one GFA, merges
their GBWTs, the indexes of their paths, with `vg gbwt -f`, which builds the GBZ over that GFA, and
keeps as reference samples those with reference paths in any input. The output's node count, edge
count and sequence length must be the sums of the inputs', and its path names exactly the union of
theirs; an output that is not is renamed to `<output>.failed`.

## Known gaps

- `replace` rewrites walks, not generic paths: a P line that visits a replaced interior stops
  `replace`, which then writes no output.
- `prepare` accepts a P-line reference path named `<sample>#<contig>`, but gives `vg snarls -P` the
  name prefix `<sample>#0#<contig>`, which that name lacks, so the snarl tree is not oriented along
  it.
- The alignment's rows are upper case, so a replaced region loses any lower-case bases the input
  had. `replace` and `verify` compare upper-cased sequences, so neither reports this.
- The realigned graph lacks the samples `prepare` dropped, and `verify`, which compares it with
  prepare's graph, does not report it.
- When several contigs share a connected component, prepare's GBZ holds them all, but `regions`
  anchors only the targets of the contig named. Preparing another contig of the component cuts the
  same component again, and `merge` stops on the overlapping ids.
- A fragment found in no allele is aligned to the centre end to end: its side plays no part, and the
  stretch of the region it does not cover is scored as an ordinary gap.
- A node between the anchors that no recorded run visits is not in the interior, so `replace` keeps
  it, with whatever edges it has to nodes that stay.
- Two alleles that both lack k-mers are at distance 0 when the medoid is chosen exactly, and 1 when
  it is chosen from sketches.
- Around the reference, two other alleles' projected alignment is composed through the reference,
  which aligns them less well than composing through the medoid did (pair F1 0.9942 against
  0.9949, and much less at a repeat expanded far beyond the reference; see The centre).
- Under the same `--max-cells`, `qc` compares a pair of the centre and an allele that `realign`
  aligned with a band with that same banded alignment, so it agrees by construction, and the region's
  pairs over `--max-cells` are measured against an alignment that may not be optimal. `qc` has no
  limit on a region's work or time, so a region can take longer in `qc` than it did in `realign`,
  and it skips a pair over its `--max-memory`, so its measures can differ with the budget.
- `run` judges a stage done by its output alone, whatever options made it, and `merge` refuses a
  contig whose reuse numbering reached past other contigs' ids
  ([reference.md](reference.md#known-gaps)).
