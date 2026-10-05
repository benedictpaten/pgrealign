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
(`--repeatmasker`), its FASTA (`--reference`), or both, on the contigs that `--contigs` names. The
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

### Merging records into entries

The records of each contig, from both sources, are taken in order of start, then end. A record joins
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
forward orientation, a reverse run's steps reversed and each flipped. An internal run is always `+`,
spelled as its walk runs (see Known gaps).

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
| `complex` | a walk has a complex run, or an anchor-free walk makes the region complex as above |
| `no_runs` | no walk has a run in the region |
| `invariant` | there is one allele and no fragment, so there is nothing to realign |
| `ok` | any other region: two or more alleles, or at least one fragment |

A reference path that is a W line crosses every region `regions` made `ok`, so `no_runs` means that
the graph given to `extract` has a P-line reference path or is not the graph the regions were
anchored in. `extract` writes a package per region, holding its status, members and reference
allele; every run that is not complex, which `replace` rewrites; and each region's interior nodes,
which `replace` removes ([packages](reference.md#packagesjsonlgz),
[runs](reference.md#occurrencestsvgz), [interiors](reference.md#interiorstsvgz)).

## Realign: the medoid star

`realign` builds one multiple alignment for each `ok` region that it selects (see Which regions are
realigned). A multiple alignment here is a set of rows of equal length, one per member, each row
being the member's sequence with gap characters (`-`) inserted. A column is one position across all
the rows.

The method is a star alignment around a medoid, the allele with the least weighted distance to the
others, which is called the centre:
1. Centre. The centre is chosen from the alleles.
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
from each sequence, at about 15 bytes per cell. A band limits it to the cells near the alignment's
path. abPOA's own band is adaptive, reaching about 10 cells plus 1% of the sequence's length either
side of the path, and `-b -1` turns it off. The aligner is told, for each alignment, whether a band
is wanted: alignments to the centre are unbanded unless too large, and alignments of inserted
sequences to one another are banded, as the sections below say.

abPOA is given A, C, G, T and N only. Every sequence is upper-cased first, and any other character is
given to abPOA as N and restored in the returned rows. Every row is checked to spell its member's
upper-cased sequence, and all rows to have one length. The centre is chosen from the same sequences.

### Choosing the centre

The distance between two alleles comes from their 15-mers. Each allele's k-mers, with k = 15, form a
multiset: a k-mer that occurs three times counts three times. k-mers that contain a character other
than A, C, G or T are skipped. The distance between alleles i and j is one minus the Jaccard index
of their multisets:

    J(i, j) = shared / (|i| + |j| - shared)

where `shared` sums, over every k-mer, the smaller of its two copy numbers, and |i| is the number of
k-mers in i. When neither allele has a k-mer, as when both are shorter than 15 bp, J is 1, so the two
are at distance 0; an allele without k-mers is at distance 1 from one that has them. The centre is
the allele i of positive weight with the least weighted sum of distances to the others, the sum over
every other allele j of w_j (1 - J(i, j)), where w_j is j's weight. Ties go to the earlier allele.
If no allele has a positive weight, every weighted sum is 0 and the centre is the first allele.

The multiset, rather than the set, carries copy number: tandem-repeat alleles differ mostly in how
many copies of the motif they hold, and two that differ only in copy number often contain the same
distinct k-mers. Because the distances are weighted, the centre lies near the sequences that most
runs spell. An allele can have weight 0 when `union` adds it (see Realigning on a panel).

### Estimating distances from sketches

Comparing one pair exactly costs about as many steps as the smaller allele has distinct k-mers. When
the sum over all pairs exceeds `star.MEDOID_WORK`, 10^9 steps, every distance is estimated from
sketches instead, which makes each pair cost about the same however long the alleles are. An
allele's sketch is the `star.SKETCH`, 1,000, smallest hash values of its k-mer multiset, each copy
of a k-mer hashed as an element of its own, so that copy number still counts. J(i, j) is estimated
as the fraction of the 1,000 smallest values of the two sketches' union that lie in both sketches;
two empty sketches are at distance 1 (see Known gaps). The hash is a fixed function, so the
sketches, and the centre, are the same in every run.

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
| `too_big` | even a banded alignment of the region's longest allele or fragment would exceed `--max-cells` cells, or its alleles hold more than 40,000,000 bases in all |
| `timeout` | the region ran longer than `--timeout` seconds |
| `aligner_failed` | abPOA exited with an error, left out a row, or returned pairwise rows that do not spell their inputs |
| `no_alleles` | the region has fragments but no allele: no walk spans it |
| `error` | any other failure, recorded with its message: a bug, a slot alignment whose rows the final check rejects, or abPOA could not be started |

Both size limits are checked before any work. A banded alignment of a sequence of length n to one
of the same length is taken to fill the smaller of n² and 2n (10 + ⌊n / 100⌋) cells, abPOA's
adaptive band, n being the length of the region's longest allele or fragment. For long sequences
this is about n² / 50, so the limit admits sequences about seven times as long as an unbanded
alignment of two equal lengths would. The slot alignments are not counted. The limit on total bases
exists because choosing the centre holds every allele's k-mer counts at once. The time limit is
checked while choosing the centre and before each call to abPOA, and abPOA itself is stopped at the
deadline. Each region's failure is recorded in its own result, so one region cannot stop a contig.

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
`realign` banded is compared with the very alignment `realign` used. An aligned pair is a base of
one sequence and a base of the other in one column, whether they match or not.

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
2. `realign` aligns the union's packages, so the panel's alleles choose each centre.
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
- An internal fragment from a walk that runs against the region's forward orientation is the reverse
  complement of what the alleles hold, so it is usually found in no allele and aligns poorly to the
  centre, though its walk still spells what it did.
- A fragment found in no allele is aligned to the centre end to end: its side plays no part, and the
  stretch of the region it does not cover is scored as an ordinary gap.
- A node between the anchors that no recorded run visits is not in the interior, so `replace` keeps
  it, with whatever edges it has to nodes that stay.
- Two alleles that both lack k-mers are at distance 0 when the centre is chosen exactly, and 1 when
  it is chosen from sketches.
- Under the same `--max-cells`, `qc` compares a pair of the centre and an allele that `realign`
  aligned with a band with that same banded alignment, so it agrees by construction, and the region's
  pairs over `--max-cells` are measured against an alignment that may not be optimal. `qc`'s
  alignments have no time limit, so a region can take longer in `qc` than it did in `realign`.
- `run` judges a stage done by its output alone, whatever options made it, and `merge` refuses a
  contig whose reuse numbering reached past other contigs' ids
  ([reference.md](reference.md#known-gaps)).
