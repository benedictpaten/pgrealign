# pgrealign

pgrealign realigns chosen regions of a pangenome graph, such as its tandem repeats.

A pangenome graph holds many haplotypes, each one or more walks through one shared graph of
sequence nodes. A target is an interval of the reference genome to realign, and its region is the
part of the graph between two anchor nodes that every walk crossing the target must pass through.
In a region, a walk spells an allele when it crosses the region, and a fragment when it starts or
ends inside it.

For each region, pgrealign:

1. extracts the alleles and fragments that the walks spell there;
2. builds one multiple alignment of them, a star: the alleles are aligned to a centre allele, the
   reference's, so that each allele's alignment to the reference is its optimal one, or, in a region
   the reference does not span or where most haplotypes are much longer than the reference, the one
   closest in k-mer content to the others, each counted once for every walk that spells it; the pairwise alignments are merged into one, a fragment that occurs
   in an allele takes that allele's row, and any other fragment is aligned to the centre with the
   alleles;
3. replaces the nodes between the anchors with a graph built from the alignment, and rewrites each
   walk through it.

Every walk spells the same sequence afterwards as before, up to letter case
([Known gaps](#known-gaps)). Only the nodes change: how the walks share them inside realigned
regions and, under replace's default numbering, every node's id.

## Install

pgrealign needs Python 3.9 or later and uses the standard library only. From the repository root:

```
pip install -e .
```

This installs the `pgrealign` command. It runs these programs, found on `PATH` unless an option
names the binary:

- vg, the pangenome toolkit (`--vg`), in prepare, regions, verify, merge and run, in extract when it
  reads a GBZ, vg's compressed format for a graph and its paths, or a graph in another of vg's
  formats, and in replace when it reads or writes one;
- abPOA, a partial-order sequence aligner (`--abpoa`), version 1.5 or later (1.5.7 is tested), in
  realign, qc and run;
- a C compiler, for `catalog --reference`, which compiles its repeat finders on first use, trfind
  and, with `--strf`, strf, unless `--trfind` and `--strf` name ones built elsewhere
  ([docs/reference.md](docs/reference.md#catalog) says which compiler and where the build is kept);
- TRF, Tandem Repeats Finder, for `catalog --trf`, which names the binary. `catalog --strf` runs
  strf, pgrealign's own reimplementation of the part of TRF that catalog uses, in its place.

## Quick start

`pgrealign run` takes one contig of a whole-genome GBZ through every stage from prepare to qc:

```
pgrealign run --gbz graph.gbz --contig chr20 \
    --targets catalogs/chm13v2.0_tr.bed.gz -o out -j 8 -t 8
```

The reference sample (`--ref-sample`) gives the graph its coordinates, and the GBZ must hold its
walk on the contig, `<sample>#0#chr20`. `--targets` is a BED of intervals on that walk:
`catalogs/chm13v2.0_tr.bed.gz` lists the tandem repeats of CHM13 v2.0 on chr1-22, chrX and chrY,
and `pgrealign catalog` builds such a catalogue for another reference. `-j` sets the most regions
realign and qc work on at once, and can be the number of cores: they start a region only while the
memory estimates of those running fit in `--max-memory`, half the machine's memory by default.
Which regions realign aligns does not depend on that budget or on the machine: it decides from the
size of each region's alignments, in cells, before aligning anything, and has no time limit. `-t`
sets how many threads vg uses. The outputs go to `out/chr20/`:

| output | stage | contents |
|---|---|---|
| `chr20.gbz`, `chr20.snarls.pb`, `chr20.ref.tsv.gz` | prepare | the contig's graph, its snarls (the parts of the graph that two nodes separate from the rest, as vg finds them), and the reference walk's node positions |
| `regions.tsv.gz` | regions | each region's anchors and status |
| `extract/` | extract | the alleles and fragments of each region that regions made `ok`, and each walk's passage through it |
| `msas.jsonl.gz` | realign | the alignment of each region that extract passed as `ok`, or why realign could not build one |
| `chr20.realigned.gbz`, `replace.json` | replace | the realigned contig, and a summary that gives its node-id range |
| `verify.json` | verify | the comparison of every path with its counterpart in `chr20.gbz` |
| `qc.tsv`, `qc.tsv.json` | qc | per region and in total, how well the alignments agree with optimal pairwise alignments, before and after realignment |

### Resuming and options

Each stage also writes a manifest beside its output, recording its command, inputs, parameters and
tools. run runs each stage whose output is missing and skips each stage whose output exists, so the
same command again carries on where it stopped; an interrupted realign resumes from
`msas.jsonl.gz.partial`. run stops with an error when verify fails, and again on every later run
until `verify.json` and the realigned GBZ are deleted.

run passes on most stage options. To change one it lacks, such as regions' `--pad`, run that stage
and those before it as subcommands and then `pgrealign run` for the rest
([docs/whole-genome.md](docs/whole-genome.md#options-that-run-does-not-pass);
[docs/reference.md](docs/reference.md#run) lists the options).

## A whole genome

Each contig is a separate `pgrealign run`, and contigs are independent until they are merged, so
they can run at once. `pgrealign merge` then joins the realigned contigs and the parts of the graph
that were not realigned, each a GBZ of its own, into one GBZ. merge renumbers no node, and stops when
two inputs' node-id ranges overlap. [docs/whole-genome.md](docs/whole-genome.md) gives the
commands, how to resolve an overlap, and what each contig costs in memory and disk.

## The stages

Each stage is a subcommand that reads and writes files, so any one can run on its own;
`pgrealign <stage> --help` lists its options. regions gives every region a status, extract gives one
to each region that regions made `ok`, and realign gives one to each `ok` package that `--select`
picks. Only a region that is `ok` at all three is replaced; every other region keeps its subgraph.

| stage | what it does |
|---|---|
| `catalog` | builds a reference's tandem-repeat catalogue from its RepeatMasker annotation, its FASTA, or both; `--trf` runs TRF on the FASTA as RepeatMasker would, for a reference with no annotation, and `--strf` runs pgrealign's reimplementation of it instead |
| `prepare` | cuts one contig's component, the connected part of the graph that holds its reference walk, out of a whole-genome GBZ, with its snarls and the reference walk's node positions, and by default drops the gref_CHM13 cover ([Known gaps](#known-gaps)); with `--graph`, writes the snarls and positions of a graph that is already the contig's, such as a Minigraph-Cactus chromosome `.vg` |
| `regions` | places each target between two anchor nodes, and merges regions that would interfere |
| `extract` | reads the contig's walks once, and writes the alleles and fragments of each `ok` region as a package |
| `realign` | builds the star alignment of each `ok` package, several at once with `-j`; `--plan` says, without aligning, which it would align and how much memory they need |
| `replace` | puts each alignment's graph in place of its region's nodes, and rewrites the walks through it |
| `verify` | checks that every path of the realigned contig spells what it did before, up to letter case |
| `qc` | measures how well each alignment, and the original graph's, agrees with optimal pairwise alignments |
| `run` | runs prepare to qc for one contig, resuming where an earlier run stopped |
| `merge` | joins per-contig GBZs whose node ids do not overlap into one whole-genome GBZ |

## Realigning on a panel

union and project realign regions on the haplotypes of one graph, the panel, and apply the result
to another graph, the recipient, whose walks use the panel graph's nodes, such as the panel graph
with some samples removed. The recipient's regions are then aligned as in the panel, whose alleles
alone choose each centre. extract runs at the same anchors on both graphs, and then:

| stage | what it does |
|---|---|
| `union` | joins the panel's packages with the recipient's, region by region, for realign |
| `project` | restricts the union's alignments to the recipient's alleles and fragments, for the recipient's replace |

[docs/reference.md](docs/reference.md#union-and-project) gives the steps.

## Calling the alignment from Python

`pgrealign/star.py` holds the alignment method on its own. `star_msa` takes sequences and returns
rows; it reads no graph and no file:

```python
from pgrealign.star import Abpoa, star_msa

alleles = [('a1', 'GATTACAGATTACAGATTACA', 5),          # (id, sequence, weight)
           ('a2', 'GATTACAGATTACA', 2),
           ('a3', 'GATTACAGATTCCAGATTACAGATTACA', 1)]
fragments = [('f1', 'TTACAGATTACA', 'right')]           # (id, sequence, side)

with Abpoa() as align:
    msa = star_msa(alleles, fragments, align=align)
print(msa.centre)                                       # the centre allele's id
for name, row in msa.rows.items():                      # one gapped row per allele and fragment
    print(name, row)
```

- A weight is any number. Only an allele of positive weight can be the centre, unless none has one;
  in the pipeline, the weight is the number of walk crossings that spell the allele.
- A fragment's side is `left` when it touches the region's left end, `right` when it touches the
  right end, and `internal` when it touches neither.
- Every row has the same length and, once its `-` gaps are removed, spells its input upper-cased.
- `align` is an `Abpoa` inside its `with` block, which runs abpoa with the method's scores, or any
  callable `align(seqs, banded=False, seeded=False)` that returns a global alignment of `seqs`, one
  row per sequence in order. Pass the aligner's scores as `scores`.
- `star_msa` asks for `banded=True` when it aligns to each other the distinct insertions that fall
  between the same two centre bases. With `max_cells`, it aligns a pair whose dynamic-programming
  matrix, the product of the two lengths, would exceed it twice, with `banded` the band's margin in
  cells and with `banded=True, seeded=True`, and keeps the one that scores higher under `scores`.

[docs/reference.md](docs/reference.md#the-python-api) describes the arguments and the aligner contract.

## Known gaps

- Nothing schedules the contigs of a genome: there is no workflow-manager packaging.
- In GFA, the text format for a graph, a walk is a W line and a generic path a P line. replace
  rewrites walks only: a P line that visits a node replace would remove makes replace stop with an
  error.
- A replaced region's bases come out upper-case, so any soft-masking there (lower-case bases marking
  repeats) is lost.
- By default the realigned graph has no gref_CHM13 cover: the copy of the reference path, plus the
  pieces of other haplotypes that hang off it, that `vg paths -u` computes. prepare drops it because
  those pieces start and end inside regions.
  [docs/whole-genome.md](docs/whole-genome.md#the-gref_chm13-cover) gives the vg commands that
  rebuild it once verify has passed.
- run takes a stage as done when its output exists, whatever options made it: after changing an
  option, delete the outputs of that stage and of every stage after it.
- abPOA reserves address space for the whole matrix of every alignment, even a banded one: 256 GB
  for two 100 kb alleles, of which a band writes a few GB. The abPOA that cactus builds maps it so
  that Linux does not count it; with another abPOA, Linux refuses a reservation larger than its
  memory and swap, and realign stops on such a region, to be resumed elsewhere
  ([docs/reference.md](docs/reference.md#memory)).

## Documentation

- [docs/method.md](docs/method.md): the method, stage by stage, with the reasoning behind each step.
- [docs/reference.md](docs/reference.md): every subcommand, the files it reads and writes, and the
  Python API.
- [docs/whole-genome.md](docs/whole-genome.md): building a realigned whole-genome graph from the
  contigs, and what each contig costs in memory and disk.
- The docstring at the top of each stage's module, which states that stage's rules in full.

## Tests

```
python3 -m unittest discover -s tests -v       # from the repository root
```

A test that calls abpoa, vg or a C compiler is skipped, and reported as skipped, when that program
is not on `PATH`; the compiler tests look for `cc`, `clang` or `gcc` and ignore `$CC`. The
end-to-end test takes a two-contig genome through run and merge, and needs both vg and abpoa.

## Layout

```
pyproject.toml    the package's build settings and the pgrealign command
pgrealign/        the package: one module per stage, except run (pipeline.py) and project (union.py)
  star.py         the star alignment, which reads no graph and no file
  induce.py       the graph of an alignment, used by replace
  graph.py        reading a graph as a stream of GFA lines
  memory.py       abPOA's memory, and running regions in parallel within a memory budget
  provenance.py   the manifest written beside each stage's output
  cli.py          the pgrealign command line
  c/trfind.c      the tandem-repeat finder that catalog compiles
  c/strf.c        the reimplementation of TRF's short repeats that catalog compiles for --strf
catalogs/         the CHM13 v2.0 tandem-repeat catalogue and its manifest
docs/             method.md, reference.md and whole-genome.md
tests/            unit and end-to-end tests
```

## Licence

MIT; see [LICENSE](LICENSE).
