"""pgrealign: realign regions of a pangenome graph.

The stages are subcommands of the `pgrealign` command line and talk to each other through files, so any
one of them can be run, replaced or parallelised on its own (docs/method.md, docs/reference.md). The
alignment method itself is star.star_msa, which takes sequences and returns rows and can be called on its
own.
"""

__version__ = "0.2.0"
