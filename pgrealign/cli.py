"""The `pgrealign` command line: one subcommand per stage (docs/pipeline-design.md)."""
import argparse
import json
import sys

from . import __version__


def cmd_catalog(a):
    from . import catalog, provenance
    contigs = None if a.contigs == 'all' else (a.contigs.split(',') if a.contigs else catalog.DEFAULT_CONTIGS)
    summary = catalog.build(a.out, repeatmasker=a.repeatmasker, fastas=a.reference or (), gap=a.gap,
                            contigs=contigs, minp=a.min_period, maxp=a.max_period, jobs=a.jobs)
    inputs = {'repeatmasker': a.repeatmasker} if a.repeatmasker else {}
    inputs.update({'reference_%d' % i: f for i, f in enumerate(a.reference or ())})
    provenance.write_manifest(a.out, 'catalog', inputs,
                              {'gap': a.gap, 'contigs': contigs, 'min_period': a.min_period,
                               'max_period': a.max_period}, extra={'summary': summary})
    json.dump(summary, sys.stdout, indent=1)
    print()


def main(argv=None):
    p = argparse.ArgumentParser(prog='pgrealign', description=__doc__)
    p.add_argument('--version', action='version', version='pgrealign ' + __version__)
    sub = p.add_subparsers(dest='stage', required=True)

    c = sub.add_parser('catalog', help='tandem-repeat catalogue of the reference (stage 0)')
    c.add_argument('-o', '--out', required=True, help='output catalogue (.bed or .bed.gz)')
    c.add_argument('--repeatmasker', help='RepeatMasker annotation as BED (10 columns)')
    c.add_argument('--reference', nargs='+', metavar='FASTA',
                   help='reference FASTA file(s) to run trfind on (one per contig runs them in parallel)')
    c.add_argument('--gap', type=int, default=50,
                   help='merge records that overlap or lie within this many bp [50]')
    c.add_argument('--min-period', type=int, default=7, help='trfind smallest period [7]')
    c.add_argument('--max-period', type=int, default=2000, help='trfind largest period [2000]')
    c.add_argument('--contigs', default=None, help="comma-separated contigs, or 'all' [chr1-22,chrX]")
    c.add_argument('-j', '--jobs', type=int, default=1, help='trfind processes [1]')
    c.set_defaults(func=cmd_catalog)

    a = p.parse_args(argv)
    a.func(a)


if __name__ == '__main__':
    main()
