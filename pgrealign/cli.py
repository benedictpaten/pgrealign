"""The `pgrealign` command line: one subcommand per stage, plus `run` for every stage of one contig and
`merge` for a whole genome. docs/reference.md describes each one, and the files they read and write."""
import argparse
import json
import os
import sys
import time

from . import __version__


def _tools(a, *names):
    """{tool: the binary this command ran}, for the manifest."""
    return {n: getattr(a, n, None) for n in names}


def _print(obj):
    json.dump(obj, sys.stdout, indent=1)
    print()


def cmd_catalog(a):
    from . import catalog, provenance
    started = time.time()
    contigs = None if a.contigs == 'all' else (a.contigs.split(',') if a.contigs else catalog.DEFAULT_CONTIGS)
    summary = catalog.build(a.out, repeatmasker=a.repeatmasker, fastas=a.reference or (), gap=a.gap,
                            contigs=contigs, minp=a.min_period, maxp=a.max_period, jobs=a.jobs)
    inputs = {'repeatmasker': a.repeatmasker} if a.repeatmasker else {}
    inputs.update({'reference_%d' % i: f for i, f in enumerate(a.reference or ())})
    provenance.write_manifest(a.out, 'catalog', inputs,
                              {'gap': a.gap, 'contigs': contigs, 'min_period': a.min_period,
                               'max_period': a.max_period}, extra={'summary': summary}, started=started,
                              checksums=True)
    _print(summary)


def drop_samples_of(a):
    """The samples prepare removes: none with --keep-all-samples, else --drop-sample, else the default."""
    from . import prepare
    if a.keep_all_samples:
        return ()
    return tuple(a.drop_sample) if a.drop_sample else prepare.DROP_SAMPLES


def cmd_prepare(a):
    from . import prepare, provenance
    started = time.time()
    drop = drop_samples_of(a)
    out = prepare.prepare(a.gbz, a.contig, a.out, ref_sample=a.ref_sample, drop_samples=drop, vg=a.vg,
                          threads=a.threads)
    provenance.write_manifest(os.path.join(a.out, a.contig), 'prepare', {'gbz': a.gbz},
                              {'contig': a.contig, 'ref_sample': a.ref_sample, 'drop_samples': list(drop)},
                              tools=_tools(a, 'vg'), extra={'outputs': out}, started=started)
    _print(out)


def cmd_regions(a):
    from . import provenance, regions
    started = time.time()
    summary = regions.build(a.targets, a.contig, a.snarls, a.ref, a.out, vg=a.vg, pad=a.pad, max_span=a.max_span)
    provenance.write_manifest(a.out, 'regions', {'targets': a.targets, 'snarls': a.snarls, 'ref': a.ref},
                              {'contig': a.contig, 'pad': a.pad, 'max_span': a.max_span}, tools=_tools(a, 'vg'),
                              extra={'summary': summary}, started=started)
    _print(summary)


def cmd_extract(a):
    from . import extract, provenance, regions
    started = time.time()
    rs = list(regions.read_regions(a.regions))
    summary = extract.extract(a.graph, rs, a.out, vg=a.vg, threads=a.threads, ref_sample=a.ref_sample,
                              max_internal_bp=a.max_internal_bp)
    provenance.write_manifest(os.path.join(a.out, 'extract'), 'extract', {'graph': a.graph, 'regions': a.regions},
                              {'ref_sample': a.ref_sample, 'max_internal_bp': a.max_internal_bp},
                              tools=_tools(a, 'vg'), extra={'summary': summary}, started=started)
    _print(summary)


def scores_of(a):
    from . import star
    return star.Scores.parse(a.scores) if a.scores is not None else star.SCORES


def cmd_realign(a):
    from . import provenance, realign
    started = time.time()
    i, n = (int(x) for x in a.shard.split('/'))
    scores = scores_of(a)
    summary = realign.run(a.packages, a.out, select=a.select, shard=(i, n), jobs=a.jobs, abpoa=a.abpoa,
                          timeout=a.timeout, max_cells=a.max_cells, scores=scores)
    provenance.write_manifest(a.out, 'realign', {'packages': a.packages},
                              {'method': 'medoid_star', 'select': a.select, 'shard': a.shard, 'timeout': a.timeout,
                               'max_cells': a.max_cells, 'scores': ' '.join(scores.flags())},
                              tools=_tools(a, 'abpoa'), extra={'summary': summary}, started=started)
    _print(summary)


def msas_inputs(msas):
    return {'msas': msas[0]} if len(msas) == 1 else {'msas_%d' % i: m for i, m in enumerate(msas)}


def cmd_replace(a):
    from . import provenance, replace
    started = time.time()
    summary = replace.replace(a.graph, a.extract, a.msas, a.out, a.id_start, vg=a.vg, threads=a.threads,
                              id_mode=a.id_mode)
    if a.gbz:
        replace.to_gbz(a.out, a.gbz, vg=a.vg, keep_gfa=a.keep_gfa)
        summary['gbz'] = a.gbz
    provenance.write_manifest(a.gbz or a.out, 'replace', {'graph': a.graph, **msas_inputs(a.msas)},
                              {'extract': a.extract, 'id_mode': a.id_mode, 'id_start': a.id_start},
                              tools=_tools(a, 'vg'), extra={'summary': summary}, started=started)
    _print(summary)


def cmd_verify(a):
    from . import provenance, verify
    started = time.time()
    res = verify.verify(a.before, a.after, vg=a.vg)
    if a.out:
        with open(a.out, 'w') as f:
            json.dump(res, f, indent=1)
        provenance.write_manifest(a.out, 'verify', {'before': a.before, 'after': a.after}, {},
                                  tools=_tools(a, 'vg'), started=started)
    _print(res)
    if not res['ok']:
        sys.exit(1)


def cmd_qc(a):
    from . import provenance, qc
    started = time.time()
    scores = scores_of(a)
    summary = qc.run(a.extract, a.msas, a.out, pairs=a.pairs, refs=a.refs, sample=a.sample, jobs=a.jobs,
                     abpoa=a.abpoa, scores=scores, ref_sample=a.ref_sample, max_cells=a.max_cells)
    provenance.write_manifest(a.out, 'qc', {'msas': a.msas},
                              {'extract': a.extract, 'pairs': a.pairs, 'refs': a.refs, 'sample': a.sample,
                               'scores': ' '.join(scores.flags())},
                              tools=_tools(a, 'abpoa'), extra={'summary': summary}, started=started)
    _print(summary)


def cmd_run(a):
    from . import pipeline
    summary = pipeline.run_contig(a.gbz, a.contig, a.targets, a.out, jobs=a.jobs, threads=a.threads, vg=a.vg,
                                  abpoa=a.abpoa, ref_sample=a.ref_sample, drop_samples=drop_samples_of(a),
                                  scores=scores_of(a), timeout=a.timeout, max_cells=a.max_cells,
                                  max_internal_bp=a.max_internal_bp, id_start=a.id_start, do_qc=not a.no_qc,
                                  log=lambda m: print(m, file=sys.stderr))
    _print(summary)


def cmd_merge(a):
    from . import merge, provenance
    started = time.time()
    summary = merge.merge(a.gbz, a.out, vg=a.vg)
    provenance.write_manifest(a.out, 'merge', {'gbz_%d' % i: g for i, g in enumerate(a.gbz)}, {},
                              tools=_tools(a, 'vg'), extra={'summary': summary}, started=started)
    _print(summary)


def cmd_union(a):
    from . import provenance, union
    started = time.time()
    summary = union.union_packages(a.base, a.add, a.out)
    provenance.write_manifest(a.out, 'union', {'base': a.base, 'add': a.add}, {}, extra={'summary': summary},
                              started=started)
    _print(summary)


def cmd_project(a):
    from . import provenance, union
    started = time.time()
    summary = union.project(a.union, a.msas, a.out)
    provenance.write_manifest(a.out, 'project', {'union': a.union, 'msas': a.msas}, {},
                              extra={'summary': summary}, started=started)
    _print(summary)


def main(argv=None):
    from .prepare import DROP_SAMPLES
    from .realign import MAX_CELLS
    from .star import SCORES
    default_drop = ','.join(DROP_SAMPLES) or 'none'
    scores_help = ("abPOA's scores for every alignment, as its flags -M, -X, -O and -E; a flag left out takes "
                   "abPOA's default, so '' gives abPOA's defaults [%s]" % ' '.join(SCORES.flags()))
    p = argparse.ArgumentParser(prog='pgrealign', description=__doc__)
    p.add_argument('--version', action='version', version='pgrealign ' + __version__)
    sub = p.add_subparsers(dest='stage', required=True)

    def vg_arg(c):
        c.add_argument('--vg', help='vg binary [vg on PATH]')

    def drop_args(c):
        g = c.add_mutually_exclusive_group()
        g.add_argument('--drop-sample', action='append',
                       help='a sample to remove from the contig before realigning (repeatable) [%s]' % default_drop)
        g.add_argument('--keep-all-samples', action='store_true', help='remove no sample')

    c = sub.add_parser('catalog', help='the tandem-repeat catalogue of a reference')
    c.add_argument('-o', '--out', required=True, help='output catalogue (.bed or .bed.gz)')
    c.add_argument('--repeatmasker', help='RepeatMasker annotation as BED (10 columns)')
    c.add_argument('--reference', nargs='+', metavar='FASTA',
                   help='reference FASTA file(s) to run trfind on (one per contig runs them in parallel)')
    c.add_argument('--gap', type=int, default=50, help='merge records that overlap or lie within this many bp [50]')
    c.add_argument('--min-period', type=int, default=7, help='trfind smallest period [7]')
    c.add_argument('--max-period', type=int, default=2000, help='trfind largest period [2000]')
    c.add_argument('--contigs', default=None, help="comma-separated contigs, or 'all' [chr1-22,chrX]")
    c.add_argument('-j', '--jobs', type=int, default=1, help='trfind processes [1]')
    c.set_defaults(func=cmd_catalog)

    c = sub.add_parser('prepare', help='cut one contig from a whole-genome GBZ; its snarls and reference index')
    c.add_argument('--gbz', required=True, help='whole-genome GBZ')
    c.add_argument('--contig', required=True)
    c.add_argument('-o', '--out', required=True, help='output directory')
    c.add_argument('--ref-sample', default='CHM13', help='the reference sample [CHM13]')
    drop_args(c)
    vg_arg(c)
    c.add_argument('-t', '--threads', type=int, default=4)
    c.set_defaults(func=cmd_prepare)

    c = sub.add_parser('regions', help="anchor every target in the contig's snarls, and merge those that touch")
    c.add_argument('--targets', required=True, help='catalogue or BED of target intervals')
    c.add_argument('--contig', required=True)
    c.add_argument('--snarls', required=True, help='<contig>.snarls.pb from prepare')
    c.add_argument('--ref', required=True, help='<contig>.ref.tsv.gz from prepare')
    c.add_argument('-o', '--out', required=True, help='regions.tsv(.gz)')
    c.add_argument('--pad', type=int, default=200, help='pad each target by this many bp [200]')
    c.add_argument('--max-span', type=int, default=100000,
                   help='leave a region whose anchors are further apart than this as it is [100000]')
    vg_arg(c)
    c.set_defaults(func=cmd_regions)

    c = sub.add_parser('extract', help="every walk's run through every region, in one streaming pass")
    c.add_argument('--graph', required=True, help='<contig>.gbz from prepare, or a GFA')
    c.add_argument('--regions', required=True, help='regions.tsv(.gz) from regions')
    c.add_argument('-o', '--out', required=True, help='output directory')
    c.add_argument('--ref-sample', default='CHM13', help='the reference sample, whose allele each package names [CHM13]')
    c.add_argument('--max-internal-bp', type=int, default=200000,
                   help="an anchor-free walk longer than this inside a region's interior makes the region complex, "
                        "left as it is [200000]")
    vg_arg(c)
    c.add_argument('-t', '--threads', type=int, default=4)
    c.set_defaults(func=cmd_extract)

    c = sub.add_parser('realign', help='the medoid-star alignment of every selected region')
    c.add_argument('--packages', required=True, help='packages.jsonl.gz from extract')
    c.add_argument('-o', '--out', required=True, help='msas.jsonl.gz')
    c.add_argument('--select', choices=['any', 'length'], default='any',
                   help="realign a region when its alleles differ at all ('any') or in length [any]")
    c.add_argument('--shard', default='0/1', help='do shard i of n, by package order [0/1]')
    c.add_argument('-j', '--jobs', type=int, default=1, help='regions aligned in parallel [1]')
    c.add_argument('--abpoa', default='abpoa', help='abPOA binary [abpoa on PATH]')
    c.add_argument('--timeout', type=float, default=900, help='seconds per region [900]')
    c.add_argument('--max-cells', type=int, default=MAX_CELLS,
                   help='band a pairwise alignment that would exceed this many cells unbanded, and leave a region as '
                        'it is when even a banded one would [%d]' % MAX_CELLS)
    c.add_argument('--scores', metavar='FLAGS', help=scores_help)
    c.set_defaults(func=cmd_realign)

    c = sub.add_parser('replace', help="splice the realigned regions into the contig's graph")
    c.add_argument('--graph', required=True, help='<contig>.gbz from prepare, or a GFA')
    c.add_argument('--extract', required=True, help='the extract output directory')
    c.add_argument('--msas', required=True, nargs='+', help='msas.jsonl.gz from realign, or one per shard')
    c.add_argument('-o', '--out', required=True, help='output GFA')
    c.add_argument('--id-mode', choices=['dense', 'reuse'], default='dense',
                   help="'dense' renumbers every node, new ones next to their regions; 'reuse' keeps old ids "
                        "and puts overflow at --id-start [dense]")
    c.add_argument('--id-start', type=int, help="dense: the first id [the contig's smallest id]; reuse: the first "
                   "id for new nodes that the region's freed ids do not cover (required)")
    c.add_argument('--gbz', help='also build this GBZ from the output GFA (vg gbwt -G)')
    c.add_argument('--keep-gfa', action='store_true', help='keep the GFA after --gbz')
    vg_arg(c)
    c.add_argument('-t', '--threads', type=int, default=4)
    c.set_defaults(func=cmd_replace)

    c = sub.add_parser('verify', help='every path of the realigned graph spells what it did before')
    c.add_argument('--before', required=True, help="the graph given to replace --graph: prepare's <contig>.gbz")
    c.add_argument('--after', required=True, help='the realigned contig GBZ')
    c.add_argument('-o', '--out', help='write the report here too')
    vg_arg(c)
    c.set_defaults(func=cmd_verify)

    c = sub.add_parser('qc', help="how well each region's alignment agrees with the pairwise alignments it implies")
    c.add_argument('--extract', required=True, help='the extract output directory')
    c.add_argument('--msas', required=True, help='msas.jsonl.gz from realign')
    c.add_argument('-o', '--out', required=True, help='per-region TSV; the summary goes to OUT.json')
    c.add_argument('--pairs', type=int, default=20, help='random allele pairs per region [20]')
    c.add_argument('--refs', type=int, default=20, help='heaviest alleles compared with the reference allele [20]')
    c.add_argument('--sample', type=int, help='measure this many random regions rather than all')
    c.add_argument('--ref-sample', default='CHM13', help='the reference sample, for packages that do not name it [CHM13]')
    c.add_argument('--scores', metavar='FLAGS', help="abPOA's scores for the optimal pairwise alignments; use the "
                   "scores the alignment was built with [%s]" % ' '.join(SCORES.flags()))
    c.add_argument('--max-cells', type=int, default=MAX_CELLS,
                   help='skip a pair whose unbanded optimal alignment would exceed this many cells [%d]' % MAX_CELLS)
    c.add_argument('-j', '--jobs', type=int, default=1, help='regions measured in parallel [1]')
    c.add_argument('--abpoa', default='abpoa', help='abPOA binary [abpoa on PATH]')
    c.set_defaults(func=cmd_qc)

    c = sub.add_parser('run', help='every stage for one contig, resuming where an earlier run stopped')
    c.add_argument('--gbz', required=True, help='whole-genome GBZ')
    c.add_argument('--contig', required=True)
    c.add_argument('--targets', required=True, help='catalogue or BED of target intervals')
    c.add_argument('-o', '--out', required=True, help='output directory; the contig goes in OUT/CONTIG')
    c.add_argument('-j', '--jobs', type=int, default=1, help='regions aligned in parallel, in realign and qc [1]')
    c.add_argument('-t', '--threads', type=int, default=4, help='threads for vg [4]')
    c.add_argument('--ref-sample', default='CHM13', help='the reference sample [CHM13]')
    drop_args(c)
    c.add_argument('--scores', metavar='FLAGS', help=scores_help)
    c.add_argument('--timeout', type=float, default=900, help='seconds per region [900]')
    c.add_argument('--max-cells', type=int, default=MAX_CELLS, help='as for realign [%d]' % MAX_CELLS)
    c.add_argument('--max-internal-bp', type=int, default=200000, help='as for extract [200000]')
    c.add_argument('--id-start', type=int, help="replace's first node id [the contig's smallest id]")
    c.add_argument('--no-qc', action='store_true', help='skip the qc stage')
    c.add_argument('--abpoa', default='abpoa', help='abPOA binary [abpoa on PATH]')
    vg_arg(c)
    c.set_defaults(func=cmd_run)

    c = sub.add_parser('merge', help='one whole-genome GBZ from per-contig GBZs whose node ids do not overlap')
    c.add_argument('gbz', nargs='+', help='the per-contig GBZs: realigned contigs and the parts left as they are')
    c.add_argument('-o', '--out', required=True, help='the merged GBZ')
    vg_arg(c)
    c.set_defaults(func=cmd_merge)

    c = sub.add_parser('union', help="merge a panel's packages with a target graph's, region by region")
    c.add_argument('--base', required=True, help="the panel graph's packages.jsonl.gz")
    c.add_argument('--add', required=True, help="the target graph's packages.jsonl.gz (same regions)")
    c.add_argument('-o', '--out', required=True, help='union packages.jsonl.gz (realign this)')
    c.set_defaults(func=cmd_union)

    c = sub.add_parser('project', help="the union's alignments restricted to the target's members")
    c.add_argument('--union', required=True, help='union packages.jsonl.gz')
    c.add_argument('--msas', required=True, help="realign's msas.jsonl.gz for the union")
    c.add_argument('-o', '--out', required=True, help="msas.jsonl.gz for the target's replace")
    c.set_defaults(func=cmd_project)

    a = p.parse_args(argv)
    a.func(a)


if __name__ == '__main__':
    main()
