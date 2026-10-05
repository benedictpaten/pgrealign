"""The `pgrealign` command line: one subcommand per stage (docs/pipeline-design.md)."""
import argparse
import json
import os
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


def cmd_prepare(a):
    from . import prepare, provenance
    out = prepare.prepare(a.gbz, a.contig, a.out, ref_sample=a.ref_sample, drop_samples=a.drop_sample,
                          vg=a.vg, threads=a.threads)
    provenance.write_manifest(os.path.join(a.out, a.contig), 'prepare', {'gbz': a.gbz},
                              {'contig': a.contig, 'ref_sample': a.ref_sample, 'drop_sample': a.drop_sample},
                              tools=['vg'], extra={'outputs': out})
    json.dump(out, sys.stdout, indent=1)
    print()


def cmd_regions(a):
    from . import provenance, regions
    summary = regions.build(a.targets, a.contig, a.snarls, a.ref, a.out, vg=a.vg, pad=a.pad,
                            max_span=a.max_span)
    provenance.write_manifest(a.out, 'regions', {'targets': a.targets, 'snarls': a.snarls, 'ref': a.ref},
                              {'contig': a.contig, 'pad': a.pad, 'max_span': a.max_span}, tools=['vg'],
                              extra={'summary': summary})
    json.dump(summary, sys.stdout, indent=1)
    print()


def cmd_extract(a):
    from . import extract, provenance, regions
    rs = list(regions.read_regions(a.regions))
    summary = extract.extract(a.graph, rs, a.out, vg=a.vg, threads=a.threads, ref_sample=a.ref_sample)
    provenance.write_manifest(os.path.join(a.out, 'extract'), 'extract',
                              {'graph': a.graph, 'regions': a.regions}, {'ref_sample': a.ref_sample}, tools=['vg'],
                              extra={'summary': summary})
    json.dump(summary, sys.stdout, indent=1)
    print()


def cmd_realign(a):
    from . import provenance, realign, star
    i, n = (int(x) for x in a.shard.split('/'))
    scores = star.Scores.parse(a.scores) if a.scores is not None else star.SCORES
    summary = realign.run(a.packages, a.out, select=a.select, shard=(i, n), jobs=a.jobs, abpoa=a.abpoa,
                          timeout=a.timeout, max_cells=a.max_cells, scores=scores)
    provenance.write_manifest(a.out, 'realign', {'packages': a.packages},
                              {'select': a.select, 'shard': a.shard, 'timeout': a.timeout,
                               'max_cells': a.max_cells, 'method': 'medoid_star', 'scores': ' '.join(scores.flags())},
                              tools=[a.abpoa], extra={'summary': summary})
    json.dump(summary, sys.stdout, indent=1)
    print()


def cmd_qc(a):
    from . import provenance, qc, star
    scores = star.Scores.parse(a.scores) if a.scores is not None else star.SCORES
    summary = qc.run(a.extract, a.msas, a.out, pairs=a.pairs, refs=a.refs, sample=a.sample, jobs=a.jobs,
                     abpoa=a.abpoa, scores=scores, ref_sample=a.ref_sample)
    provenance.write_manifest(a.out, 'qc', {'msas': a.msas}, {'extract': a.extract, 'pairs': a.pairs,
                              'refs': a.refs, 'sample': a.sample, 'scores': ' '.join(scores.flags())},
                              tools=[a.abpoa], extra={'summary': summary})
    json.dump(summary, sys.stdout, indent=1)
    print()


def cmd_replace(a):
    from . import provenance, replace
    summary = replace.replace(a.graph, a.extract, a.msas, a.out, a.id_start, vg=a.vg, threads=a.threads,
                              id_mode=a.id_mode)
    if a.gbz:
        replace.to_gbz(a.out, a.gbz, vg=a.vg, keep_gfa=a.keep_gfa)
        summary['gbz'] = a.gbz
    provenance.write_manifest(a.gbz or a.out, 'replace', {'graph': a.graph, 'msas': a.msas},
                              {'id_start': a.id_start}, tools=['vg'], extra={'summary': summary})
    json.dump(summary, sys.stdout, indent=1)
    print()


def cmd_run(a):
    from . import pipeline, star
    scores = star.Scores.parse(a.scores) if a.scores is not None else star.SCORES
    drop = a.drop_sample if a.drop_sample is not None else ['gref_CHM13']
    summary = pipeline.run_contig(a.gbz, a.contig, a.targets, a.out, jobs=a.jobs, threads=a.threads, vg=a.vg,
                                  abpoa=a.abpoa, ref_sample=a.ref_sample, drop_samples=drop, scores=scores,
                                  timeout=a.timeout, max_cells=a.max_cells, id_start=a.id_start, do_qc=not a.no_qc,
                                  log=lambda m: print(m, file=sys.stderr))
    json.dump(summary, sys.stdout, indent=1)
    print()


def cmd_merge(a):
    from . import merge, provenance
    summary = merge.merge(a.gbz, a.out, vg=a.vg)
    provenance.write_manifest(a.out, 'merge', {'gbz_%d' % i: g for i, g in enumerate(a.gbz)}, {}, tools=['vg'],
                              extra={'summary': summary})
    json.dump(summary, sys.stdout, indent=1)
    print()


def cmd_union(a):
    from . import provenance, union
    summary = union.union_packages(a.base, a.add, a.out)
    provenance.write_manifest(a.out, 'union', {'base': a.base, 'add': a.add}, {}, extra={'summary': summary})
    json.dump(summary, sys.stdout, indent=1)
    print()


def cmd_project(a):
    from . import provenance, union
    summary = union.project(a.union, a.msas, a.out)
    provenance.write_manifest(a.out, 'project', {'union': a.union, 'msas': a.msas}, {}, extra={'summary': summary})
    json.dump(summary, sys.stdout, indent=1)
    print()


def cmd_verify(a):
    from . import provenance, verify
    res = verify.verify(a.before, a.after, vg=a.vg)
    if a.out:
        with open(a.out, 'w') as f:
            json.dump(res, f, indent=1)
        provenance.write_manifest(a.out, 'verify', {'before': a.before, 'after': a.after}, {}, tools=['vg'])
    json.dump(res, sys.stdout, indent=1)
    print()
    if not res['ok']:
        sys.exit(1)


def main(argv=None):
    from .star import SCORES as star_scores
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

    c = sub.add_parser('prepare', help="cut one contig from a whole-genome GBZ; snarls; reference index")
    c.add_argument('--gbz', required=True, help='whole-genome GBZ')
    c.add_argument('--contig', required=True)
    c.add_argument('-o', '--out', required=True, help='output directory')
    c.add_argument('--ref-sample', default='CHM13')
    c.add_argument('--drop-sample', action='append', default=None,
                   help='sample to remove before realigning (repeatable) [gref_CHM13]')
    c.add_argument('--vg', help='vg binary [vg on PATH]')
    c.add_argument('-t', '--threads', type=int, default=4)
    c.set_defaults(func=cmd_prepare)

    c = sub.add_parser('regions', help='anchor every target in the snarl tree; merge; size policy (stage 1)')
    c.add_argument('--targets', required=True, help='catalogue or BED of target intervals')
    c.add_argument('--contig', required=True)
    c.add_argument('--snarls', required=True, help='<contig>.snarls.pb from prepare')
    c.add_argument('--ref', required=True, help='<contig>.ref.tsv.gz from prepare')
    c.add_argument('-o', '--out', required=True, help='regions.tsv(.gz)')
    c.add_argument('--pad', type=int, default=200, help='pad each target by this many bp [200]')
    c.add_argument('--max-span', type=int, default=100000,
                   help='keep regions whose anchor span exceeds this as they are [100000]')
    c.add_argument('--vg', help='vg binary [vg on PATH]')
    c.set_defaults(func=cmd_regions)

    c = sub.add_parser('extract', help="every walk's run through every region, one streaming pass (stage 2)")
    c.add_argument('--graph', required=True, help='<contig>.gbz from prepare, or a GFA')
    c.add_argument('--regions', required=True)
    c.add_argument('-o', '--out', required=True, help='output directory')
    c.add_argument('--ref-sample', default='CHM13', help='the reference sample, whose allele each package names [CHM13]')
    c.add_argument('--vg', help='vg binary [vg on PATH]')
    c.add_argument('-t', '--threads', type=int, default=4)
    c.set_defaults(func=cmd_extract)

    c = sub.add_parser('realign', help='one MSA per selected region (stage 3; default the medoid star)')
    c.add_argument('--packages', required=True, help='packages.jsonl.gz from extract')
    c.add_argument('-o', '--out', required=True, help='msas.jsonl.gz')
    c.add_argument('--select', choices=['any', 'length'], default='any',
                   help="realign a region when its alleles differ at all ('any') or in length [any]")
    c.add_argument('--shard', default='0/1', help='do shard i of n (by package order) [0/1]')
    c.add_argument('-j', '--jobs', type=int, default=1)
    c.add_argument('--abpoa', default='abpoa')
    c.add_argument('--timeout', type=float, default=900, help='seconds per region [900]')
    c.add_argument('--max-cells', type=int, default=450_000_000,
                   help='keep a region whose centre x longest allele exceeds this [450,000,000]')
    c.add_argument('--scores', metavar='FLAGS',
                   help="abPOA's scores for every alignment, as its flags -M, -X, -O and -E; a flag left out "
                        "takes abPOA's default, so '' gives abPOA's defaults [%s]" % ' '.join(star_scores.flags()))
    c.set_defaults(func=cmd_realign)

    c = sub.add_parser('qc', help="how well each region's alignment agrees with the pairwise alignments it implies")
    c.add_argument('--extract', required=True, help='the extract output directory')
    c.add_argument('--msas', required=True, help='msas.jsonl.gz from realign')
    c.add_argument('-o', '--out', required=True, help='per-region TSV; the summary goes to OUT.json')
    c.add_argument('--pairs', type=int, default=20, help='random allele pairs per region [20]')
    c.add_argument('--refs', type=int, default=20, help='heaviest alleles compared with the reference allele [20]')
    c.add_argument('--sample', type=int, help='measure this many random regions rather than all')
    c.add_argument('--ref-sample', default='CHM13', help='the reference sample, for packages that do not name it [CHM13]')
    c.add_argument('--scores', metavar='FLAGS', help="abPOA's scores for the optimal pairwise alignments, as for realign; "
                   "use the scores the alignment was built with [the realign default]")
    c.add_argument('-j', '--jobs', type=int, default=1)
    c.add_argument('--abpoa', default='abpoa')
    c.set_defaults(func=cmd_qc)

    c = sub.add_parser('replace', help="splice the realigned regions into the contig's graph (stage 5)")
    c.add_argument('--graph', required=True, help='<contig>.gbz from prepare, or a GFA')
    c.add_argument('--extract', required=True, help='the extract output directory')
    c.add_argument('--msas', required=True, help='msas.jsonl.gz from realign')
    c.add_argument('-o', '--out', required=True, help='output GFA')
    c.add_argument('--id-mode', choices=['dense', 'reuse'], default='dense',
                   help="'dense' renumbers every node, new ones next to their regions; 'reuse' keeps old ids "
                        "and puts overflow at --id-start [dense]")
    c.add_argument('--id-start', type=int, help="dense: the first id [the contig's smallest id]; reuse: the first "
                   "id for new nodes that the region's freed ids do not cover (required)")
    c.add_argument('--gbz', help='also build this GBZ from the output GFA (vg gbwt -G)')
    c.add_argument('--keep-gfa', action='store_true', help='keep the GFA after --gbz')
    c.add_argument('--vg', help='vg binary [vg on PATH]')
    c.add_argument('-t', '--threads', type=int, default=4)
    c.set_defaults(func=cmd_replace)

    c = sub.add_parser('run', help='every stage for one contig, resuming where an earlier run stopped')
    c.add_argument('--gbz', required=True, help='whole-genome GBZ')
    c.add_argument('--contig', required=True)
    c.add_argument('--targets', required=True, help='catalogue or BED of target intervals')
    c.add_argument('-o', '--out', required=True, help='output directory; the contig goes in OUT/CONTIG')
    c.add_argument('-j', '--jobs', type=int, default=1, help='parallel regions in realign and qc [1]')
    c.add_argument('-t', '--threads', type=int, default=4, help='threads for vg [4]')
    c.add_argument('--ref-sample', default='CHM13')
    c.add_argument('--drop-sample', action='append', default=None,
                   help='sample to remove before realigning (repeatable) [gref_CHM13]')
    c.add_argument('--scores', metavar='FLAGS', help="abPOA's scores, as for realign")
    c.add_argument('--timeout', type=float, default=900, help='seconds per region [900]')
    c.add_argument('--max-cells', type=int, default=450_000_000, help='as for realign [450,000,000]')
    c.add_argument('--id-start', type=int, help="replace's first node id [the contig's smallest id]")
    c.add_argument('--no-qc', action='store_true', help='skip the qc stage')
    c.add_argument('--abpoa', default='abpoa')
    c.add_argument('--vg', help='vg binary [vg on PATH]')
    c.set_defaults(func=cmd_run)

    c = sub.add_parser('merge', help='one whole-genome GBZ from per-contig GBZs whose node ids do not overlap')
    c.add_argument('gbz', nargs='+', help='the per-contig GBZs: realigned contigs and the parts left as they are')
    c.add_argument('-o', '--out', required=True, help='the merged GBZ')
    c.add_argument('--vg', help='vg binary [vg on PATH]')
    c.set_defaults(func=cmd_merge)

    c = sub.add_parser('union', help="merge a panel's packages with a target graph's, region by region")
    c.add_argument('--base', required=True, help="the panel graph's packages.jsonl.gz")
    c.add_argument('--add', required=True, help="the target graph's packages.jsonl.gz (same regions)")
    c.add_argument('-o', '--out', required=True, help='union packages.jsonl.gz (realign this)')
    c.set_defaults(func=cmd_union)

    c = sub.add_parser('project', help="the union's MSAs restricted to the target's members")
    c.add_argument('--union', required=True, help='union packages.jsonl.gz')
    c.add_argument('--msas', required=True, help="realign's msas.jsonl.gz for the union")
    c.add_argument('-o', '--out', required=True, help="msas.jsonl.gz for the target's replace")
    c.set_defaults(func=cmd_project)

    c = sub.add_parser('verify', help='every path of the realigned graph spells what it did before (stage 6)')
    c.add_argument('--before', required=True, help='the original contig GBZ')
    c.add_argument('--after', required=True, help='the realigned contig GBZ')
    c.add_argument('-o', '--out', help='write the report here too')
    c.add_argument('--vg', help='vg binary [vg on PATH]')
    c.set_defaults(func=cmd_verify)

    a = p.parse_args(argv)
    if getattr(a, 'drop_sample', 'unset') is None:
        a.drop_sample = ['gref_CHM13']
    a.func(a)


if __name__ == '__main__':
    main()
