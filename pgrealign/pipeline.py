"""`run`: every stage for one contig, in order, resuming where an earlier run stopped.

For a contig C of a whole-genome GBZ, `run_contig` writes, under OUT/C/:

| stage | output | skipped when |
|---|---|---|
| prepare | C.gbz, C.snarls.pb, C.ref.tsv.gz | each exists |
| regions | regions.tsv.gz | it exists |
| extract | extract/ | extract/summary.json exists |
| realign | msas.jsonl.gz | it exists; an interrupted realign resumes from msas.jsonl.gz.partial |
| replace | C.realigned.gbz, and replace.json | the GBZ exists |
| verify | verify.json | it exists; run stops if it records a failure |
| qc | qc.tsv and qc.tsv.json | it exists, or qc is off |

The file by which run judges a stage done is written after the stage's other files, under a temporary
name renamed into place once complete, so a run killed partway through a stage redoes that stage. Each
stage then writes its manifest, <output>.manifest.json, as the subcommands do. The realigned contigs of a genome are
then combined with `pgrealign merge`.
"""
import json
import os
import time

from . import extract, merge, prepare, provenance, qc, realign, regions, replace, verify
from .graph import find_vg
from .star import SCORES


def _done(path):
    return os.path.exists(path)


def _write_json(path, obj):
    with open(path + '.tmp', 'w') as f:
        json.dump(obj, f, indent=1)
        f.write('\n')
    os.replace(path + '.tmp', path)


def run_contig(gbz, contig, targets, out_dir, jobs=1, threads=4, vg=None, abpoa='abpoa', ref_sample='CHM13',
               drop_samples=None, scores=SCORES, timeout=900, max_cells=realign.MAX_CELLS, max_internal_bp=200000,
               id_start=None, do_qc=True, log=print):
    """Run every stage for `contig` (see the module docstring). `drop_samples` None means prepare's
    default. Returns {stage: summary} for the stages it ran, and the contig's node-id range."""
    d = os.path.join(out_dir, contig)
    os.makedirs(d, exist_ok=True)
    drop = prepare.DROP_SAMPLES if drop_samples is None else tuple(drop_samples)
    vg_tool, abpoa_tool = {'vg': vg}, {'abpoa': abpoa}
    out = {}

    log('[%s] prepare' % contig)
    t = time.time()
    out['prepare'] = prepare.prepare(gbz, contig, d, ref_sample=ref_sample, drop_samples=drop, vg=vg,
                                     threads=threads)
    cgbz, snarls, ref = out['prepare']['gbz'], out['prepare']['snarls'], out['prepare']['ref']
    if not _done(cgbz + '.manifest.json'):
        provenance.write_manifest(cgbz, 'prepare', {'gbz': gbz},
                                  {'contig': contig, 'ref_sample': ref_sample, 'drop_samples': list(drop)},
                                  tools=vg_tool, extra={'outputs': out['prepare']}, started=t)

    reg = os.path.join(d, 'regions.tsv.gz')
    if not _done(reg):
        log('[%s] regions' % contig)
        t = time.time()
        s = regions.build(targets, contig, snarls, ref, reg, vg=vg)
        provenance.write_manifest(reg, 'regions', {'targets': targets, 'snarls': snarls, 'ref': ref},
                                  {'contig': contig}, tools=vg_tool, extra={'summary': s}, started=t)
        out['regions'] = s

    ex = os.path.join(d, 'extract')
    if not _done(os.path.join(ex, 'summary.json')):
        log('[%s] extract' % contig)
        t = time.time()
        s = extract.extract(cgbz, list(regions.read_regions(reg)), ex, vg=vg, threads=threads,
                            ref_sample=ref_sample, max_internal_bp=max_internal_bp)
        provenance.write_manifest(os.path.join(ex, 'extract'), 'extract', {'graph': cgbz, 'regions': reg},
                                  {'ref_sample': ref_sample, 'max_internal_bp': max_internal_bp}, tools=vg_tool,
                                  extra={'summary': s}, started=t)
        out['extract'] = s

    msas = os.path.join(d, 'msas.jsonl.gz')
    if not _done(msas):
        log('[%s] realign' % contig)
        t = time.time()
        s = realign.run(os.path.join(ex, 'packages.jsonl.gz'), msas, jobs=jobs, abpoa=abpoa, timeout=timeout,
                        max_cells=max_cells, scores=scores)
        provenance.write_manifest(msas, 'realign', {'packages': os.path.join(ex, 'packages.jsonl.gz')},
                                  {'method': 'medoid_star', 'scores': ' '.join(scores.flags()),
                                   'timeout': timeout, 'max_cells': max_cells}, tools=abpoa_tool,
                                  extra={'summary': s}, started=t)
        out['realign'] = s

    rgbz = os.path.join(d, contig + '.realigned.gbz')
    if not _done(rgbz):
        log('[%s] replace' % contig)
        t = time.time()
        gfa = os.path.join(d, contig + '.realigned.gfa')
        s = replace.replace(cgbz, ex, msas, gfa, id_start=id_start, vg=vg, threads=threads)
        _write_json(os.path.join(d, 'replace.json'), s)
        replace.to_gbz(gfa, rgbz, vg=vg)
        provenance.write_manifest(rgbz, 'replace', {'graph': cgbz, 'msas': msas},
                                  {'extract': ex, 'id_mode': 'dense', 'id_start': id_start}, tools=vg_tool,
                                  extra={'summary': s}, started=t)
        out['replace'] = s

    ver = os.path.join(d, 'verify.json')
    if not _done(ver):
        log('[%s] verify' % contig)
        t = time.time()
        s = verify.verify(cgbz, rgbz, vg=vg)
        _write_json(ver, s)
        provenance.write_manifest(ver, 'verify', {'before': cgbz, 'after': rgbz}, {}, tools=vg_tool, started=t)
        out['verify'] = {k: s[k] for k in ('ok', 'paths', 'bases', 'metadata_identical')}
    with open(ver) as f:
        if not json.load(f)['ok']:
            raise RuntimeError('verify failed for %s: see %s' % (contig, ver))

    q = os.path.join(d, 'qc.tsv')
    if do_qc and not _done(q):
        log('[%s] qc' % contig)
        t = time.time()
        s = qc.run(ex, msas, q, jobs=jobs, abpoa=abpoa, scores=scores, ref_sample=ref_sample, max_cells=max_cells)
        provenance.write_manifest(q, 'qc', {'msas': msas}, {'extract': ex, 'scores': ' '.join(scores.flags())},
                                  tools=abpoa_tool, extra={'summary': s}, started=t)
        out['qc'] = s

    out['id_range'] = list(merge.node_range(rgbz, find_vg(vg)))
    return out
