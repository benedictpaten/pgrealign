"""`run`: every stage for one contig, in order, resuming where an earlier run stopped.

For a contig C of a whole-genome GBZ, `run_contig` writes, under OUT/C/:

| stage | output | skipped when |
|---|---|---|
| prepare | C.gbz, C.snarls.pb, C.ref.tsv.gz | each exists |
| regions | regions.tsv.gz | it exists |
| extract | extract/ | extract/summary.json exists |
| realign | msas.jsonl.gz | it exists; an interrupted realign resumes from msas.jsonl.gz.partial |
| replace | C.realigned.gbz, and replace.json | the GBZ exists |
| verify | verify.json | it exists |
| qc | qc.tsv and qc.tsv.json | it exists, or qc is off |

Each stage also writes its manifest, as the subcommands do. The realigned contigs of a genome are then
combined with `pgrealign merge`.
"""
import json
import os

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
               drop_samples=('gref_CHM13',), scores=SCORES, timeout=900, max_cells=realign.MAX_CELLS,
               id_start=None, do_qc=True, log=print):
    """Run every stage for `contig` (see the module docstring). Returns {stage: summary}."""
    d = os.path.join(out_dir, contig)
    os.makedirs(d, exist_ok=True)
    out = {}

    log('[%s] prepare' % contig)
    out['prepare'] = prepare.prepare(gbz, contig, d, ref_sample=ref_sample, drop_samples=drop_samples, vg=vg,
                                     threads=threads)
    cgbz, snarls, ref = out['prepare']['gbz'], out['prepare']['snarls'], out['prepare']['ref']

    reg = os.path.join(d, 'regions.tsv.gz')
    if not _done(reg):
        log('[%s] regions' % contig)
        s = regions.build(targets, contig, snarls, ref, reg, vg=vg)
        provenance.write_manifest(reg, 'regions', {'targets': targets, 'snarls': snarls, 'ref': ref},
                                  {'contig': contig}, tools=['vg'], extra={'summary': s})
        out['regions'] = s

    ex = os.path.join(d, 'extract')
    if not _done(os.path.join(ex, 'summary.json')):
        log('[%s] extract' % contig)
        out['extract'] = extract.extract(cgbz, [r for r in regions.read_regions(reg)], ex, vg=vg,
                                         threads=threads, ref_sample=ref_sample)

    msas = os.path.join(d, 'msas.jsonl.gz')
    if not _done(msas):
        log('[%s] realign' % contig)
        s = realign.run(os.path.join(ex, 'packages.jsonl.gz'), msas, jobs=jobs, abpoa=abpoa, timeout=timeout,
                        max_cells=max_cells, scores=scores)
        provenance.write_manifest(msas, 'realign', {'packages': os.path.join(ex, 'packages.jsonl.gz')},
                                  {'method': 'medoid_star', 'scores': ' '.join(scores.flags()),
                                   'timeout': timeout, 'max_cells': max_cells}, tools=[abpoa], extra={'summary': s})
        out['realign'] = s

    rgbz = os.path.join(d, contig + '.realigned.gbz')
    if not _done(rgbz):
        log('[%s] replace' % contig)
        gfa = os.path.join(d, contig + '.realigned.gfa')
        s = replace.replace(cgbz, ex, msas, gfa, id_start=id_start, vg=vg, threads=threads)
        replace.to_gbz(gfa, rgbz, vg=vg)
        _write_json(os.path.join(d, 'replace.json'), s)
        provenance.write_manifest(rgbz, 'replace', {'graph': cgbz, 'msas': msas}, {'id_start': id_start},
                                  tools=['vg'], extra={'summary': s})
        out['replace'] = s

    ver = os.path.join(d, 'verify.json')
    if not _done(ver):
        log('[%s] verify' % contig)
        s = verify.verify(cgbz, rgbz, vg=vg)
        _write_json(ver, s)
        if not s['ok']:
            raise RuntimeError('verify failed for %s: see %s' % (contig, ver))
        out['verify'] = {k: s[k] for k in ('ok', 'paths', 'bases', 'metadata_identical')}

    q = os.path.join(d, 'qc.tsv')
    if do_qc and not _done(q):
        log('[%s] qc' % contig)
        out['qc'] = qc.run(ex, msas, q, jobs=jobs, abpoa=abpoa, scores=scores, ref_sample=ref_sample)

    out['id_range'] = list(merge.node_range(rgbz, find_vg(vg)))
    return out
