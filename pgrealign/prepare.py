"""Stage `prepare` (per contig): the inputs the later stages need, cut from a whole-genome GBZ.

1. `<contig>.gbz`: the contig's connected component (`vg chunk --gbz --contig`). Node ids stay
   genome-wide.
2. The samples in `drop_samples` (by default DROP_SAMPLES, the gref_CHM13 cover) are removed from it
   with `vg gbwt -R`. Each of their paths starts and ends inside regions, which would make every
   region a fragment case. The GBZ's reference samples are then set to those that remain, since
   `vg gbwt -R` leaves the reference_samples tag naming the removed ones. The reference sample itself
   cannot be dropped.
3. `<contig>.snarls.pb`: `vg snarls -T -P <ref>#0#<contig>`, the decomposition vg call itself uses.
4. `<contig>.ref.tsv.gz`: the reference path's nodes, as node, orient, start and end (0-based).

With the default DROP_SAMPLES, the output graph, and so the realigned graph built from it, has no gref
cover. To give the realigned graph one, compute it on a mutable copy and rebuild the GBZ from that:

    vg convert -p <contig>.realigned.gbz > g.pg
    vg paths -x g.pg -u -Q <ref_sample> > g.gref.pg
    vg convert -f g.gref.pg > g.gref.gfa
    vg gbwt -G g.gref.gfa --gbz-format -g <contig>.realigned.gref.gbz
"""
import gzip
import os
import shutil
import subprocess
import tempfile

from .graph import GfaStream, find_vg, parse_w, steps

DROP_SAMPLES = ('gref_CHM13',)


def run(cmd, **kw):
    return subprocess.run([str(c) for c in cmd], check=True, **kw)


def prepare(gbz, contig, out_dir, ref_sample='CHM13', drop_samples=None, vg=None, threads=4):
    """Write the module docstring's outputs for `contig` of `gbz` to `out_dir`, skipping each that
    exists. `drop_samples` None means DROP_SAMPLES; an empty sequence drops nothing. Returns their paths."""
    drop_samples = tuple(DROP_SAMPLES if drop_samples is None else drop_samples)
    if ref_sample in drop_samples:
        raise ValueError('the reference sample %s cannot be dropped' % ref_sample)
    vg = find_vg(vg)
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, contig + '.gbz')
    tmp = tempfile.mkdtemp(prefix='prepare.', dir=out_dir)
    try:
        if not os.path.exists(out):
            run([vg, 'chunk', '-x', gbz, '--gbz', '--contig', contig, '-b', os.path.join(tmp, 'chunk')])
            chunks = [f for f in os.listdir(tmp) if f.startswith('chunk_') and f.endswith('.gbz')]
            if len(chunks) != 1:
                raise RuntimeError('vg chunk --contig %s wrote %d components' % (contig, len(chunks)))
            chunk = os.path.join(tmp, chunks[0])
            if drop_samples:
                drop_samples_from(chunk, drop_samples, out + '.tmp', tmp, vg)
                os.replace(out + '.tmp', out)
            else:
                os.replace(chunk, out)
        ref_path = '%s#0#%s' % (ref_sample, contig)
        snarls = os.path.join(out_dir, contig + '.snarls.pb')
        if not os.path.exists(snarls):
            with open(snarls + '.tmp', 'wb') as f:
                run([vg, 'snarls', '-T', '-P', ref_path, '-t', threads, out], stdout=f)
            os.replace(snarls + '.tmp', snarls)
        ref = os.path.join(out_dir, contig + '.ref.tsv.gz')
        if not os.path.exists(ref):
            write_ref_index(out, ref_sample, contig, ref, vg=vg, threads=threads)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return {'gbz': out, 'snarls': snarls, 'ref': ref}


def reference_samples(gbwt, vg):
    """The samples a GBWT's reference_samples tag names, in its order."""
    p = run([vg, 'gbwt', '--tags', gbwt], capture_output=True, text=True)
    for line in p.stdout.splitlines():
        k, _, v = line.partition('\t')
        if k == 'reference_samples':
            return v.split()
    return []


def drop_samples_from(gbz, samples, out, tmp, vg):
    """Write `gbz` without the paths of `samples` to `out`, with the reference samples that remain.
    `tmp` is a directory for the GBWTs in between."""
    g1, g2 = os.path.join(tmp, 'in.gbwt'), os.path.join(tmp, 'out.gbwt')
    run([vg, 'gbwt', '-Z', gbz, '-o', g1])
    refs = reference_samples(g1, vg)
    keep = [s for s in refs if s not in samples]
    if refs and not keep:
        raise ValueError('dropping %s would leave no reference sample' % ', '.join(samples))
    cmd = [vg, 'gbwt', '-o', g2]
    for s in samples:
        cmd += ['-R', s]
    for s in keep:
        cmd += ['--set-reference', s]
    run(cmd + [g1])
    run([vg, 'gbwt', '-x', gbz, g2, '--gbz-format', '-g', out])


def write_ref_index(graph, ref_sample, contig, out, vg=None, threads=4):
    """Write the reference path's nodes with their 0-based positions to `out`. The path is the W line of
    `ref_sample` on `contig`, or a P line named <ref_sample>#0#<contig> or <ref_sample>#<contig>. A GBZ
    is read with `vg convert -f -H`, which leaves out the haplotype walks."""
    names = ('%s#0#%s' % (ref_sample, contig), '%s#%s' % (ref_sample, contig))
    lens, found = {}, []
    with GfaStream(graph, vg, threads, drop_haplotypes=True) as f:
        for line in f:
            t = line[0]
            if t == 'S':
                x = line.split('\t', 3)
                lens[int(x[1])] = len(x[2].rstrip('\n'))
            elif t == 'W':
                w = parse_w(line)
                if w.sample == ref_sample and w.contig == contig:
                    found.append((int(w.start), w.steps))
            elif t == 'P':
                x = line.rstrip('\n').split('\t')
                if x[1] in names:
                    found.append((0, ''.join(('>' if s[-1] == '+' else '<') + s[:-1] for s in x[2].split(','))))
    if len(found) != 1:
        raise RuntimeError('%d %s paths for %s in %s, not 1' % (len(found), ref_sample, contig, graph))
    pos, walk = found[0]
    with gzip.open(out + '.tmp', 'wt') as f:
        f.write('#node\torient\tstart\tend\n')
        for o, n in steps(walk):
            f.write('%d\t%s\t%d\t%d\n' % (n, '+' if o == '>' else '-', pos, pos + lens[n]))
            pos += lens[n]
    os.replace(out + '.tmp', out)


def read_ref_index(path):
    with gzip.open(path, 'rt') as f:
        return [(int(x[0]), x[1], int(x[2]), int(x[3])) for x in (l.rstrip('\n').split('\t') for l in f)
                if not x[0].startswith('#')]
