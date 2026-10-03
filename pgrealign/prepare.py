"""Stage `prepare` (per contig): the inputs the later stages need, cut from a whole-genome GBZ.

1. `<contig>.gbz`: the contig's connected component (`vg chunk --gbz --contig`). Node ids stay
   genome-wide.
2. The named samples (e.g. the gref_CHM13 cover) are removed from it. Each of their paths starts and
   ends inside regions, which would make every region a fragment case, so the cover is recomputed
   after `replace` instead (`vg paths -u`).
3. `<contig>.snarls.pb`: `vg snarls -T -P <ref>#0#<contig>`, the decomposition vg call itself uses.
4. `<contig>.ref.tsv.gz`: the reference path's nodes, as node, orient, start and end (0-based).
"""
import gzip
import os
import shutil
import subprocess
import tempfile

from .graph import find_vg, steps


def run(cmd, **kw):
    subprocess.run([str(c) for c in cmd], check=True, **kw)


def prepare(gbz, contig, out_dir, ref_sample='CHM13', drop_samples=('gref_CHM13',), vg=None, threads=4):
    vg = find_vg(vg)
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, contig + '.gbz')
    tmp = tempfile.mkdtemp(prefix='prepare.', dir=out_dir)
    try:
        if not os.path.exists(out):
            run([vg, 'chunk', '-x', gbz, '--gbz', '--contig', contig, '-b', os.path.join(tmp, 'chunk'),
                 '-t', threads])
            chunks = [f for f in os.listdir(tmp) if f.startswith('chunk_') and f.endswith('.gbz')]
            if len(chunks) != 1:
                raise RuntimeError('vg chunk --contig %s wrote %d components' % (contig, len(chunks)))
            chunk = os.path.join(tmp, chunks[0])
            if drop_samples:
                g1, g2 = os.path.join(tmp, 'a.gbwt'), os.path.join(tmp, 'b.gbwt')
                run([vg, 'gbwt', '-Z', chunk, '-o', g1])
                run([vg, 'gbwt', '-o', g2, *sum([['-R', s] for s in drop_samples], []), g1])
                run([vg, 'gbwt', '-x', chunk, g2, '--gbz-format', '-g', out + '.tmp'])
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
            write_ref_index(vg, out, ref_sample, contig, ref, tmp, threads)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return {'gbz': out, 'snarls': snarls, 'ref': ref}


def write_ref_index(vg, gbz, ref_sample, contig, out, tmp, threads):
    """The reference walk with node positions, from `vg convert -f -H` (no haplotype paths, so it does
    not stream every haplotype)."""
    lens, walk, start = {}, None, 0
    p = subprocess.Popen([vg, 'convert', '-f', '-H', '-t', str(threads), gbz], stdout=subprocess.PIPE,
                         text=True, bufsize=1 << 20)
    for line in p.stdout:
        t = line[0]
        if t == 'S':
            x = line.split('\t', 3)
            lens[int(x[1])] = len(x[2].rstrip('\n'))
        elif t == 'W':
            x = line.rstrip('\n').split('\t')
            if x[1] == ref_sample and x[3] == contig:
                if walk is not None:
                    raise RuntimeError('more than one %s W line for %s' % (ref_sample, contig))
                walk, start = x[6], int(x[4])
        elif t == 'P':
            x = line.rstrip('\n').split('\t')
            if x[1] in ('%s#0#%s' % (ref_sample, contig), '%s#%s' % (ref_sample, contig)):
                walk = ''.join(('>' if s[-1] == '+' else '<') + s[:-1] for s in x[2].split(','))
    if p.wait():
        raise RuntimeError('vg convert -f -H %s failed' % gbz)
    if walk is None:
        raise RuntimeError('no %s path for %s in %s' % (ref_sample, contig, gbz))
    pos = start
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
