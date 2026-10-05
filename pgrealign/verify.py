"""Stage 6, `verify`: every path of the realigned graph spells exactly what it spelled before.

`replace` gives every run through a replaced region new steps, and under its default dense numbering
every walk's node ids change, so walks cannot be compared step by step. `verify` compares two graphs
(GBZ or GFA) by what their paths spell and how they are named:
- The path metadata (`vg paths -M`: name, sense, sample, haplotype, locus, phase block and subrange)
  must be identical.
- Every path's sequence (`vg paths -F`), upper-cased, must have the same length and MD5 in both.
- The node and edge counts and sequence length of both graphs (`vg stats -lz`) are reported.

The report counts the paths missing from `after`, the paths extra in it, and the paths in both whose
sequences differ (`paths_missing`, `paths_extra`, `paths_differ`), and names up to SHOW of each, the
first in sorted order. It is `ok` when the metadata is identical and all three counts are 0.

`replace` checks each run as it rewrites it. This is the end-to-end check, on the GBZ built from its
output, and it streams every base of every path once per graph.
"""
import hashlib
import subprocess

from .graph import find_vg

SHOW = 20


def path_digests(graph, vg=None):
    """{path name: (length, md5)} from `vg paths -F`."""
    vg = find_vg(vg)
    out, name, h, n = {}, None, None, 0
    with subprocess.Popen([vg, 'paths', '-F', '-x', graph], stdout=subprocess.PIPE, bufsize=1 << 20) as p:
        for line in p.stdout:
            if line.startswith(b'>'):
                if name is not None:
                    out[name] = (n, h.hexdigest())
                name, h, n = line[1:].split()[0].decode(), hashlib.md5(), 0
            else:
                s = line.strip().upper()
                h.update(s)
                n += len(s)
        if name is not None:
            out[name] = (n, h.hexdigest())
    if p.returncode:
        raise RuntimeError('vg paths -F %s failed' % graph)
    return out


def metadata(graph, vg=None):
    p = subprocess.run([find_vg(vg), 'paths', '-M', '-x', graph], capture_output=True, text=True, check=True)
    return sorted(p.stdout.splitlines())


def stats(graph, vg=None):
    p = subprocess.run([find_vg(vg), 'stats', '-lz', graph], capture_output=True, text=True, check=True)
    return dict((k, int(v)) for k, v in (l.split('\t') for l in p.stdout.splitlines() if '\t' in l))


def verify(before, after, vg=None):
    """The report described in the module docstring, for `after` against `before`."""
    res = {'before': stats(before, vg), 'after': stats(after, vg)}
    res['metadata_identical'] = metadata(before, vg) == metadata(after, vg)
    a, b = path_digests(before, vg), path_digests(after, vg)
    res['paths'] = len(a)
    res['bases'] = sum(v[0] for v in a.values())
    diff = {'missing': set(a) - set(b), 'extra': set(b) - set(a),
            'differ': {k for k in set(a) & set(b) if a[k] != b[k]}}
    for kind, names in diff.items():
        res['paths_' + kind] = len(names)
        res['paths_%s_names' % kind] = sorted(names)[:SHOW]
    res['ok'] = res['metadata_identical'] and not any(diff.values())
    return res
