"""Stage 6, `verify`: the realigned graph spells exactly what the original did.

The check compares two GBZs (or GFAs):
- The path metadata (`vg paths -M`) must be identical.
- Every path's full sequence (`vg paths -F`) must have the same length and MD5.
- The graph sizes before and after are reported.

`replace` already checks every run it rewrites, and walks outside the regions are copied unchanged.
So this is the end-to-end confirmation, and it streams every base of every path once per graph.
"""
import hashlib
import subprocess

from .graph import find_vg


def path_digests(graph, vg=None):
    """{path name: (length, md5)} from `vg paths -F`."""
    vg = find_vg(vg)
    out, name, h, n = {}, None, None, 0
    p = subprocess.Popen([vg, 'paths', '-F', '-x', graph], stdout=subprocess.PIPE, bufsize=1 << 20)
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
    if p.wait():
        raise RuntimeError('vg paths -F %s failed' % graph)
    return out


def metadata(graph, vg=None):
    p = subprocess.run([find_vg(vg), 'paths', '-M', '-x', graph], capture_output=True, text=True, check=True)
    return sorted(p.stdout.splitlines())


def stats(graph, vg=None):
    p = subprocess.run([find_vg(vg), 'stats', '-lz', graph], capture_output=True, text=True, check=True)
    return dict((k, int(v)) for k, v in (l.split('\t') for l in p.stdout.splitlines() if '\t' in l))


def verify(before, after, vg=None):
    res = {'before': stats(before, vg), 'after': stats(after, vg)}
    res['metadata_identical'] = metadata(before, vg) == metadata(after, vg)
    a, b = path_digests(before, vg), path_digests(after, vg)
    res['paths'] = len(a)
    res['paths_missing'] = sorted(set(a) - set(b))[:20]
    res['paths_extra'] = sorted(set(b) - set(a))[:20]
    res['paths_differ'] = sorted(k for k in set(a) & set(b) if a[k] != b[k])[:20]
    res['bases'] = sum(v[0] for v in a.values())
    res['ok'] = (res['metadata_identical'] and not res['paths_missing'] and not res['paths_extra']
                 and not res['paths_differ'])
    return res
