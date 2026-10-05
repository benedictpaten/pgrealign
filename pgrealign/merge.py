"""Stage `merge`: one whole-genome GBZ from per-contig GBZs.

The inputs are GBZs whose node ids do not overlap: each realigned contig from replace, and every part of
the genome that was not realigned (contigs without targets, such as chrM, and components with no
reference contig), cut out as they are. Nothing is renumbered.
1. Each input's node-id range is read with `vg stats -r`. If two ranges overlap, merge stops and names
   them. replace's dense numbering starts at the contig's smallest input id, so a contig whose graph
   grew past the ids its input left free runs into the next contig; replace that contig again with
   `--id-start` above every other input's ids.
2. Each input's GBWT is taken out with `vg gbwt -Z`.
3. The inputs' nodes and edges are concatenated into one GFA of S and L lines, from
   `vg convert -f -H`.
4. `vg gbwt -f` merges the GBWTs, which needs disjoint node ids, and builds the GBZ over that GFA. The
   reference samples are those with reference paths in any input; `vg gbwt` would otherwise drop the
   tag and read them as haplotypes.
5. The output is checked: its node count, edge count and length are the sums of the inputs', and its
   path names are exactly the union of theirs. If not, it is renamed to <output>.failed.
"""
import os
import shutil
import subprocess
import tempfile

from .graph import find_vg


def _run(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw).stdout


def node_range(gbz, vg):
    """(smallest, largest) node id of a GBZ."""
    out = _run([vg, 'stats', '-r', gbz])
    lo, hi = out.split()[-1].split(':')
    return int(lo), int(hi)


def graph_size(gbz, vg):
    """{'nodes', 'edges', 'length'} of a GBZ."""
    out = _run([vg, 'stats', '-z', '-l', gbz]).split()
    vals = dict(zip(out[0::2], out[1::2]))
    return {'nodes': int(vals['nodes']), 'edges': int(vals['edges']), 'length': int(vals['length'])}


def path_names(gbz, vg):
    """(set of path names, set of samples with a reference path) of a GBZ."""
    names, refs = set(), set()
    for line in _run([vg, 'paths', '-M', '-x', gbz]).splitlines():
        if line.startswith('#'):
            continue
        f = line.split('\t')
        names.add(f[0])
        if f[1] == 'REFERENCE':
            refs.add(f[2])
    return names, refs


def check_disjoint(ranges):
    """`ranges` is {input: (lo, hi)}. Raises ValueError naming the first two that overlap."""
    order = sorted(ranges.items(), key=lambda kv: kv[1])
    for (a, (_, a_hi)), (b, (b_lo, _)) in zip(order, order[1:]):
        if b_lo <= a_hi:
            raise ValueError('node ids overlap: %s ends at %d, %s starts at %d; replace one of them again '
                             'with --id-start above every other input\'s ids' % (a, a_hi, b, b_lo))


def merge(gbzs, out_gbz, vg=None, workdir=None):
    """Merge `gbzs` into `out_gbz` (see the module docstring). Returns a summary."""
    vg = find_vg(vg)
    ranges = {g: node_range(g, vg) for g in gbzs}
    check_disjoint(ranges)
    wd = tempfile.mkdtemp(prefix='merge.', dir=workdir)
    try:
        gbwts, names, refs, size = [], set(), set(), {'nodes': 0, 'edges': 0, 'length': 0}
        nodes_gfa = os.path.join(wd, 'nodes.gfa')
        with open(nodes_gfa, 'w') as out:
            out.write('H\tVN:Z:1.1\n')
            for i, g in enumerate(sorted(gbzs, key=lambda x: ranges[x])):
                gw = os.path.join(wd, 'in%d.gbwt' % i)
                _run([vg, 'gbwt', '-Z', g, '-o', gw])
                gbwts.append(gw)
                with subprocess.Popen([vg, 'convert', '-f', '-H', g], stdout=subprocess.PIPE, text=True,
                                      bufsize=1 << 20) as p:
                    for line in p.stdout:
                        if line[0] in 'SL':
                            out.write(line)
                if p.returncode:
                    raise subprocess.CalledProcessError(p.returncode, p.args)
                n, r = path_names(g, vg)
                names |= n
                refs |= r
                for k, v in graph_size(g, vg).items():
                    size[k] += v
        cmd = [vg, 'gbwt', '-x', nodes_gfa, '-f', '-g', out_gbz]
        for s in sorted(refs):
            cmd += ['--set-reference', s]
        _run(cmd + gbwts)
    finally:
        shutil.rmtree(wd, ignore_errors=True)
    got = graph_size(out_gbz, vg)
    got_names, got_refs = path_names(out_gbz, vg)
    problem = None
    if got != size:
        problem = 'merged graph has %s, the inputs sum to %s' % (got, size)
    elif got_names != names:
        problem = 'merged graph has %d paths, the inputs %d' % (len(got_names), len(names))
    if problem:
        os.replace(out_gbz, out_gbz + '.failed')
        raise ValueError(problem + '; the graph is left in %s.failed' % out_gbz)
    return {'inputs': len(gbzs), 'id_ranges': {g: list(r) for g, r in ranges.items()}, **got,
            'paths': len(got_names), 'reference_samples': sorted(got_refs)}
