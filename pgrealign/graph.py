"""Reading a graph as a stream of GFA lines, and the walk and reference-path helpers the stages share.

The stages never hold a whole graph's walks in memory. They read `vg convert -f <graph.gbz>` (or a GFA
file) once, line by line:
- S lines: segments (nodes), stored as a dict of sequences, or of lengths only;
- L lines: edges, ignored unless asked for;
- W lines: walks, one per haplotype fragment or reference path.

A walk's identity is its W-line metadata, sample#haplotype#contig#start, which is unique in a GBZ.
"""
import collections
import gzip
import re
import shutil
import subprocess

STEP = re.compile(r'([<>])(\d+)')
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

Walk = collections.namedtuple('Walk', 'sample hap contig start end steps')


def revcomp(s):
    return s.translate(COMP)[::-1]


def walk_key(w):
    """The unique name of a W line: sample#hap#contig#start."""
    return '%s#%s#%s#%s' % (w.sample, w.hap, w.contig, w.start)


def find_vg(vg=None):
    vg = vg or shutil.which('vg')
    if vg is None:
        raise RuntimeError('vg is not on PATH (or pass --vg)')
    return vg


class GfaStream:
    """Iterate over the lines of a GFA file or of `vg convert -f GBZ`, as a context manager."""

    def __init__(self, path, vg=None, threads=4):
        self.path, self.vg, self.threads, self.proc, self.fh = path, vg, threads, None, None

    def __enter__(self):
        if self.path.endswith('.gbz'):
            cmd = [find_vg(self.vg), 'convert', '-f', '-t', str(self.threads), self.path]
            self.proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
            self.fh = self.proc.stdout
        else:
            self.fh = gzip.open(self.path, 'rt') if self.path.endswith('.gz') else open(self.path)
        return self.fh

    def __exit__(self, *exc):
        if self.proc is not None:
            self.proc.stdout.close()
            rc = self.proc.wait()
            if rc and exc[0] is None:
                raise RuntimeError('vg convert -f %s failed (%d)' % (self.path, rc))
        elif self.fh is not None:
            self.fh.close()
        return False


def parse_w(line):
    """A W line as a Walk whose `steps` is still the raw walk string."""
    x = line.rstrip('\n').split('\t')
    return Walk(x[1], x[2], x[3], x[4], x[5], x[6])


def steps(walk_string):
    """[(orient, node)] of a walk string such as '>12<13>14'."""
    return [(o, int(n)) for o, n in STEP.findall(walk_string)]


def spell(path_steps, seq):
    """The sequence a list of (orient, node) steps spells, given node sequences."""
    return ''.join(seq[n] if o == '>' else revcomp(seq[n]) for o, n in path_steps)


def reference_walk(stream_path, sample, contig, vg=None):
    """[(node, orient, start, end)] of the reference W line `sample`/`contig`, 0-based, plus node lengths
    for every segment. Reads the whole stream once."""
    lens, ref = {}, None
    with GfaStream(stream_path, vg) as f:
        for line in f:
            t = line[0]
            if t == 'S':
                x = line.split('\t', 3)
                lens[int(x[1])] = len(x[2].rstrip('\n'))
            elif t == 'W':
                head = line.split('\t', 4)
                if head[1] == sample and head[3] == contig:
                    if ref is not None:
                        raise ValueError('more than one %s W line for %s' % (sample, contig))
                    ref = parse_w(line)
    if ref is None:
        raise ValueError('no %s W line for %s in %s' % (sample, contig, stream_path))
    pos, out = int(ref.start), []
    for o, n in steps(ref.steps):
        out.append((n, '+' if o == '>' else '-', pos, pos + lens[n]))
        pos += lens[n]
    return out
