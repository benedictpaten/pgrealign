"""Reading a graph as a stream of GFA lines, and the walk helpers the stages share.

The stages never hold a whole graph's walks in memory. They read `vg convert -f <graph>` of a graph in one
of vg's own formats (VG_SUFFIXES: a GBZ, or a .vg such as a Minigraph-Cactus chromosome graph), or a GFA
file, plain or gzipped, once, line by line:
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
# graphs read through `vg convert -f`; a file with any other name is read as GFA
VG_SUFFIXES = ('.gbz', '.vg', '.pg', '.hg', '.xg')

Walk = collections.namedtuple('Walk', 'sample hap contig start end steps')


def revcomp(s):
    return s.translate(COMP)[::-1]


def walk_key(w):
    """The unique name of a W line: sample#hap#contig#start."""
    return '%s#%s#%s#%s' % (w.sample, w.hap, w.contig, w.start)


def open_text(path, mode='rt'):
    """A text file, gzip-compressed when its name ends in .gz. `mode` is 'rt', 'wt' or 'at'."""
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


def find_vg(vg=None):
    vg = vg or shutil.which('vg')
    if vg is None:
        raise RuntimeError('vg is not on PATH (or pass --vg)')
    return vg


class GfaStream:
    """Iterate over the lines of a GFA file (plain or gzipped), or of `vg convert -f` of a graph whose name
    ends in one of VG_SUFFIXES, as a context manager. `drop_haplotypes` passes `-H` to vg, which leaves
    out a GBZ's haplotype walks. A caller must still expect them: vg keeps them for its other formats,
    and a GFA file is read as it is."""

    def __init__(self, path, vg=None, threads=4, drop_haplotypes=False):
        self.path, self.vg, self.threads, self.proc, self.fh = path, vg, threads, None, None
        self.drop_haplotypes = drop_haplotypes

    def __enter__(self):
        if self.path.endswith(VG_SUFFIXES):
            cmd = [find_vg(self.vg), 'convert', '-f', '-t', str(self.threads)]
            if self.drop_haplotypes:
                cmd.append('-H')
            cmd.append(self.path)
            self.proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
            self.fh = self.proc.stdout
        else:
            self.fh = open_text(self.path)
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
