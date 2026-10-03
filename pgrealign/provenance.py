"""A manifest written next to every stage's output: what ran, on what, with which tools."""
import datetime
import json
import os
import platform
import shutil
import subprocess
import sys

from . import __version__


def tool_version(name, args=("--version",)):
    """The first line a tool prints for its version, or None if it is not on PATH."""
    path = shutil.which(name)
    if path is None:
        return None
    try:
        p = subprocess.run([path, *args], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return path
    text = (p.stdout or p.stderr).strip().splitlines()
    return '%s (%s)' % (text[0] if text else '?', path)


def file_info(path):
    st = os.stat(path)
    return {'path': os.path.abspath(path), 'bytes': st.st_size,
            'mtime': datetime.datetime.fromtimestamp(st.st_mtime).isoformat(timespec='seconds')}


def write_manifest(out_path, stage, inputs, params, tools=(), extra=None):
    """Write <out_path>.manifest.json: the stage, the command line, the inputs, the parameters and the
    versions of pgrealign, Python and the named external tools."""
    m = {
        'stage': stage,
        'pgrealign': __version__,
        'python': platform.python_version(),
        'command': sys.argv,
        'started': datetime.datetime.now().isoformat(timespec='seconds'),
        'inputs': {k: file_info(v) for k, v in inputs.items()},
        'params': params,
        'tools': {t: tool_version(t) for t in tools},
    }
    if extra:
        m.update(extra)
    with open(out_path + '.manifest.json', 'w') as f:
        json.dump(m, f, indent=1, sort_keys=True)
        f.write('\n')
    return m
