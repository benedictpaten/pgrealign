"""A manifest written next to every stage's output: what ran, on what, with which tools, and when."""
import collections.abc
import datetime
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys

from . import __version__

# The arguments that make a tool print its version, by tool name, where they are not --version.
VERSION_ARGS = {'vg': ('version',)}


def tool_version(name, args=None, binary=None):
    """The binary a stage ran for the tool `name` and the version it reports, as {'path', 'version', 'error'}.

    `binary` is a path or a name looked up on PATH, and defaults to `name`. 'path' is its absolute path, or
    None if it cannot be found. 'version' is the first line that `<binary> <args>` prints, on stdout or
    else on stderr, where `args` defaults to VERSION_ARGS for the name and to --version otherwise. It is None
    when the binary cannot be run, exits non-zero or prints nothing, and 'error' then says which; otherwise
    'error' is None."""
    found = shutil.which(binary or name)
    if found is None:
        return {'path': None, 'version': None, 'error': 'not found: %s' % (binary or name)}
    path = os.path.abspath(found)
    if args is None:
        args = VERSION_ARGS.get(os.path.basename(name), ('--version',))
    try:
        p = subprocess.run([path, *args], stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, text=True, errors='replace', timeout=60)
    except (OSError, subprocess.TimeoutExpired) as e:
        return {'path': path, 'version': None, 'error': str(e)}
    if p.returncode != 0:
        lines = (p.stderr.strip() or p.stdout.strip()).splitlines()
        return {'path': path, 'version': None,
                'error': 'exit status %d' % p.returncode + (': ' + lines[0] if lines else '')}
    lines = (p.stdout.strip() or p.stderr.strip()).splitlines()
    if not lines:
        return {'path': path, 'version': None, 'error': 'no output'}
    return {'path': path, 'version': lines[0].strip(), 'error': None}


def _tools(tools):
    """{name: binary} from a mapping of the same (a binary of None means the name on PATH) or from a
    sequence of names or paths, each named by its basename."""
    if isinstance(tools, collections.abc.Mapping):
        return {name: binary or name for name, binary in tools.items()}
    if isinstance(tools, str):
        tools = [tools]
    out = {}
    for t in tools:
        name = os.path.basename(t)
        if out.get(name, t) != t:
            raise ValueError('two tools named %s: %s and %s' % (name, out[name], t))
        out[name] = t
    return out


def _timestamp(t):
    """`t`, a datetime or seconds since the epoch, as local time in ISO 8601 with its UTC offset."""
    if not isinstance(t, datetime.datetime):
        t = datetime.datetime.fromtimestamp(t)
    return t.astimezone().isoformat(timespec='seconds')


def file_info(path, checksum=False):
    """An input's path, size and modification time, and with `checksum` its MD5."""
    st = os.stat(path)
    out = {'path': os.path.abspath(path), 'bytes': st.st_size, 'mtime': _timestamp(st.st_mtime)}
    if checksum:
        h = hashlib.md5()
        with open(path, 'rb') as f:
            for block in iter(lambda: f.read(1 << 20), b''):
                h.update(block)
        out['md5'] = h.hexdigest()
    return out


def write_manifest(out_path, stage, inputs, params, tools=(), extra=None, started=None, checksums=False):
    """Write <out_path>.manifest.json: the stage, the command line, the inputs, the parameters, when the
    stage started and finished, and the versions of pgrealign, Python and the external tools.

    `tools` maps each tool's name to the binary the stage ran (None for the name on PATH), or is a sequence
    of names or paths, each named by its basename; the manifest records each as tool_version does.
    `started` is when the stage began, as a datetime or a time.time() value taken before it ran, and is
    null in the manifest if not given. 'finished' is when write_manifest is called. With `checksums`,
    each input's MD5 is recorded too, so that someone else can confirm they have the same inputs.
    `extra` is merged in
    last. Returns the manifest."""
    finished = datetime.datetime.now()
    m = {
        'stage': stage,
        'pgrealign': __version__,
        'python': platform.python_version(),
        'command': list(sys.argv),
        'started': None if started is None else _timestamp(started),
        'finished': _timestamp(finished),
        'inputs': {k: file_info(v, checksums) for k, v in inputs.items()},
        'params': params,
        'tools': {name: tool_version(name, binary=binary) for name, binary in _tools(tools or ()).items()},
    }
    if extra:
        m.update(extra)
    with open(out_path + '.manifest.json.tmp', 'w') as f:
        json.dump(m, f, indent=1, sort_keys=True)
        f.write('\n')
    os.replace(out_path + '.manifest.json.tmp', out_path + '.manifest.json')
    return m
