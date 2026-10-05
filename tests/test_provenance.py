import datetime
import json
import os
import shutil
import tempfile
import time
import unittest
from unittest import mock

from pgrealign import provenance

# A stand-in for vg: it rejects --version, as vg does, and answers `version`.
FAKE_VG = '''#!/bin/sh
if [ "$1" = version ]; then echo "vg version v0.0.0-test"; echo "Compiled with nothing"; exit 0; fi
echo "error[vg] command $1 not found" >&2
exit 1
'''


class TestProvenance(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d)

    def script(self, name, text, mode=0o755):
        path = os.path.join(self.d, name)
        with open(path, 'w') as f:
            f.write(text)
        os.chmod(path, mode)
        return path

    def test_vg_is_asked_with_version(self):
        vg = self.script('vg', FAKE_VG)
        self.assertEqual(provenance.tool_version('vg', binary=vg),
                         {'path': vg, 'version': 'vg version v0.0.0-test', 'error': None})
        # under another file name it is still vg, when named so
        other = self.script('vg-test', FAKE_VG)
        self.assertEqual(provenance.tool_version('vg', binary=other)['version'], 'vg version v0.0.0-test')

    def test_other_tools_are_asked_with_double_dash_version(self):
        tool = self.script('tool', '#!/bin/sh\n[ "$1" = --version ] && echo "tool 1.2" >&2 && exit 0\nexit 3\n')
        self.assertEqual(provenance.tool_version('tool', binary=tool)['version'], 'tool 1.2')
        self.assertEqual(provenance.tool_version('tool', args=('-v',), binary=tool),
                         {'path': tool, 'version': None, 'error': 'exit status 3'})

    def test_a_failure_is_never_a_version(self):
        vg = self.script('vg', FAKE_VG)
        r = provenance.tool_version('vg', args=('--version',), binary=vg)
        self.assertEqual(r, {'path': vg, 'version': None,
                             'error': 'exit status 1: error[vg] command --version not found'})
        silent = self.script('silent', '#!/bin/sh\nexit 0\n')
        self.assertEqual(provenance.tool_version('silent', binary=silent)['error'], 'no output')
        # executable but not a program: the same keys, with the path, not a bare path
        junk = self.script('junk', '\x00\x01\x02')
        r = provenance.tool_version('junk', binary=junk)
        self.assertEqual((r['path'], r['version']), (junk, None))
        self.assertTrue(r['error'])

    def test_not_found(self):
        missing = os.path.join(self.d, 'missing')
        self.assertEqual(provenance.tool_version('missing', binary=missing),
                         {'path': None, 'version': None, 'error': 'not found: ' + missing})

    def test_name_on_path_is_resolved_to_an_absolute_path(self):
        vg = self.script('vg', FAKE_VG)
        with mock.patch.dict(os.environ, {'PATH': self.d}):
            self.assertEqual(provenance.tool_version('vg')['path'], vg)

    def test_tools_by_name_path_or_mapping(self):
        vg = self.script('vg', FAKE_VG)
        other = self.script('vg-test', FAKE_VG)
        self.assertEqual(provenance._tools(['vg', '/opt/x/abpoa']), {'vg': 'vg', 'abpoa': '/opt/x/abpoa'})
        self.assertEqual(provenance._tools({'vg': other, 'abpoa': None}), {'vg': other, 'abpoa': 'abpoa'})
        self.assertEqual(provenance._tools('vg'), {'vg': 'vg'})
        with self.assertRaises(ValueError):
            provenance._tools([vg, '/elsewhere/vg'])

    def test_manifest(self):
        vg = self.script('vg', FAKE_VG)
        other = self.script('vg-test', FAKE_VG)
        inp = os.path.join(self.d, 'in.txt')
        with open(inp, 'w') as f:
            f.write('abc')
        out = os.path.join(self.d, 'out.tsv')
        t0 = time.time() - 5

        # the positional form callers already use, with the start time as a time.time() value
        m = provenance.write_manifest(out, 'stage', {'x': inp}, {'p': 1}, [vg], {'summary': 2}, started=t0)
        with open(out + '.manifest.json') as f:
            self.assertEqual(json.load(f), m)
        self.assertEqual((m['stage'], m['params'], m['summary']), ('stage', {'p': 1}, 2))
        self.assertEqual(m['inputs']['x']['bytes'], 3)
        self.assertEqual(m['tools'], {'vg': {'path': vg, 'version': 'vg version v0.0.0-test', 'error': None}})
        started = datetime.datetime.fromisoformat(m['started'])
        finished = datetime.datetime.fromisoformat(m['finished'])
        self.assertIsNotNone(started.utcoffset())
        self.assertLessEqual(abs(started.timestamp() - t0), 1)
        self.assertGreaterEqual((finished - started).total_seconds(), 4)

        # a mapping names the binary; a datetime is a start time too
        when = datetime.datetime(2026, 1, 2, 3, 4, 5)
        m = provenance.write_manifest(out, 'stage', {}, {}, tools={'vg': other}, started=when)
        self.assertEqual(m['tools']['vg']['path'], other)
        self.assertEqual(datetime.datetime.fromisoformat(m['started']).replace(tzinfo=None), when)

        # without a start time it is null, never the time of writing
        m = provenance.write_manifest(out, 'stage', {}, {})
        self.assertIsNone(m['started'])
        self.assertEqual(m['tools'], {})


@unittest.skipUnless(shutil.which('vg'), 'vg not on PATH')
class TestRealVg(unittest.TestCase):
    def test_vg_reports_its_version(self):
        r = provenance.tool_version('vg')
        self.assertIsNone(r['error'])
        self.assertTrue(r['version'].startswith('vg version '), r['version'])
        self.assertEqual(r['path'], os.path.abspath(shutil.which('vg')))
        self.assertEqual(provenance.tool_version('vg', binary=r['path']), r)


@unittest.skipUnless(shutil.which('abpoa'), 'abpoa not on PATH')
class TestRealAbpoa(unittest.TestCase):
    def test_abpoa_reports_its_version(self):
        r = provenance.tool_version('abpoa')
        self.assertIsNone(r['error'])
        self.assertTrue(r['version'])


if __name__ == '__main__':
    unittest.main()
