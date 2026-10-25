#!/usr/bin/env python3
"""Reconstruct the published atlas and copy its small companion files."""
import base64
import gzip
import hashlib
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '_site'
BASE_COMMIT = '189a315dd1d2d85ea6175807e3fbe4ecd1b1b22c'
BASE_DIGEST = 'f2cc92f2bbe989dde49bbfec73b81f57ecc9273a55ef834bc5d06c1959e096aa'
PATCHES = ['compact-spreads.patch', 'storytelling.patch', 'galton-reader.patch',
           'huygens-axis.patch', 'testing-story.patch']


def apply_patch(content):
    subprocess.run(['patch', '--batch', '--forward', str(OUT / 'index.html')],
                   input=content, check=True, stdout=subprocess.PIPE)


def main():
    OUT.mkdir(exist_ok=True)
    original = subprocess.check_output(['git', 'show', f'{BASE_COMMIT}:index.html'], cwd=ROOT)
    (OUT / 'index.html').write_bytes(original)
    encoded = b''.join(p.read_bytes() for p in sorted((ROOT / 'updates').glob('archive-*.b64')))
    apply_patch(gzip.decompress(base64.b64decode(encoded)))
    for name in PATCHES:
        apply_patch((ROOT / 'updates' / name).read_bytes())
    actual = hashlib.sha256((OUT / 'index.html').read_bytes()).hexdigest()
    if actual != BASE_DIGEST:
        raise SystemExit(f'Historical source verification failed: {actual}')
    apply_patch((ROOT / 'updates' / 'statistics-editorial.patch').read_bytes())
    for name in ['.nojekyll', 'NOTICE.txt', 'atlas-editorial.js',
                 'atlas-statistics.js', 'atlas-statistics.css']:
        shutil.copyfile(ROOT / name, OUT / name)
    print(f'Built verified atlas in {OUT}')


if __name__ == '__main__':
    main()
