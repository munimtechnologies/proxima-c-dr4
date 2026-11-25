"""DR4-day launcher for the pre-registered Proxima c test.

MV2, Munim, Inc. (mv2@munimtech.com). Apache-2.0. Made with AI agents working in MV2's name.

Downloads this repository at a fixed tag (v0.2.1: the frozen v0.1.0 analysis code with the v0.2.0 DataLink
I/O of 2026-11-22), installs its requirements, runs its unit tests and then runs proxima_job.py with the
remaining arguments, for example:

    python dr4_bootstrap.py --tag v0.2.1 -- --data-dir /data --out-dir /out --validated

Standard library only, so the job's workdir needs just this one file.
"""
import argparse
import io
import os
import subprocess
import sys
import tempfile
import urllib.request
import zipfile

REPO = "munimtechnologies/proxima-c-dr4"


def archive_url(tag):
    return "https://github.com/%s/archive/refs/tags/%s.zip" % (REPO, tag)


def safe_names(names):
    """Raise if any archive member would land outside the extraction folder."""
    for n in names:
        parts = n.replace("\\", "/").split("/")
        if n.startswith("/") or n.startswith("\\") or ".." in parts or (parts and ":" in parts[0]):
            raise ValueError("unsafe path in archive: %r" % n)
    return True


def extract(blob, dest):
    """Extract a GitHub archive zip into dest and return the path of its single top folder."""
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        names = z.namelist()
        if not names:
            raise ValueError("empty archive")
        safe_names(names)
        tops = {n.replace("\\", "/").split("/")[0] for n in names}
        if len(tops) != 1:
            raise ValueError("archive has %d top folders, expected 1" % len(tops))
        z.extractall(dest)
    return os.path.join(dest, tops.pop())


def split_args(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="v0.2.1")
    ap.add_argument("--skip-tests", action="store_true")
    ap.add_argument("--skip-install", action="store_true")
    if "--" in argv:
        i = argv.index("--")
        own, rest = argv[:i], argv[i + 1:]
    else:
        own, rest = argv, []
    return ap.parse_args(own), rest


def main(argv=None):
    args, rest = split_args(sys.argv[1:] if argv is None else argv)
    url = archive_url(args.tag)
    print("[bootstrap] downloading %s" % url, flush=True)
    blob = urllib.request.urlopen(url, timeout=180).read()
    root = extract(blob, tempfile.mkdtemp(prefix="proxima-"))
    print("[bootstrap] extracted to %s" % root, flush=True)
    if not args.skip_install and os.path.exists(os.path.join(root, "requirements.txt")):
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", "requirements.txt"], cwd=root)
    if not args.skip_tests:
        subprocess.check_call([sys.executable, "-m", "unittest", "-v", "test_gaiaorb"], cwd=root)
    return subprocess.call([sys.executable, "proxima_job.py"] + rest, cwd=root)


if __name__ == "__main__":
    sys.exit(main())
