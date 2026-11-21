"""Downloads and file readers for Gaia DR4 epoch astrometry (I/O only).

Part of proxima-c-dr4 (MV2, Munim, Inc.). Apache-2.0. Made with AI agents.
numpy/astropy are imported lazily so that the tested core stays stdlib-only.
"""
import glob
import io
import os
import time
import urllib.parse
import urllib.request
import zipfile

from gaiaorb import epochs

PRERELEASE_URLS = (
    "https://anonftp.cosmos.esa.int/pub/GAIA_PUBLIC_DATA/Gaia_DR4/dr4-prerelease/"
    "gaia-dr4-prerelease-epoch-astrometry_2026-06-26.zip",
    "http://anonftp.cosmos.esa.int/pub/GAIA_PUBLIC_DATA/Gaia_DR4/dr4-prerelease/"
    "gaia-dr4-prerelease-epoch-astrometry_2026-06-26.zip",
)
UA = {"User-Agent": "proxima-c-dr4/0.1 (open research code)"}


def fetch(url, timeout=600):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def download(urls, dest, tries=4, log=print):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    last = None
    for attempt in range(tries):
        for url in urls:
            try:
                blob = fetch(url)
                tmp = dest + ".part"
                with open(tmp, "wb") as f:
                    f.write(blob)
                os.replace(tmp, dest)
                log("downloaded %s (%d bytes)" % (url, len(blob)))
                return dest
            except Exception as e:  # try the next URL, then back off
                last = e
                log("download failed %s: %s" % (url, e))
        time.sleep(15 * (attempt + 1))
    raise RuntimeError("download failed: %s" % last)


def extract(zip_path, dest_dir):
    os.makedirs(dest_dir, exist_ok=True)
    with zipfile.ZipFile(zip_path) as z:
        for name in z.namelist():
            if name.endswith("/"):
                continue
            target = os.path.join(dest_dir, os.path.basename(name))
            with open(target, "wb") as f:
                f.write(z.read(name))
    return sorted(glob.glob(os.path.join(dest_dir, "*")))


def _py(v):
    import numpy as np
    if v is np.ma.masked:
        return None
    if isinstance(v, (bytes, np.bytes_)):
        return v.decode()
    if isinstance(v, str):
        return epochs.parse_cell(v) if v[:1] in "([" else v
    if np.ndim(v) > 0:
        a = np.ma.asarray(v)
        mask = np.ma.getmaskarray(a).tolist()
        vals = a.data.tolist()
        return [None if m else x for x, m in zip(vals, mask)]
    if isinstance(v, np.generic):
        return v.item()
    return v


def read_table_rows(path):
    ext = path.lower().rsplit(".", 1)[-1]
    if ext == "csv":
        return epochs.read_csv_rows(path)
    if ext in ("xml", "vot", "votable"):
        from astropy.io.votable import parse
        tab = parse(path).get_first_table().to_table(use_names_over_ids=True)
    else:
        from astropy.table import Table
        tab = Table.read(path)
    names = list(tab.colnames)
    rows = []
    for r in tab:
        rows.append({n: _py(r[n]) for n in names})
    return rows


def tap_csv(query, url="https://gea.esac.esa.int/tap-server/tap/sync"):
    q = urllib.parse.urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query))
    text = fetch(url + "?" + q).decode("utf-8", "replace")
    import csv
    return list(csv.DictReader(io.StringIO(text)))


def datalink_urls(source_id, release="Gaia DR4"):
    base = "https://gea.esac.esa.int/data-server/data?"
    out = []
    for ident in ("%s %d" % (release, source_id), str(source_id)):
        for fmt in ("CSV", "votable"):
            out.append(base + urllib.parse.urlencode(dict(
                RETRIEVAL_TYPE="EPOCH_ASTROMETRY", ID=ident, FORMAT=fmt, RELEASE=release,
                DATA_STRUCTURE="INDIVIDUAL")))
    return out


def save_epochs(blob, dest_dir, stem):
    """Store a DataLink response (zip, csv or votable); return table paths."""
    os.makedirs(dest_dir, exist_ok=True)
    if blob[:2] == b"PK":
        zp = os.path.join(dest_dir, stem + ".zip")
        with open(zp, "wb") as f:
            f.write(blob)
        return [p for p in extract(zp, os.path.join(dest_dir, stem))
                if p.lower().rsplit(".", 1)[-1] in ("csv", "xml", "vot", "ecsv", "fits")]
    head = blob[:200].lstrip().lower()
    ext = "xml" if head.startswith(b"<?xml") or head.startswith(b"<votable") else "csv"
    p = os.path.join(dest_dir, stem + "." + ext)
    with open(p, "wb") as f:
        f.write(blob)
    return [p]
