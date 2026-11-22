"""Downloads and file readers for Gaia DR4 epoch astrometry (I/O only).

Shared by proxima-c-dr4, radio-planets-dr4, pceb-planets-dr4 and gj3512-giants-dr4
(MV2, Munim, Inc.). Apache-2.0. Made with AI agents working in MV2's name.
numpy/astropy are imported lazily so that the tested core stays stdlib-only.

v0.2.0 (2026-11-22, before Gaia DR4; I/O only, no analysis change):
- Gaia Archive v3.15 (2026-10-14) renamed the DataLink data structures: INDIVIDUAL is now
  DATAMODEL_STANDARD and RAW is now DATAMODEL_GAIA. datalink_urls() now asks for
  DATAMODEL_STANDARD first, then the old INDIVIDUAL, then no DATA_STRUCTURE at all.
- save_epochs() now keeps only files that parse into epoch-astrometry rows (a centroid_pos_al
  column and at least one row). An archive error page, an empty table or a file of another
  product is rejected, so the calling job moves on to the next URL instead of crashing later.
"""
import glob
import io
import os
import time
import urllib.parse
import urllib.request
import zipfile

from gaiaorb import epochs

VERSION = "0.2.0"
PRERELEASE_URLS = (
    "https://anonftp.cosmos.esa.int/pub/GAIA_PUBLIC_DATA/Gaia_DR4/dr4-prerelease/"
    "gaia-dr4-prerelease-epoch-astrometry_2026-06-26.zip",
    "http://anonftp.cosmos.esa.int/pub/GAIA_PUBLIC_DATA/Gaia_DR4/dr4-prerelease/"
    "gaia-dr4-prerelease-epoch-astrometry_2026-06-26.zip",
)
UA = {"User-Agent": "gaiaorb-dr4/0.2 (open research code)"}
DATALINK_BASE = "https://gea.esac.esa.int/data-server/data?"
# Newest name first (Gaia Archive release 3.15, 2026-10-14); None = parameter left out.
DATA_STRUCTURES = ("DATAMODEL_STANDARD", "INDIVIDUAL", None)
TABLE_EXTS = ("csv", "xml", "vot", "votable", "ecsv", "fits")
REQUIRED_COLUMN = "centroid_pos_al"


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
    """Candidate DataLink URLs for one source's epoch astrometry, best first."""
    out = []
    for structure in DATA_STRUCTURES:
        for ident in ("%s %d" % (release, source_id), str(source_id)):
            for fmt in ("CSV", "votable"):
                params = dict(RETRIEVAL_TYPE="EPOCH_ASTROMETRY", ID=ident, FORMAT=fmt, RELEASE=release)
                if structure is not None:
                    params["DATA_STRUCTURE"] = structure
                out.append(DATALINK_BASE + urllib.parse.urlencode(params))
    return out


def looks_like_epochs(rows):
    """True if rows are epoch astrometry: at least one row with a centroid_pos_al column."""
    return bool(rows) and any(REQUIRED_COLUMN in r for r in rows)


def usable_tables(paths, log=print):
    """Keep only the paths that parse into epoch-astrometry rows."""
    keep = []
    for p in paths:
        try:
            rows = read_table_rows(p)
        except Exception as e:  # an error page or a product we cannot read
            log("rejected %s: unreadable (%s)" % (p, str(e)[:200]))
            continue
        if looks_like_epochs(rows):
            keep.append(p)
        else:
            cols = sorted(rows[0].keys())[:12] if rows else []
            log("rejected %s: %d rows, columns %s" % (p, len(rows), cols))
    return keep


def save_epochs(blob, dest_dir, stem, log=print):
    """Store a DataLink response (zip, csv or votable); return the paths of usable epoch tables
    (an empty list if the response is an error page or holds no epoch astrometry)."""
    os.makedirs(dest_dir, exist_ok=True)
    if blob[:2] == b"PK":
        zp = os.path.join(dest_dir, stem + ".zip")
        with open(zp, "wb") as f:
            f.write(blob)
        try:
            paths = extract(zp, os.path.join(dest_dir, stem))
        except zipfile.BadZipFile as e:
            log("rejected %s: bad zip (%s)" % (zp, e))
            return []
        paths = [p for p in paths if p.lower().rsplit(".", 1)[-1] in TABLE_EXTS]
        return usable_tables(paths, log)
    head = blob[:200].lstrip().lower()
    ext = "xml" if head.startswith(b"<?xml") or head.startswith(b"<votable") else "csv"
    if head.startswith(b"<!doctype html") or head.startswith(b"<html"):
        log("rejected %s response: HTML page" % stem)
        return []
    p = os.path.join(dest_dir, stem + "." + ext)
    with open(p, "wb") as f:
        f.write(blob)
    return usable_tables([p], log)
