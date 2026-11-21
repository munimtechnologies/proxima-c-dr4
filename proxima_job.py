#!/usr/bin/env python3
"""The pre-registered Proxima c test on Gaia DR4 epoch astrometry.
To be run on or after the Gaia DR4 release (2026-12-02).

MV2 / Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name.
Rules: PREREGISTRATION.md. Only the I/O (source lookup, download URLs and
file format) may be adapted after the release; the analysis may not.

Usage:
  python proxima_job.py --data-dir /data --out-dir /out [--validated]
        [--source-id N] [--epochs-file PATH] [--epochs-url URL]
"""
import argparse
import json
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import io_gaia  # noqa: E402
from gaiaorb import epochs, prereg  # noqa: E402
from validate_job import jsonable, write_csv  # noqa: E402

T0 = time.time()
# Proxima at J2017.5 (DR3 position at J2016.0 moved by its proper motion)
RA_2017_5, DEC_2017_5 = 217.38889, -62.67575


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def resolve_source(log=log):
    for table in ("gaiadr4.gaia_source", "gaiadr4.gaia_source_lite"):
        q = ("SELECT source_id, ra, dec, parallax, phot_g_mean_mag FROM %s WHERE "
             "1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', %.5f, %.5f, 0.02)) "
             "AND parallax > 700" % (table, RA_2017_5, DEC_2017_5))
        try:
            rows = io_gaia.tap_csv(q)
            log("%s: %s" % (table, rows))
            if rows:
                return int(rows[0]["source_id"]), rows[0]
        except Exception as e:
            log("TAP %s failed: %s" % (table, e))
    return prereg.PROXIMA["dr3_source_id"], dict(note="DR4 lookup failed; using the DR3 source_id")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--validated", action="store_true",
                    help="set only if the validation job's verdict was VALIDATED")
    ap.add_argument("--source-id", type=int, default=None)
    ap.add_argument("--epochs-file", default=None)
    ap.add_argument("--epochs-url", default=None)
    args = ap.parse_args()
    os.makedirs(args.data_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)

    if args.source_id is not None:
        sid, src = args.source_id, dict(note="given on the command line")
    else:
        sid, src = resolve_source()
    log("source_id %s" % sid)

    paths = []
    if args.epochs_file:
        paths = [args.epochs_file]
    else:
        urls = [args.epochs_url] if args.epochs_url else io_gaia.datalink_urls(sid)
        for url in urls:
            try:
                blob = io_gaia.fetch(url)
                if len(blob) < 200:
                    log("short response from %s: %r" % (url, blob[:200]))
                    continue
                paths = io_gaia.save_epochs(blob, args.data_dir, "proxima-epochs")
                log("epochs from %s -> %s" % (url, paths))
                if paths:
                    break
            except Exception as e:
                log("fetch failed %s: %s" % (url, e))
    if not paths:
        raise SystemExit("no epoch astrometry found for source %s" % sid)

    rows = []
    for p in paths:
        rows.extend(io_gaia.read_table_rows(p))
    by = epochs.rows_by_source(rows)
    mine = by.get(sid) or (list(by.values())[0] if len(by) == 1 else [])
    points, counts = epochs.normal_points(mine, inflate=True)
    log("rows %d -> points %d (%s)" % (len(mine), len(points), counts))

    try:
        res = prereg.run_test(points, prereg.PROXIMA, log=log)
    except Exception:
        res = dict(status="error", error=traceback.format_exc())
    res["verdict"] = prereg.verdict(res, args.validated)
    res.update(project="proxima-c-dr4", version="0.1.0", source_id=sid, source=src,
               validated_flag=bool(args.validated), counts=counts, runtime_s=time.time() - T0)
    with open(os.path.join(args.out_dir, "results.json"), "w") as f:
        json.dump(jsonable(res), f, indent=1)
    if res.get("status") == "ok":
        write_csv(os.path.join(args.out_dir, "injections.csv"), res["injections"])
        write_csv(os.path.join(args.out_dir, "periodogram.csv"), res["periodogram"])
        write_csv(os.path.join(args.out_dir, "points.csv"), points)
    log("VERDICT %s" % res["verdict"])


if __name__ == "__main__":
    main()
