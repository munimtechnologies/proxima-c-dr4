#!/usr/bin/env python3
"""Validation of the gaiaorb fitter on the Gaia DR4 pre-release epoch
astrometry, and a dress rehearsal of the pre-registered Proxima c test on
the G ~ 9 pre-release star HD 183633 (a sensitivity forecast).

MV2 / Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name.
Rules: PREREGISTRATION.md (sections 3 and 4).

Usage: python validate_job.py --data-dir /data --out-dir /out
"""
import argparse
import csv
import glob
import json
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import io_gaia  # noqa: E402
from gaiaorb import epochs, prereg, validate, kepler  # noqa: E402

T0 = time.time()


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def write_csv(path, rows):
    if not rows:
        open(path, "w").close()
        return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, float) and o != o:
        return None
    return o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--skip-rehearsal", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.data_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)

    zp = io_gaia.download(io_gaia.PRERELEASE_URLS, os.path.join(args.data_dir, "dr4-prerelease.zip"), log=log)
    files = io_gaia.extract(zp, os.path.join(args.data_dir, "dr4-prerelease"))
    tables = [p for p in files if p.lower().endswith((".xml", ".vot", ".csv", ".ecsv", ".fits"))]
    tables.sort(key=lambda p: (0 if "epoch_astrometry" in os.path.basename(p).lower() else 1, p))
    log("files: %s" % [os.path.basename(p) for p in files])
    rows = io_gaia.read_table_rows(tables[0])
    log("rows: %d, columns: %s" % (len(rows), sorted(rows[0].keys()) if rows else []))
    by = epochs.rows_by_source(rows)

    summary = []
    prod, raw = {}, {}
    for sid, rs in sorted(by.items(), key=lambda kv: str(kv[0])):
        p_prod, c_prod = epochs.normal_points(rs, inflate=True)
        p_raw, _ = epochs.normal_points(rs, inflate=False)
        prod[sid], raw[sid] = p_prod, p_raw
        sig = sorted(p["sig"] for p in p_prod)
        summary.append(dict(source_id=sid, n_rows=len(rs), n_transits=c_prod["transits"],
                            n_points=c_prod["kept"], ccd_total=c_prod["ccd_total"], ccd_used=c_prod["ccd_used"],
                            too_few_ccd=c_prod["too_few_ccd"],
                            span_yr=(p_prod[-1]["t"] - p_prod[0]["t"]) if p_prod else None,
                            median_sigma_mas=sig[len(sig) // 2] if sig else None,
                            median_inflation=sorted(p["infl"] for p in p_prod)[len(p_prod) // 2] if p_prod else None))
        log("source %s: %s" % (sid, summary[-1]))
    write_csv(os.path.join(args.out_dir, "sources.csv"), summary)

    gates = []
    details = {}

    def run(name, fn):
        try:
            log("gate block " + name)
            g, d = fn()
            gates.extend(g)
            details[name] = d
        except Exception as e:  # a crash is a failed gate, never a pass
            gates.append(dict(gate=name, value=None, reference=None, tolerance=None, passed=False,
                              error="%s: %s" % (type(e).__name__, e)))
            details[name] = dict(error=traceback.format_exc())
            log("ERROR in %s: %s" % (name, e))

    def need(sid):
        if sid not in raw or len(raw[sid]) < 10:
            raise RuntimeError("source %s missing or too few points" % sid)
        return sid

    run("bh3_5p", lambda: (validate.check_bh3_5p(raw[need(validate.BH3)]), None))
    run("bh3_orbit", lambda: validate.check_bh3_orbit(raw[need(validate.BH3)]))
    run("hd114762_orbit", lambda: validate.check_dr3_orbit(raw[need(validate.HD114762)], validate.DR3_HD114762,
                                                           "hd114762", 1.0, 1.35, 2.25, e_tol=0.13))
    run("gaia4_orbit", lambda: validate.check_dr3_orbit(raw[need(validate.GAIA4)], validate.DR3_GAIA4,
                                                        "gaia4", 36.0, 0.16, 0.47))
    for k, q in enumerate(validate.QSOS):
        run("qso%d" % (k + 1), lambda q=q, k=k: validate.check_quasar(prod[need(q)], "qso%d" % (k + 1)))
    run("injection_recovery", lambda: validate.check_injection_recovery(prod[need(validate.HD183633)]))

    validated = bool(gates) and all(g["passed"] for g in gates)
    out = dict(project="proxima-c-dr4", part="validation on the Gaia DR4 pre-release", version="0.1.0",
               verdict="VALIDATED" if validated else "FAILED",
               n_gates=len(gates), n_passed=sum(1 for g in gates if g["passed"]),
               gates=gates, details=details, sources=summary, runtime_s=time.time() - T0)
    with open(os.path.join(args.out_dir, "validation.json"), "w") as f:
        json.dump(jsonable(out), f, indent=1)
    write_csv(os.path.join(args.out_dir, "gates.csv"),
              [dict(gate=g["gate"], value=g["value"], reference=g["reference"], tolerance=g["tolerance"],
                    passed=g["passed"], error=g.get("error", "")) for g in gates])
    log("VALIDATION %s (%d/%d gates)" % (out["verdict"], out["n_passed"], out["n_gates"]))

    if args.skip_rehearsal:
        return
    # Dress rehearsal / forecast: the exact pre-registered Proxima c test run on
    # HD 183633 (G ~ 9, like Proxima). Its own perspective acceleration is
    # negligible (plx ~ 1 mas), so mu_r = 0; injections use Proxima's parallax
    # and mass, so the result is the detectability of Proxima c at G ~ 9.
    star = dict(prereg.PROXIMA, rv_kms=0.0, name="HD 183633 (rehearsal with Proxima c injections)")
    try:
        res = prereg.run_test(prod[need(validate.HD183633)], star, log=log)
        res["verdict_if_this_were_proxima"] = prereg.verdict(res, validated, star)
    except Exception as e:
        res = dict(status="error", error=traceback.format_exc())
        log("ERROR in rehearsal: %s" % e)
    res["note"] = ("Rehearsal on HD 183633 epochs, not Proxima. Proxima's DR4 transit count and "
                   "noise will differ; this is a forecast of sensitivity at G ~ 9.")
    with open(os.path.join(args.out_dir, "rehearsal_hd183633.json"), "w") as f:
        json.dump(jsonable(res), f, indent=1)
    if res.get("status") == "ok":
        write_csv(os.path.join(args.out_dir, "rehearsal_injections.csv"), res["injections"])
        write_csv(os.path.join(args.out_dir, "rehearsal_periodogram.csv"), res["periodogram"])
        log("rehearsal: S_obs=%.2f thr=%.2f limits=%s" % (res["S_obs"], res["threshold"], json.dumps(res["limits"])))
    log("done")


if __name__ == "__main__":
    main()
