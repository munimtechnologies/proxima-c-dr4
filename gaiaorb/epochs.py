"""Gaia DR4 epoch astrometry rows -> per-transit normal points.

Part of proxima-c-dr4 (MV2, Munim, Inc.). Apache-2.0. Made with AI agents.

Format (Gaia Archive / DataLink, as read by ESA's gaiasupdate package):
one row per transit; CCD-level arrays (about 10 entries: SM + 9 AF) for
obs_time_tcb (int ns since 2010-01-01T00:00 TCB), centroid_pos_al (mas),
centroid_pos_error_al (mas), scan_pos_angle (deg), used_by_agis_al (bool);
per-transit scalars obs_time_bary_corr (ns), parallax_factor_al.
Exploded tables (one row per CCD) are accepted too.

Time (gaiasupdate.set_relative_time): t [Julian yr] = 2010.0
  + (obs_time_tcb + obs_time_bary_corr) ns -> yr - 2017.5.
Light-time term (gaiasupdate perspective acceleration):
  fr [au] = -obs_time_bary_corr / 1e9 / 499.004783836156.
"""
import csv
import math

NS_PER_JYEAR = 365.25 * 86400.0 * 1e9
TCB_REF_JYEAR = 2010.0      # 2010-01-01T00:00:00 TCB as a Julian year
DR4_REF_JYEAR = 2017.5
AU_LIGHT_TIME_S = 499.004783836156

ARRAY_FIELDS = ("obs_time_tcb", "centroid_pos_al", "centroid_pos_error_al", "scan_pos_angle",
                "used_by_agis_al")


def rel_time_years(tcb_ns, bary_ns):
    return TCB_REF_JYEAR + float(tcb_ns) / NS_PER_JYEAR + float(bary_ns) / NS_PER_JYEAR - DR4_REF_JYEAR


def _scalar(tok):
    t = tok.strip().strip('"').strip("'")
    lo = t.lower()
    if lo in ("true", "t"):
        return True
    if lo in ("false", "f"):
        return False
    if lo in ("", "nan", "null", "none", "--"):
        return None
    try:
        return int(t)
    except ValueError:
        pass
    try:
        return float(t)
    except ValueError:
        return t


def parse_cell(s):
    """'(1.0, NaN, 3)' or '[true, false]' -> list; '0.68' -> 0.68; '' -> None."""
    if s is None:
        return None
    t = s.strip()
    if len(t) >= 2 and t[0] in "([{" and t[-1] in ")]}":
        inner = t[1:-1].strip()
        if not inner:
            return []
        return [_scalar(x) for x in inner.split(",")]
    return _scalar(t)


def read_csv_rows(path):
    rows = []
    with open(path, newline="") as fh:
        for rec in csv.DictReader(fh):
            rows.append({k.strip(): parse_cell(v) for k, v in rec.items() if k is not None})
    return rows


def _num(x):
    if x is None or isinstance(x, bool):
        return None
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _is_list(x):
    return isinstance(x, (list, tuple))


def ccd_groups(rows):
    """Group CCD-level measurements by transit. Returns a list of lists of
    dicts with keys tcb, bary, w, sw, psi, f, used."""
    groups = {}
    order = []
    for k, r in enumerate(rows):
        t = r.get("obs_time_tcb")
        if _is_list(t):
            n = len(t)
            key = ("row", k)
            items = []
            for j in range(n):
                def get(name, j=j):
                    v = r.get(name)
                    if _is_list(v):
                        return v[j] if j < len(v) else None
                    return v
                items.append(dict(tcb=get("obs_time_tcb"), bary=get("obs_time_bary_corr"),
                                  w=get("centroid_pos_al"), sw=get("centroid_pos_error_al"),
                                  psi=get("scan_pos_angle"), f=get("parallax_factor_al"),
                                  used=get("used_by_agis_al")))
        else:
            tid = r.get("transit_id")
            key = ("tid", tid) if tid is not None else ("bf", r.get("obs_time_bary_corr"),
                                                       r.get("parallax_factor_al"))
            items = [dict(tcb=t, bary=r.get("obs_time_bary_corr"), w=r.get("centroid_pos_al"),
                          sw=r.get("centroid_pos_error_al"), psi=r.get("scan_pos_angle"),
                          f=r.get("parallax_factor_al"), used=r.get("used_by_agis_al"))]
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].extend(items)
    return [groups[k] for k in order]


def normal_points(rows, min_ccd=3, inflate=True):
    """Collapse each transit's usable AF CCDs (used_by_agis_al True, finite
    values, time > 0) to one weighted-mean AL point. sigma is the error of
    the weighted mean, inflated by sqrt(chi2_within/(n-1)) when > 1 and
    `inflate` is set. Returns (points sorted by time, counts)."""
    counts = dict(transits=0, kept=0, too_few_ccd=0, ccd_total=0, ccd_used=0)
    pts = []
    for g in ccd_groups(rows):
        counts["transits"] += 1
        good = []
        for m in g:
            counts["ccd_total"] += 1
            if m["used"] is not True:
                continue
            tcb, bary, w, sw = _num(m["tcb"]), _num(m["bary"]), _num(m["w"]), _num(m["sw"])
            psi, f = _num(m["psi"]), _num(m["f"])
            if None in (tcb, bary, w, sw, psi, f) or tcb <= 0.0 or sw <= 0.0:
                continue
            good.append((rel_time_years(m["tcb"], bary), w, sw, math.radians(psi), f, bary))
        counts["ccd_used"] += len(good)
        if len(good) < min_ccd:
            counts["too_few_ccd"] += 1
            continue
        W = sum(1.0 / (g_[2] ** 2) for g_ in good)
        wm = sum(g_[1] / g_[2] ** 2 for g_ in good) / W
        tm = sum(g_[0] / g_[2] ** 2 for g_ in good) / W
        sm = sum(math.sin(g_[3]) / g_[2] ** 2 for g_ in good)
        cm = sum(math.cos(g_[3]) / g_[2] ** 2 for g_ in good)
        psi_m = math.atan2(sm, cm)
        fm = sum(g_[4] / g_[2] ** 2 for g_ in good) / W
        bm = sum(g_[5] / g_[2] ** 2 for g_ in good) / W
        chi2w = sum((g_[1] - wm) ** 2 / g_[2] ** 2 for g_ in good)
        infl = max(1.0, chi2w / (len(good) - 1)) if inflate else 1.0
        sig = math.sqrt(infl / W)
        pts.append(dict(t=tm, w=wm, sig=sig, s=math.sin(psi_m), c=math.cos(psi_m), f=fm,
                        fr=-bm / 1e9 / AU_LIGHT_TIME_S, n=len(good), infl=infl))
        counts["kept"] += 1
    pts.sort(key=lambda p: p["t"])
    return pts, counts


def rows_by_source(rows):
    out = {}
    for r in rows:
        sid = r.get("source_id")
        try:
            sid = int(sid)
        except (TypeError, ValueError):
            pass
        out.setdefault(sid, []).append(r)
    return out
