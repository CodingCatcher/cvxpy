#!/usr/bin/env python3
"""按笔记 14 的"打孔为主"框架，给定坝顶高程，把坝线按交付物 1 的 12 段汇总：
每段的坝高、CFRD / RCC 的上下游坝趾距离、能否闭合、下游坝趾出宗地的比例、平均断面法方量。
用法：python3 scripts/balanced_reach_table.py A_parcel_ring_3700 3850 [3900 ...]
输出：outputs/balanced_reaches_<line>_<crest>.csv 并打印。"""
import csv, json, os, sys
import numpy as np
from shapely.geometry import LineString, Polygon, Point
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from terrain_model import (Terrain, SECTIONS, FT, DATA, OUT, cross_section, dam_section,
                           inward_normals, longitudinal_profile)  # noqa

LINE = sys.argv[1] if len(sys.argv) > 1 else "A_parcel_ring_3700"
CRESTS = [float(x) for x in sys.argv[2:]] or [3850.0]
TYPES = ("cfrd", "rcc")
t = Terrain(os.path.join(DATA, "dem_3dep_1m_utm10.tif"), os.path.join(DATA, "dem_3dep_1m_utm10.json"))
cand = json.load(open(os.path.join(DATA, "dam_line_candidates.json")))[LINE]
line = LineString(cand["coords"]); ring = Polygon(cand["coords"])
from terrain_model import load_parcel  # noqa
parcel = load_parcel()
s, pts, zg = longitudinal_profile(t, line, 10.0); nrm = inward_normals(pts, ring, None)
secs = [cross_section(t, pts[i], nrm[i]) for i in range(len(pts))]
reaches = list(csv.DictReader(open(os.path.join(OUT, "D1_reaches.csv")))) if LINE == "A_parcel_ring_3700" else None

for crest_ft in CRESTS:
    crest = crest_ft * FT; nwl = crest - 20 * FT
    per = {ty: [dam_section(*secs[i], crest, SECTIONS[ty], nwl) for i in range(len(pts))] for ty in TYPES}
    H = crest - zg
    if reaches is None:  # 无预定义分段：按 400 m 均分
        reaches = [{"reach": k + 1, "sta_from": k * 400.0, "sta_to": min((k + 1) * 400.0, s[-1])} for k in range(int(np.ceil(s[-1] / 400.0)))]
    rows = []
    for r in reaches:
        a, b = float(r["sta_from"]), float(r["sta_to"]); idx = np.where((s >= a - 1e-6) & (s <= b + 1e-6))[0]
        if len(idx) == 0: continue
        rec = {"reach": r["reach"], "sta_from": a, "sta_to": b, "length_m": float(s[idx[-1]] - s[idx[0]]) + 10,
               "H_min": float(np.nanmin(H[idx])), "H_max": float(np.nanmax(H[idx]))}
        for ty in TYPES:
            dd = [per[ty][i] for i in idx]
            areas = np.array([d["area"] if d["area"] is not None else np.nan for d in dd])
            tu = np.array([d["toe_us"] if d["toe_us"] is not None else np.nan for d in dd])
            td = np.array([d["toe_ds"] if d["toe_ds"] is not None else np.nan for d in dd])
            closes = np.isfinite(areas)
            # 下游坝趾是否在宗地内（负 offset = 下游）
            out = 0; n = 0; over = 0.0
            for j, i in enumerate(idx):
                if not closes[j] or H[i] <= 0: continue
                p = pts[i] + td[j] * nrm[i]; n += 1
                if not parcel.contains(Point(p)):
                    out += 1; over = max(over, parcel.exterior.distance(Point(p)))
            V = float(np.nansum([0.5 * (areas[j] + areas[j + 1]) * 10 for j in range(len(areas) - 1) if closes[j] and closes[j + 1]]))
            rec.update({f"{ty}_closes_frac": float(closes[H[idx] > 0].mean()) if (H[idx] > 0).any() else 1.0,
                        f"{ty}_toe_us_max": float(np.nanmax(tu)) if np.isfinite(tu).any() else np.nan,
                        f"{ty}_toe_ds_max": float(np.nanmin(td)) if np.isfinite(td).any() else np.nan,
                        f"{ty}_base_max": float(np.nanmax(tu - td)) if np.isfinite(tu).any() else np.nan,
                        f"{ty}_vol_Mm3": V / 1e6,
                        f"{ty}_toe_ds_out_frac": (out / n) if n else 0.0, f"{ty}_toe_ds_over_m": over})
        rows.append(rec)
    fn = os.path.join(OUT, f"balanced_reaches_{LINE}_{crest_ft:.0f}.csv")
    with open(fn, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"\n== {LINE} crest {crest_ft:.0f} ft (nwl {crest_ft-20:.0f}) ==")
    print("段  桩号        长  H范围      CFRD:闭合 上趾 下趾 底宽 方量 出界%/m   | RCC:闭合 上趾 下趾 底宽 方量 出界%/m")
    for r in rows:
        print(f"{r['reach']:>2} {r['sta_from']:5.0f}-{r['sta_to']:5.0f} {r['length_m']:4.0f} {r['H_min']:5.0f}-{r['H_max']:4.0f}  "
              f"{r['cfrd_closes_frac']:4.2f} {r['cfrd_toe_us_max']:5.0f} {r['cfrd_toe_ds_max']:5.0f} {r['cfrd_base_max']:5.0f} {r['cfrd_vol_Mm3']:5.2f} {100*r['cfrd_toe_ds_out_frac']:3.0f}/{r['cfrd_toe_ds_over_m']:3.0f} | "
              f"{r['rcc_closes_frac']:4.2f} {r['rcc_toe_us_max']:5.0f} {r['rcc_toe_ds_max']:5.0f} {r['rcc_base_max']:5.0f} {r['rcc_vol_Mm3']:5.2f} {100*r['rcc_toe_ds_out_frac']:3.0f}/{r['rcc_toe_ds_over_m']:3.0f}")
    print("saved", fn)
