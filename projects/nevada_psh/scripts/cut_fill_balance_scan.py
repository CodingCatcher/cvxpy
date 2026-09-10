#!/usr/bin/env python3
"""目标"上面尽可能高 + 挖出来的合格石料就近堆成坝"的供需扫描（不选型，只给曲线与平衡点）。

供给侧：环内开挖到底面 F 的挖方（原岩方），按高程分成"≥3,950 ft（台地面，可能是 Mehrten 泥流盖层）"
        与"<3,950 ft（可能是 Calaveras 基底）"两份；考虑坝体上游坝趾的退让（RCC 30 m、CFRD 250 m、组合 150 m）。
需求侧：各断面策略的坝体压实方（网格法）换成原岩方（SF 0.95）；RCC 的骨料原岩方（0.9 m³/m³）。
库容：开挖后的毛库容（水位 4,100）减上游楔（取 types 汇总的值，随底面下降略有变化，作为近似）。
另：各站下游坝趾是否落在宗地界外。

用法：python3 cut_fill_balance_scan.py [A_parcel_ring_3700]
"""
from __future__ import annotations
import csv, json, os, sys
import numpy as np
from shapely.geometry import LineString, Point, Polygon
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from terrain_model import Terrain, FT, DATA, OUT, load_parcel, longitudinal_profile, inward_normals  # noqa

LINE = sys.argv[1] if len(sys.argv) > 1 else "A_parcel_ring_3700"
NWL = 4100 * FT; MEHRTEN_TOP_FT = 3950.0
SF, AGG = 0.95, 0.9
STRATS = {  # 压实方 Mm³（网格法）来自 types 汇总 / run 汇总
    "全环 RCC": {"setback": 30.0},
    "组合 II（陡坡 RCC + 其余 CFRD）": {"setback": 150.0},
    "RCC 高段(1-4) + 堆石边界段(5-12)": {"setback": 60.0},
    "全环 CFRD": {"setback": 250.0},
}
HIGH_REACH_END_M = 2695.0   # 段 1–4 的终点（沿 3,700 线的陡坡高段）

t = Terrain(os.path.join(DATA, "dem_3dep_1m_utm10.tif"), os.path.join(DATA, "dem_3dep_1m_utm10.json"))
cand = json.load(open(os.path.join(DATA, "dam_line_candidates.json")))[LINE]
ring = Polygon(cand["coords"]); line = LineString(cand["coords"])
types = json.load(open(os.path.join(OUT, f"types_{LINE}_summary.json")))["totals"]
combo = json.load(open(os.path.join(OUT, f"terrain_{LINE}_auto_close_summary.json")))
reach_rows = list(csv.DictReader(open(os.path.join(OUT, f"types_{LINE}_stations.csv"))))
def _vol(ty, lo, hi):
    A = [float(r[f"{ty}_area"]) if r[f"{ty}_area"] not in ("",) else float("nan") for r in reach_rows if lo <= float(r["station"]) <= hi]
    W = [float(r[f"{ty}_us_wedge"]) if r[f"{ty}_us_wedge"] not in ("",) else 0.0 for r in reach_rows if lo <= float(r["station"]) <= hi]
    V = sum(0.5 * (A[i] + A[i + 1]) * 10 for i in range(len(A) - 1) if np.isfinite(A[i]) and np.isfinite(A[i + 1]))
    Wv = sum(0.5 * (W[i] + W[i + 1]) * 10 for i in range(len(A) - 1))
    return V, Wv
_types_tot = json.load(open(os.path.join(OUT, f"types_{LINE}_summary.json")))["totals"]
_k_rcc = _types_tot["rcc"]["grid_m3"] / _types_tot["rcc"]["end_area_m3"]; _k_cf = _types_tot["cfrd"]["grid_m3"] / _types_tot["cfrd"]["end_area_m3"]
_v_rcc_hi, _w_rcc_hi = _vol("rcc", 0, HIGH_REACH_END_M); _v_cf_lo, _w_cf_lo = _vol("cfrd", HIGH_REACH_END_M + 1, 1e9)
demand = {
    "RCC 高段(1-4) + 堆石边界段(5-12)": {"rcc_m3": _v_rcc_hi * _k_rcc, "fill_m3": _v_cf_lo * _k_cf, "wedge": _w_rcc_hi + _w_cf_lo},
    "全环 RCC": {"rcc_m3": types["rcc"]["grid_m3"], "fill_m3": 0.0, "wedge": types["rcc"]["us_wedge_m3"]},
    "组合 II（陡坡 RCC + 其余 CFRD）": {"rcc_m3": combo["dam_volume_m3"]["by_kind"].get("rcc", 0) * combo["dam_volume_m3"]["grid_method"] / combo["dam_volume_m3"]["end_area_10m"],
                                  "fill_m3": combo["dam_volume_m3"]["by_kind"].get("cfrd", 0) * combo["dam_volume_m3"]["grid_method"] / combo["dam_volume_m3"]["end_area_10m"],
                                  "wedge": combo["storage"]["V_dam_inside_below_nwl_m3"]},
    "全环 CFRD": {"rcc_m3": 0.0, "fill_m3": types["cfrd"]["grid_m3"], "wedge": types["cfrd"]["us_wedge_m3"]},
}

win, m_ring = t.mask_inside(ring); z_ring = t.z[win]
st = list(csv.DictReader(open(os.path.join(OUT, f"types_{LINE}_stations.csv"))))
s_, pts, zg = longitudinal_profile(t, line, 10.0); nrm = inward_normals(pts, ring, None)
combo_st = {float(r["station_m"]): r["kind"] for r in csv.DictReader(open(os.path.join(OUT, f"terrain_{LINE}_auto_close_stations.csv")))}
def pit_polygon(strat):
    """库内可开挖区 = 上游坝趾连成的多边形（水面下坝体以内），再按切坡退让。"""
    pts_in = []
    for i, r in enumerate(st):
        if strat == "全环 RCC": ty = "rcc"
        elif strat == "全环 CFRD": ty = "cfrd"
        elif strat.startswith("RCC 高段"): ty = "rcc" if float(r["station"]) <= HIGH_REACH_END_M else "cfrd"
        else: ty = "rcc" if combo_st.get(float(r["station"]), "cfrd") == "rcc" else "cfrd"
        v = r[f"{ty}_toe_us"]; off = float(v) if v not in ("", None) else 0.0
        pts_in.append(pts[i] + max(off, 0.0) * nrm[i])
    pg = Polygon(pts_in).buffer(0)
    if pg.geom_type != "Polygon": pg = max(pg.geoms, key=lambda g: g.area)
    return pg
rows = []
floors = [4050, 4000, 3950, 3900, 3850, 3800, 3750, 3700, 3650, 3600, 3550, 3500]
X, Y = np.meshgrid(t.xc[win[1]], t.yc[win[0]])
import shapely
for strat, cfg in STRATS.items():
    pit0 = pit_polygon(strat)
    d = demand[strat]
    need_bank = d["fill_m3"] / SF + d["rcc_m3"] * AGG
    for F in floors:
        Fm = F * FT
        # 切坡退让：底面低于轴线最低地面时，按 0.75:1 切坡从坝趾内缘退让
        setback = max(0.0, (float(np.nanmin(zg)) - Fm)) * 0.75
        pit = pit0.buffer(-setback)
        if pit.is_empty: pit = pit0.buffer(-setback * 0.5)
        m_pit = shapely.contains_xy(pit, X.ravel(), Y.ravel()).reshape(X.shape)
        cut = np.where(m_pit, np.maximum(z_ring - Fm, 0.0), 0.0)
        cut_hi = float(cut[z_ring >= MEHRTEN_TOP_FT * FT].sum()); cut_lo = float(cut[z_ring < MEHRTEN_TOP_FT * FT].sum())
        # 台地盖层按厚度切分：像元高于 3,950 的部分算盖层，其下算基底
        cap = np.where(m_pit, np.clip(z_ring - MEHRTEN_TOP_FT * FT, 0, None) - np.clip(Fm - MEHRTEN_TOP_FT * FT, 0, None), 0.0)
        cap = np.clip(cap, 0, None); cap = np.minimum(cap, cut)
        cut_cap = float(cap.sum()); cut_base = float(cut.sum()) - cut_cap
        z_new = np.where(m_pit, np.minimum(z_ring, Fm), z_ring)
        gross = float(np.where(m_ring, np.maximum(NWL - z_new, 0.0), 0.0).sum())
        rows.append({"strategy": strat, "floor_ft": F, "setback_m": setback, "pit_area_acre": float(m_pit.sum() / 4046.86),
                     "cut_total_Mm3": cut.sum() / 1e6, "cut_cap_Mm3": cut_cap / 1e6, "cut_base_Mm3": cut_base / 1e6,
                     "need_bank_Mm3": need_bank / 1e6, "need_fill_bank_Mm3": d["fill_m3"] / SF / 1e6, "need_agg_bank_Mm3": d["rcc_m3"] * AGG / 1e6,
                     "balance_total": cut.sum() / need_bank if need_bank else float("nan"), "balance_base_only": cut_base / need_bank if need_bank else float("nan"),
                     "gross_Mm3": gross / 1e6, "net_Mm3": (gross - d["wedge"]) / 1e6, "net_GWh_518m": (gross - d["wedge"]) * 1000 * 9.81 * 518 * 0.85 / 3.6e12})
cols = list(rows[0].keys())
with open(os.path.join(OUT, f"cutfill_scan_{LINE}.csv"), "w") as f:
    f.write(",".join(cols) + "\n")
    for r in rows: f.write(",".join(f"{r[c]:.2f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
for r in rows:
    print(f"{r['strategy'][:10]:10s} F={r['floor_ft']} pit {r['pit_area_acre']:5.0f} ac cut {r['cut_total_Mm3']:6.1f} (cap {r['cut_cap_Mm3']:5.1f}, base {r['cut_base_Mm3']:5.1f}) need {r['need_bank_Mm3']:6.1f} bal {r['balance_total']:4.2f}/{r['balance_base_only']:4.2f} gross {r['gross_Mm3']:5.1f} net {r['net_Mm3']:5.1f} GWh {r['net_GWh_518m']:5.1f}")

# 下游坝趾出界统计（按交付物 1 的分段）
P = load_parcel()
reaches = list(csv.DictReader(open(os.path.join(OUT, "D1_reaches.csv"))))
out = []
for rc in reaches:
    a, b = float(rc["sta_from"]), float(rc["sta_to"]); rec = {"reach": int(rc["reach"]), "sta": f"{a:.0f}-{b:.0f}"}
    for ty in ("rcc", "cfrd", "ecrd"):
        n_out = 0; n = 0; over = 0.0
        for i, r in enumerate(st):
            sta = float(r["station"])
            if sta < a or sta > b: continue
            v = r[f"{ty}_toe_ds"]
            if v in ("", None): continue
            n += 1; p = Point(*(pts[i] + float(v) * nrm[i]))
            if not P.contains(p):
                n_out += 1; over = max(over, p.distance(P.exterior) if not P.contains(p) else 0.0)
        rec[f"{ty}_frac_outside"] = round(n_out / n, 2) if n else None; rec[f"{ty}_max_overshoot_m"] = round(over)
    out.append(rec)
for rec in out: print(rec)
json.dump(out, open(os.path.join(OUT, f"toe_outside_parcel_{LINE}.json"), "w"), indent=1)
