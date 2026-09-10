#!/usr/bin/env python3
"""交付物 1（坝型比较表）与交付物 2（整条坝线纵剖面与代表断面）——"打孔为主"框架（笔记 14）的 HTML 版。

框架：坝线固定（宗地内沿 3,700 ft 等高线闭合的环），坝顶高程与库底高程是变量，所有坝料来自上游坝趾以内的碗形坑
（坑壁顶线 = 上游坝趾内 20 m 平台，坑壁 0.75:1，底面平到库底）。
方案（坝顶由 outputs/balanced_<line>.csv 自动取）：
  A  组合 III：第 1、4 段 RCC，其余 CFRD，坝顶 = 孔养得起的最高坝顶
  B  全环 RCC，坝顶 4,120 ft（题设坝顶）
  C  全环 CFRD，坝顶 = 全堆石环养得起的最高坝顶
  R  参照：全环 CFRD，坝顶 4,120 ft（笔记 12/13 固定坝顶框架下的"缺口"）
每个方案两个库底：平衡点（挖方刚够坝料的最浅底面）与深挖参考点（净库容达几何极限 95% 的最浅底面）。
输入：outputs/balanced_A_parcel_ring_3700.csv（scripts/balanced_design_scan.py，碗形坑规则）、outputs/D1_reaches.csv（分段）、DEM 与坝线。
输出：deliverables/D1_dam_type_comparison.html、deliverables/D2_profile_and_sections.html（图 base64 内嵌，单文件可下载），
      outputs/D1v2_*.png/.csv、outputs/D2v2_*.png/.csv。
用法：python3 scripts/deliverables_1_2_html.py
"""
from __future__ import annotations
import base64, csv, json, math, os, sys
import numpy as np
import shapely
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from terrain_model import (Terrain, SECTIONS, FT, ACRE, DATA, OUT, Bowl, cross_section, dam_section,
                           inward_normals, longitudinal_profile, load_parcel)  # noqa
from deliverables_1_2 import geol, gravity_screen  # noqa
DEL = os.path.join(ROOT, "deliverables")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
matplotlib.rcParams["font.family"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]; matplotlib.rcParams["axes.unicode_minus"] = False
import logging; logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

LINE = "A_parcel_ring_3700"
LOWER_FT = 2400.0; ETA = 0.85
SF_ROCK = 1.2; LF_ROCK = 0.65; AGG_RCC = 0.9
CAP_FT = 3950.0; BERM = 20.0; WALL_M = 0.75
FLEET_ROCK_BCM_H = 5000.0; RCC_PLANT_M3_H = 700.0; HOURS_PER_YEAR = 4000.0
COMBO_RCC_STA = [(0.0, 606.0), (2302.8, 2695.2)]
COL = {"rcc": "#e69500", "cfrd": "#2e8b57"}
KIND_CN = {"rcc": "RCC 重力", "cfrd": "面板堆石 CFRD"}


def fnum(x):
    try:
        return float(x) if x not in ("", None) else float("nan")
    except Exception:
        return float("nan")


def gwh(v_m3, nwl_ft):
    return v_m3 * 1000 * 9.81 * (nwl_ft - LOWER_FT) * FT * ETA / 3.6e12


def station_kind(ty, st):
    if ty != "combo": return ty
    return "rcc" if any(a - 1e-6 <= st <= b + 1e-6 for a, b in COMBO_RCC_STA) else "cfrd"


# ---------------------------------------------------------------- 输入
t = Terrain(os.path.join(DATA, "dem_3dep_1m_utm10.tif"), os.path.join(DATA, "dem_3dep_1m_utm10.json"))
cand = json.load(open(os.path.join(DATA, "dam_line_candidates.json")))[LINE]
line = LineString(cand["coords"]); ring = Polygon(cand["coords"]); parcel = load_parcel()
CW = "顺时针" if not ring.exterior.is_ccw else "逆时针"
s, pts, zg = longitudinal_profile(t, line, 10.0); nrm = inward_normals(pts, ring, None)
secs = [cross_section(t, pts[i], nrm[i]) for i in range(len(pts))]
zmin_axis = float(np.nanmin(zg))
win, m_ring = t.mask_inside(ring); z_ring = t.z[win]
X, Y = np.meshgrid(t.xc[win[1]], t.yc[win[0]])
EXT_RING = [t.xc[win[1]][0] - 0.5, t.xc[win[1]][-1] + 0.5, t.yc[win[0]][-1] - 0.5, t.yc[win[0]][0] + 0.5]
reaches = list(csv.DictReader(open(os.path.join(OUT, "D1_reaches.csv"))))
for r in reaches:
    for k in r:
        if k != "geology": r[k] = fnum(r[k])
    r["reach"] = int(r["reach"])
bal = {(round(float(r["crest_ft"])), r["type"]): r for r in csv.DictReader(open(os.path.join(OUT, f"balanced_{LINE}.csv")))}
bal_uni = {(round(float(r["crest_ft"])), r["type"]): r for r in csv.DictReader(open(os.path.join(OUT, f"balanced_{LINE}_uniform.csv")))} if os.path.exists(os.path.join(OUT, f"balanced_{LINE}_uniform.csv")) else {}
srows = list(csv.DictReader(open(os.path.join(OUT, f"types_{LINE}_stations.csv"))))
slope_ds = np.array([fnum(r["slope_ds_deg"]) for r in srows]); slope_us = np.array([fnum(r["slope_us_deg"]) for r in srows])


def highest_balanced(ty):
    cs = [c for (c, tt), r in bal.items() if tt == ty and r["balanced"] == "True"]
    return max(cs) if cs else None

CREST_A = highest_balanced("combo"); CREST_C = highest_balanced("cfrd")
SCHEMES = [
    {"key": "A", "name": "方案 A：组合 III（1、4 段 RCC，其余 CFRD）", "crest": float(CREST_A), "type": "combo"},
    {"key": "B", "name": "方案 B：全环 RCC", "crest": 4120.0, "type": "rcc"},
    {"key": "C", "name": "方案 C：全环 CFRD", "crest": float(CREST_C), "type": "cfrd"},
    {"key": "R", "name": "参照 R：全环 CFRD，固定坝顶 4,120（笔记 12/13 框架）", "crest": 4120.0, "type": "cfrd"},
]
for sch in SCHEMES:
    br = bal[(round(sch["crest"]), sch["type"])]
    sch["floor_min"] = fnum(br["floor_ft"]) if br["balanced"] == "True" else None
    sch["floor_deep"] = fnum(br["maxnet_floor_ft"])
    sch["bal"] = br


def build(sch):
    """一个方案的逐站断面、坝体足迹、碗形坑。"""
    crest = sch["crest"] * FT; nwl = crest - 20 * FT
    n = len(pts); A = np.full(n, np.nan); W = np.zeros(n); TU = np.full(n, np.nan); TD = np.full(n, np.nan); ZTD = np.full(n, np.nan); KIND = []; D = []
    for i in range(n):
        k = station_kind(sch["type"], s[i]); KIND.append(k)
        d = dam_section(*secs[i], crest, SECTIONS[k], nwl); D.append(d)
        A[i] = d["area"]; W[i] = d["area_us_below_nwl"] if np.isfinite(d["area_us_below_nwl"]) else 0.0
        if d["toe_us"] is not None: TU[i] = d["toe_us"]
        if d["toe_ds"] is not None: TD[i] = d["toe_ds"]; ZTD[i] = d.get("z_toe_ds", np.nan)
    H = crest - zg; KIND = np.array(KIND); feas = np.isfinite(A)
    V_ea = float(sum(0.5 * (A[i] + A[i + 1]) * 10 for i in range(n - 1) if feas[i] and feas[i + 1]))
    wedge = float(sum(0.5 * (W[i] + W[i + 1]) * 10 for i in range(n - 1) if feas[i] and feas[i + 1]))
    toe_poly = Polygon([pts[i] + max(np.nan_to_num(TU[i], nan=0.0), 0.0) * nrm[i] for i in range(n)]).buffer(0)
    if toe_poly.geom_type != "Polygon": toe_poly = max(toe_poly.geoms, key=lambda g: g.area)
    bowl = Bowl(t, toe_poly, BERM, WALL_M); inside, wmax = bowl.field(X, Y, z_ring)
    # 坝体足迹：逐站四边形并集（鼻尖处内侧折线自相交也不会把库内填掉）
    quads = {"rcc": [], "cfrd": []}
    for i in range(n - 1):
        if all(H[j] > 0 and feas[j] and np.isfinite(TU[j]) and np.isfinite(TD[j]) for j in (i, i + 1)) and KIND[i] == KIND[i + 1]:
            q = Polygon([pts[i] + TD[i] * nrm[i], pts[i + 1] + TD[i + 1] * nrm[i + 1], pts[i + 1] + TU[i + 1] * nrm[i + 1], pts[i] + TU[i] * nrm[i]]).buffer(0)
            if not q.is_empty: quads[KIND[i]].append(q)
    foot = {k: unary_union(v) if v else None for k, v in quads.items()}
    return {"sch": sch, "crest": crest, "nwl": nwl, "A": A, "W": W, "TU": TU, "TD": TD, "ZTD": ZTD, "KIND": KIND, "H": H, "feas": feas, "D": D,
            "V_ea": V_ea, "wedge": wedge, "toe_poly": toe_poly, "bowl": bowl, "inside": inside, "wmax": wmax, "foot": foot,
            "dam_len": float((H > 0).sum() * 10), "H_max": float(np.nanmax(H))}


def excavate(b, floor_ft):
    Fm = floor_ft * FT
    z_new = np.where(b["inside"], np.minimum(z_ring, np.maximum(Fm, b["wmax"])), z_ring)
    dz = np.maximum(z_ring - z_new, 0.0); cut = float(np.nansum(dz))
    cap = float(np.nansum(np.maximum(z_ring - np.maximum(z_new, CAP_FT * FT), 0.0)))
    water = m_ring & (z_new < b["nwl"]); gross = float(np.nansum(np.where(m_ring, np.maximum(b["nwl"] - z_new, 0.0), 0.0)))
    return {"floor_ft": floor_ft, "cut": cut, "cap": cap, "base": cut - cap, "pit_acre": float((dz > 0.01).sum() * t.px ** 2 / ACRE),
            "gross": gross, "net": gross - b["wedge"], "water": water, "z_new": z_new, "dz": dz, "max_depth_m": float(np.nanmax(dz))}


def reach_rows(b):
    rows = []
    for r in reaches:
        a, c = r["sta_from"], r["sta_to"]; idx = np.where((s >= a - 1e-6) & (s <= c + 1e-6))[0]
        H = b["H"][idx]; kinds = b["KIND"][idx]; hasdam = H > 0; tu = b["TU"][idx]; td = b["TD"][idx]; A = b["A"][idx]
        out = 0; nn = 0; over = 0.0
        for j, i in enumerate(idx):
            if not (H[j] > 0 and np.isfinite(td[j])): continue
            p = pts[i] + td[j] * nrm[i]; nn += 1
            if not parcel.contains(Point(p)): out += 1; over = max(over, parcel.exterior.distance(Point(p)))
        closes = np.isfinite(A)
        rows.append({"reach": r["reach"], "sta_from": a, "sta_to": c, "length_m": float(len(idx) * 10), "H_min": float(H.min()), "H_max": float(H.max()),
                     "dam_len_m": float(hasdam.sum() * 10), "kind": "/".join(sorted(set(kinds[hasdam]))) if hasdam.any() else "无坝",
                     "closes_frac": float(closes[hasdam].mean()) if hasdam.any() else 1.0,
                     "toe_us_max": float(np.nanmax(tu)) if np.isfinite(tu).any() else np.nan, "toe_ds_max": float(np.nanmin(td)) if np.isfinite(td).any() else np.nan,
                     "base_max": float(np.nanmax(tu - td)) if np.isfinite(tu).any() else np.nan,
                     "vol_ea_Mm3": float(np.nansum([0.5 * (A[j] + A[j + 1]) * 10 for j in range(len(A) - 1) if closes[j] and closes[j + 1]])) / 1e6,
                     "toe_out_frac": (out / nn) if nn else 0.0, "toe_over_m": over})
    return rows


# ---------------------------------------------------------------- 计算
B = {sch["key"]: build(sch) for sch in SCHEMES}
B["A_cfrd"] = build({"key": "A_cfrd", "name": "", "crest": SCHEMES[0]["crest"], "type": "cfrd"})   # 用于分段表里"若此段用 CFRD/RCC"的对照
B["A_rcc"] = build({"key": "A_rcc", "name": "", "crest": SCHEMES[0]["crest"], "type": "rcc"})
EX = {}
for sch in SCHEMES:
    b = B[sch["key"]]
    EX[(sch["key"], "deep")] = excavate(b, sch["floor_deep"])
    if sch["floor_min"] is not None: EX[(sch["key"], "min")] = excavate(b, sch["floor_min"])
RR = {k: reach_rows(B[k]) for k in B}
H6 = max(20, round(RR["A"][5]["H_max"])) if RR["A"][5]["dam_len_m"] > 0 else 20   # 重力坝筛查最低取 20 m（更低的坝 FS 没有意义）
GS = {H: gravity_screen(H) for H in sorted({H6, round(B["A"]["H_max"]), 129})}


def scheme_row(sch, which):
    b = B[sch["key"]]; e = EX.get((sch["key"], which))
    if e is None: return None
    br = sch["bal"]; grid = fnum(br["dam_grid_Mm3"]); gc = fnum(br.get("grid_cfrd_Mm3", "")); gr = fnum(br.get("grid_rcc_Mm3", ""))
    if not np.isfinite(gc): gc, gr = (grid, 0.0) if sch["type"] == "cfrd" else (0.0, grid)
    need = fnum(br["need_bank_Mm3"]); net = e["net"] / 1e6; cut = e["cut"] / 1e6; nwl_ft = sch["crest"] - 20
    row = {"方案": sch["name"] + ("（深挖参考点）" if which == "deep" else "（最少开挖平衡点）"),
           "坝顶 / 水位 ft": f"{sch['crest']:.0f} / {nwl_ft:.0f}", "H 最大 m": round(b["H_max"]), "有坝长度 m": round(b["dam_len"]),
           "坝体 网格法 Mm³": f"{grid:.1f}（CFRD {gc:.1f} + RCC {gr:.1f}）" if gc > 0 and gr > 0 else f"{grid:.1f}",
           "坝体 平均断面法 Mm³": round(b["V_ea"] / 1e6, 1), "需原岩 Mm³": round(need, 1),
           "库底 ft": f"{e['floor_ft']:.0f}", "最大挖深 m": round(e["max_depth_m"]), "开挖面积 acre": round(e["pit_acre"]),
           "挖方 Mm³（盖层 + 基底）": f"{cut:.1f}（{e['cap']/1e6:.1f} + {e['base']/1e6:.1f}）",
           "挖方 ÷ 需料": round(cut / need, 2) if need > 0 else "—", "余料 Mm³ 原岩": round(cut - need, 1),
           "毛库容 Mm³": round(e["gross"] / 1e6, 1), "上游楔 Mm³": round(b["wedge"] / 1e6, 1), "净库容 Mm³": round(net, 1),
           f"GWh（下库 {LOWER_FT:.0f} ft，η {ETA}）": round(gwh(e["net"], nwl_ft), 1), "坝体(网格) / 净水": round(grid / net, 2) if net > 0 else "—",
           "堆石：压实方 / 原岩 / 松方 Mm³": f"{gc:.1f} / {gc/SF_ROCK:.1f} / {gc/SF_ROCK/LF_ROCK:.1f}" if gc > 0 else "—",
           "RCC：m³ / 骨料原岩 Mm³": f"{gr:.1f} / {gr*AGG_RCC:.1f}" if gr > 0 else "—",
           "堆石车队年数 @5,000 BCM/h": round(gc / SF_ROCK * 1e6 / FLEET_ROCK_BCM_H / HOURS_PER_YEAR, 1) if gc > 0 else "—",
           "RCC 拌合年数 @700 m³/h": round(gr * 1e6 / RCC_PLANT_M3_H / HOURS_PER_YEAR, 1) if gr > 0 else "—",
           "开挖年数 @5,000 BCM/h": round(cut * 1e6 / FLEET_ROCK_BCM_H / HOURS_PER_YEAR, 1)}
    ref_cut = fnum(br["maxnet_cut_Mm3"]) if which == "deep" else fnum(br["cut_Mm3"]); ref_net = fnum(br["maxnet_Mm3"]) if which == "deep" else fnum(br["net_Mm3"])
    row["_check"] = f"{sch['key']}/{which}: 扫描 {ref_cut:.1f}/{ref_net:.1f} vs 本脚本 {cut:.1f}/{net:.1f}"
    return row


SROWS = [r for sch in SCHEMES for which in ("min", "deep") for r in [scheme_row(sch, which)] if r]
scols = [c for c in SROWS[0].keys() if not c.startswith("_")]
with open(os.path.join(OUT, "D1v2_schemes.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=scols + ["_check"]); w.writeheader(); w.writerows(SROWS)
print("方案复核（挖方/净库容 Mm³）:"); [print("  ", r["_check"]) for r in SROWS]

# ---------------------------------------------------------------- 图 1.1 平面
def plan_panel(ax, key, which, picks=None):
    b = B[key]; e = EX[(key, which)]; sch = b["sch"]
    sl = t.slope_deg(); w2, _ = t.mask_inside(line.buffer(600))
    ext = [t.xc[w2[1]][0] - 0.5, t.xc[w2[1]][-1] + 0.5, t.yc[w2[0]][-1] - 0.5, t.yc[w2[0]][0] + 0.5]
    ax.imshow(sl[w2], cmap="Greys", extent=ext, vmin=0, vmax=55, alpha=0.75)
    ax.contour(t.xc, t.yc, t.z / FT, levels=np.arange(2600, 4400, 100), colors="#a08000", linewidths=0.3, alpha=0.7)
    # 开挖深度与水面
    dz = np.ma.masked_where(e["dz"] <= 0.01, e["dz"])
    im = ax.imshow(dz, cmap="YlOrRd", alpha=0.85, extent=EXT_RING, vmin=0, vmax=max(60, e["max_depth_m"]))
    wat = np.ma.masked_where(~e["water"], np.ones_like(z_ring))
    ax.imshow(wat, cmap=matplotlib.colors.ListedColormap(["#2c7fb8"]), alpha=0.35, extent=EXT_RING)
    if not b["bowl"].top.is_empty:
        tx, ty = b["bowl"].top.exterior.xy; ax.plot(tx, ty, "-", color="#d62728", lw=0.9)
    for kind in ("rcc", "cfrd"):
        fp = b["foot"][kind]
        if fp is None or fp.is_empty: continue
        for pp in ([fp] if fp.geom_type == "Polygon" else list(fp.geoms)):
            xx, yy = pp.exterior.xy; ax.fill(xx, yy, color=COL[kind], alpha=0.6, lw=0.5, edgecolor="k")
    ax.plot(*line.xy, "k-", lw=0.6)
    pxx, pyy = parcel.exterior.xy; ax.plot(pxx, pyy, "-", color="white", lw=1.6); ax.plot(pxx, pyy, "--", color="k", lw=0.8)
    for r in reaches:
        i = int(np.argmin(np.abs(s - 0.5 * (r["sta_from"] + r["sta_to"])))); p = pts[i] - 60 * nrm[i]
        ax.text(p[0], p[1], str(r["reach"]), fontsize=8, ha="center", va="center", bbox=dict(boxstyle="circle,pad=0.15", fc="white", ec="0.3", lw=0.5, alpha=0.9))
    if picks:
        for i, lab in picks:
            a_ = pts[i] - 300 * nrm[i]; b_ = pts[i] + 300 * nrm[i]
            ax.plot([a_[0], b_[0]], [a_[1], b_[1]], "-", color="#7b1fa2", lw=0.9); ax.text(a_[0], a_[1], f"S{lab}", fontsize=7, color="#7b1fa2", ha="right", va="center")
    x0, y0 = ext[0] + 120, ext[2] + 120; ax.plot([x0, x0 + 500], [y0, y0], "k-", lw=2); ax.text(x0 + 250, y0 + 25, "500 m", ha="center", fontsize=8)
    ax.annotate("N", xy=(ext[1] - 150, ext[3] - 120), xytext=(ext[1] - 150, ext[3] - 330), ha="center", fontsize=10, arrowprops=dict(arrowstyle="-|>", lw=1.2))
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{sch['name']}\n坝顶 {sch['crest']:.0f} / 水位 {sch['crest']-20:.0f} ft，库底 {e['floor_ft']:.0f} ft：挖 {e['cut']/1e6:.1f} Mm³（最深 {e['max_depth_m']:.0f} m），净库容 {e['net']/1e6:.1f} Mm³", fontsize=10)
    cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.01); cb.set_label("开挖深度 m", fontsize=8); cb.ax.tick_params(labelsize=7)
    ax.legend(handles=[Patch(color=COL["rcc"], alpha=0.7, label="RCC 坝体足迹（坝趾到坝趾）"), Patch(color=COL["cfrd"], alpha=0.7, label="CFRD 坝体足迹"),
                       Patch(facecolor="#fd8d3c", label="开挖深度（碗形坑，坑壁 0.75:1）"), Patch(facecolor="#2c7fb8", alpha=0.4, label="开挖后水面"),
                       plt.Line2D([], [], color="#d62728", lw=0.9, label="坑壁顶线（上游坝趾内 20 m）"),
                       plt.Line2D([], [], color="k", ls="--", lw=0.8, label="宗地界"), plt.Line2D([], [], color="k", lw=0.6, label="坝轴（3,700 ft 环）")], loc="lower right", fontsize=7)


def pick_stations(pairs):
    return [(int(np.argmin(np.abs(s - sta))), lab) for sta, lab in pairs]

PICKS_A = pick_stations([(80, "1"), (1298, "2"), (2253, "3"), (2432, "4"), (2775, "5"), (2910, "6"), (4560, "12")])
PICKS_B = pick_stations([(80, "1"), (1298, "2"), (2432, "4"), (2775, "5"), (2903, "6"), (3552, "8"), (4560, "12")])
fig, axs = plt.subplots(1, 2, figsize=(17.5, 9))
plan_panel(axs[0], "A", "deep", picks=PICKS_A); plan_panel(axs[1], "B", "deep", picks=PICKS_B)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "D1v2_plan.png"), dpi=110); plt.close(fig)

# ---------------------------------------------------------------- 图 1.2 供需柱
fig, ax = plt.subplots(figsize=(11.5, 5.2)); labels = []; xs = []; x = 0; ymax = 0
for sch in SCHEMES:
    for which in ("min", "deep"):
        e = EX.get((sch["key"], which))
        if e is None: continue
        br = sch["bal"]; gc = fnum(br.get("grid_cfrd_Mm3", "")); gr = fnum(br.get("grid_rcc_Mm3", ""))
        if not np.isfinite(gc): gc, gr = (fnum(br["dam_grid_Mm3"]), 0.0) if sch["type"] == "cfrd" else (0.0, fnum(br["dam_grid_Mm3"]))
        nf, na = gc / SF_ROCK, gr * AGG_RCC
        ax.bar(x - 0.2, nf, 0.38, color=COL["cfrd"], alpha=0.85); ax.bar(x - 0.2, na, 0.38, bottom=nf, color=COL["rcc"], alpha=0.85)
        ax.bar(x + 0.2, e["base"] / 1e6, 0.38, color="#7f7f7f"); ax.bar(x + 0.2, e["cap"] / 1e6, 0.38, bottom=e["base"] / 1e6, color="#c9a227")
        top = max(nf + na, e["cut"] / 1e6); ymax = max(ymax, top); ax.text(x, top + 1.5, f"挖÷需 {e['cut']/1e6/(nf+na):.2f}", ha="center", fontsize=8)
        labels.append(f"{sch['key']} {sch['crest']:.0f} ft\n底 {e['floor_ft']:.0f}（{'平衡' if which=='min' else '深挖'}）"); xs.append(x); x += 1
ax.set_ylim(0, ymax * 1.12); ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=8); ax.set_ylabel("Mm³ 原岩")
ax.legend(handles=[Patch(color=COL["cfrd"], label="需料：堆石原岩（压实方 ÷ 1.2）"), Patch(color=COL["rcc"], label="需料：RCC 骨料原岩（RCC × 0.9）"),
                   Patch(color="#7f7f7f", label="孔能出：基底岩（< 3,950 ft）"), Patch(color="#c9a227", label="孔能出：Mehrten 盖层（≥ 3,950 ft，能否作坝料待试验）")], fontsize=8, loc="upper left")
ax.set_title("每个方案的坝料需求（左柱）与碗形坑能出的料（右柱）", fontsize=10); ax.grid(axis="y", alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "D1v2_balance_bars.png"), dpi=120); plt.close(fig)

# ---------------------------------------------------------------- 图 1.4 坑的规则对比（旧统一退让 vs 碗形）
if bal_uni:
    fig, ax = plt.subplots(figsize=(10, 4.6))
    for ty, col, nm in (("cfrd", COL["cfrd"], "全环 CFRD"), ("combo", "#7b1fa2", "组合 III"), ("rcc", COL["rcc"], "全环 RCC")):
        cs = sorted(c for (c, tt) in bal if tt == ty)
        ax.plot(cs, [fnum(bal[(c, ty)]["need_bank_Mm3"]) for c in cs], "-", color=col, lw=1.6, label=f"{nm} 需料")
        ax.plot(cs, [fnum(bal[(c, ty)]["max_cut_Mm3"]) for c in cs], "o--", color=col, lw=1, ms=4, label=f"{nm} 孔的极限挖方（碗形坑）")
        cu = sorted(c for (c, tt) in bal_uni if tt == ty)
        ax.plot(cu, [fnum(bal_uni[(c, ty)]["max_cut_Mm3"]) for c in cu], ":", color=col, lw=1, alpha=0.7, label=f"{nm} 极限挖方（旧：统一退让）")
    ax.set_xlabel("坝顶高程 ft"); ax.set_ylabel("Mm³ 原岩"); ax.grid(alpha=0.3); ax.legend(fontsize=7, ncol=3, loc="upper left"); ax.set_ylim(0, 140)
    ax.set_title("坝料需求 vs 孔的极限挖方随坝顶高程的变化：实线在虚线之下才养得起", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "D1v2_need_vs_cut.png"), dpi=120); plt.close(fig)

# ---------------------------------------------------------------- 分段表（A、B）
def reach_text_A(k, ra, rc, rr):
    """k 段号；ra 方案 A 的段行；rc/rr 同坝顶下全环 CFRD / 全环 RCC 的段行。"""
    f0 = lambda v: f"{v:.0f}" if np.isfinite(v) else "—"
    if k == 1:
        return ("RCC", f"西坡下游 {reaches[0]['slope_ds_mean']:.0f}°：CFRD 1.4:1 在此底宽 {f0(rc['base_max'])} m、{rc['vol_ea_Mm3']:.1f} Mm³（占全 CFRD 环 {100*rc['vol_ea_Mm3']/sum(x['vol_ea_Mm3'] for x in RR['A_cfrd']):.0f}%）；RCC 底宽 {f0(rr['base_max'])} m、{rr['vol_ea_Mm3']:.1f} Mm³，骨料 {rr['vol_ea_Mm3']*AGG_RCC:.1f} Mm³ 原岩从孔里来")
    if k == 2:
        return ("CFRD", f"北鼻尖 {reaches[1]['slope_ds_mean']:.0f}°，离 1.4:1 的闭合极限 35.5° 够远，下游趾 {f0(-rc['toe_ds_max'])} m 在界内；{rc['vol_ea_Mm3']:.1f} Mm³ 是孔料的主要去处，上游坝趾以内就是坑")
    if k == 3:
        return ("CFRD", f"地面 {reaches[2]['slope_ds_mean']:.0f}° 但坡脚有台阶，CFRD 趾 {f0(-rc['toe_ds_max'])} m；若要再抬坝顶，此段先改 RCC（省 {(rc['vol_ea_Mm3']/SF_ROCK - rr['vol_ea_Mm3']*AGG_RCC):.1f} Mm³ 原岩）")
    if k == 4:
        return ("RCC", f"东北坡 31–42°：CFRD 坝趾 {f0(-rc['toe_ds_max'])} m，{100*rc['toe_out_frac']:.0f}% 桩号出宗地北界最多 {rc['toe_over_m']:.0f} m；RCC 坝趾 {f0(-rr['toe_ds_max'])} m 全在界内。不用 RCC 就得把坝轴向库内退 200 m 以上")
    if k == 5:
        return ("CFRD，坝轴内移 80–100 m", f"沿东界的过渡段，两种坝型的下游坝趾都出界（CFRD {rc['toe_over_m']:.0f} m、RCC {rr['toe_over_m']:.0f} m）；料 {rc['vol_ea_Mm3']:.1f} Mm³ 从孔里来")
    if k == 6:
        if ra["dam_len_m"] > 0:
            return (f"前 {ra['dam_len_m']:.0f} m 低坝（H ≤ {ra['H_max']:.0f} m），其余天然山体", f"地面 3,867 ft 起低于坝顶，是题设\"低段可用混凝土、按挡水坝分析\"的段落（重力坝筛查 FS {GS[H6][0]:.1f} @ {H6} m）；坝趾出东界 {max(rc['toe_over_m'], rr['toe_over_m']):.0f} m，要内移")
        return ("无坝", "地面 3,867–4,019 ft 高于坝顶，天然山体挡水")
    if k in (7, 8, 9, 10, 11):
        r = reaches[k - 1]; return ("无坝", f"地面 {r['z_ground_ft_min']:.0f}–{r['z_ground_ft_max']:.0f} ft 高于坝顶，天然山体挡水")
    if k == 12:
        return ("RCC 或改线", f"沿宗地西界顺坡而下的闭合段，任何坝型的下游坝趾都在界外（RCC {rr['toe_over_m']:.0f} m、CFRD {rc['toe_over_m']:.0f} m）；方量小（RCC {rr['vol_ea_Mm3']:.2f} Mm³）；要么取得西侧的带状用地，要么在宗地内沿更高一条等高线回折闭合")

def cell(x):
    if x["dam_len_m"] <= 0: return "无坝"
    base = f"{x['base_max']:.0f}" if np.isfinite(x["base_max"]) else "—"; tds = f"{-x['toe_ds_max']:.0f}" if np.isfinite(x["toe_ds_max"]) else "—"
    return f"H {max(x['H_min'],0):.0f}–{x['H_max']:.0f} m；{x['kind'].upper()}；下游趾 {tds} m，底宽 {base} m；{x['vol_ea_Mm3']:.2f} Mm³；出界 {x['toe_out_frac']*100:.0f}%/{x['toe_over_m']:.0f} m"

RROWS = []
for r, ra, rb, rc, rr in zip(reaches, RR["A"], RR["B"], RR["A_cfrd"], RR["A_rcc"]):
    k = r["reach"]; kind, why = reach_text_A(k, ra, rc, rr)
    if rb["dam_len_m"] <= 0: tb = "无坝（地面高于坝顶）"
    elif rb["toe_out_frac"] < 0.05: tb = "RCC；坝趾在界内"
    else: tb = "RCC；下游坝趾 %.0f%% 出界最多 %.0f m → 内移/改线/用地" % (100 * rb["toe_out_frac"], rb["toe_over_m"])
    RROWS.append({"段": k, "桩号 m": f"{r['sta_from']:.0f}–{r['sta_to']:.0f}", "长 m": f"{r['length_m']:.0f}", "地面 ft": f"{r['z_ground_ft_min']:.0f}–{r['z_ground_ft_max']:.0f}",
                  "下游坡°": f"{r['slope_ds_mean']:.0f}", "库侧坡°": f"{r['slope_us_mean']:.0f}", "地质（可能）": r["geology"],
                  f"方案 A（坝顶 {CREST_A:.0f}）": cell(ra), "方案 A 坝型": kind, "方案 A 理由 / 条件": why,
                  f"同坝顶若全段 CFRD": cell(rc) if rc["dam_len_m"] > 0 else "无坝", f"同坝顶若全段 RCC": cell(rr) if rr["dam_len_m"] > 0 else "无坝",
                  "方案 B（坝顶 4,120，全 RCC）": cell(rb), "方案 B 说明": tb})
with open(os.path.join(OUT, "D1v2_reaches.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(RROWS[0].keys())); w.writeheader(); w.writerows(RROWS)

# ---------------------------------------------------------------- 图 2.1 纵剖面
def profile_figure():
    bA, bB = B["A"], B["B"]; zft = zg / FT; cA = bA["crest"] / FT; cB = bB["crest"] / FT
    fA, fB = SCHEMES[0]["floor_deep"], SCHEMES[1]["floor_deep"]
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(17, 9.5), gridspec_kw={"height_ratios": [3, 1.2]}, sharex=True)
    colors = {"Calaveras 千枚岩/泥质岩（可能）": "#d9d2e9", "Mehrten 火山泥流（可能）": "#f4cccc", "变质火山岩（可能）": "#cfe2f3"}
    for r in reaches:
        ax.axvspan(r["sta_from"], r["sta_to"], color=colors[r["geology"]], alpha=0.45, lw=0); ax.axvline(r["sta_from"], color="0.4", lw=0.5, ls=":")
        ax.text((r["sta_from"] + r["sta_to"]) / 2, 4160, str(r["reach"]), ha="center", fontsize=8)
    for kind in ("rcc", "cfrd"):
        m = (bA["H"] > 0) & (bA["KIND"] == kind)
        ax.fill_between(s, zft, cA, where=m, color=COL[kind], alpha=0.6, lw=0, label=f"方案 A 坝体：{KIND_CN[kind]}")
    ax.fill_between(s, np.maximum(zft, cA), cB, where=bB["H"] > 0, color=COL["rcc"], alpha=0.25, hatch="..", lw=0, label=f"方案 B 在 A 之上再加的 RCC 坝体（{cA:.0f}→4,120）")
    ax.plot(s, zft, "k-", lw=1.3, label="轴线地面（1 m DEM）")
    ax.axhline(cA, color="#7b1fa2", lw=1.2, label=f"方案 A 坝顶 {cA:.0f} ft（水位 {cA-20:.0f}）"); ax.axhline(cA - 20, color="#7b1fa2", lw=0.7, ls="--")
    ax.axhline(cB, color="#d62728", lw=1.2, label="方案 B 坝顶 4,120 ft（水位 4,100）"); ax.axhline(cB - 20, color="#d62728", lw=0.7, ls="--")
    ax.axhline(fA, color="#7b1fa2", lw=0.8, ls=":", label=f"方案 A 库底 {fA:.0f} ft（深挖参考点）"); ax.axhline(fB, color="#d62728", lw=0.8, ls=":", label=f"方案 B 库底 {fB:.0f} ft（深挖参考点）")
    ax.plot(s, bA["ZTD"] / FT, ".", color="#1b5e20", ms=3, label="方案 A 下游坝趾高程"); ax.plot(s, bB["ZTD"] / FT, "x", color="#b71c1c", ms=3, label="方案 B 下游坝趾高程（RCC 0.8:1）")
    for name, c in colors.items(): ax.fill_between([], [], [], color=c, alpha=0.45, label=name)
    ax.set_ylabel("高程 ft"); ax.set_ylim(3050, 4280); ax.grid(alpha=0.3); ax.legend(loc="lower left", fontsize=7.5, ncol=3)
    ax.set_title(f"坝线 {LINE} 纵剖面（L = {s[-1]:.0f} m，{CW}）：地面、两个方案的坝顶/水位/库底、坝体与下游坝趾高程、分段号与可能地质", fontsize=10)
    ax2.plot(s, bA["H"], "-", color="#7b1fa2", lw=1, label=f"H（坝顶 {cA:.0f} − 地面）m"); ax2.plot(s, bB["H"], "-", color="#d62728", lw=1, label="H（坝顶 4,120 − 地面）m")
    ax2.plot(s, slope_ds, "m-", lw=0.6, alpha=0.8, label="下游侧地面坡角°（100 m 内）"); ax2.plot(s, slope_us, "c-", lw=0.6, alpha=0.8, label="库侧地面坡角°")
    ax2.axhline(35.5, color="green", lw=0.6, ls=":", label="1.4:1 面的坡角 35.5°（地面更陡则 CFRD 永不闭合）"); ax2.axhline(0, color="k", lw=0.4)
    ax2.set_xlabel(f"桩号 m（从南界西端起{CW}：西坡 → 北鼻尖 → 东北坡 → 东界 → 南台地 → 西南闭合）"); ax2.set_ylabel("m / °"); ax2.grid(alpha=0.3); ax2.legend(loc="upper right", fontsize=7.5, ncol=2)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "D2v2_profile.png"), dpi=110); plt.close(fig)

profile_figure()

# ---------------------------------------------------------------- 图 2.2/2.3 代表断面
def section_figure(key, picks, fname):
    b = B[key]; sch = b["sch"]; e = EX[(key, "deep")]; crest_ft = sch["crest"]; nwl_ft = crest_ft - 20; floor_ft = e["floor_ft"]
    n = len(picks); fig, axs = plt.subplots(n, 1, figsize=(14, 3.4 * n)); axs = np.atleast_1d(axs); table = []
    for ax, (i, lab) in zip(axs, picks):
        off, zsec = secs[i]; zft = zsec / FT; kind = b["KIND"][i]; d = b["D"][i]; H = b["H"][i]
        P = pts[i][None, :] + off[:, None] * nrm[i][None, :]
        znew = b["bowl"].surface(P, zsec, floor_ft * FT)
        # 只画库侧（+）的开挖；下游侧若断面穿过冲沟再进入另一叶库区，那里的坑不在本断面讨论
        znew = np.where(off >= 0, znew, zsec)
        exc = znew < zsec - 1e-6
        ax.fill_between(off, znew / FT, zft, where=exc, color="#f4a582", alpha=0.75, lw=0, label=f"开挖到库底 {floor_ft:.0f} ft（坑壁 0.75:1，平台 20 m）")
        ax.plot(off, zft, "k-", lw=1.2, label="地面"); ax.plot(off, znew / FT, "-", color="#b2182b", lw=0.8)
        bw = SECTIONS[kind]["crest_w"]
        if H > 0 and d["toe_us"] is not None and d["toe_ds"] is not None:
            xs = [d["toe_ds"], -bw / 2, bw / 2, d["toe_us"]]; zs = [d["z_toe_ds"] / FT, crest_ft, crest_ft, d["z_toe_us"] / FT]
            m = (off >= d["toe_ds"]) & (off <= d["toe_us"]); gx = off[m]; gz = zft[m]
            ax.fill(list(xs) + list(gx[::-1]), list(zs) + list(gz[::-1]), color=COL[kind], alpha=0.65, lw=0.8, edgecolor="k",
                    label=f"{kind.upper()}：底宽 {d['base_width']:.0f} m，面积 {d['area']:.0f} m²，上游楔 {d['area_us_below_nwl']:.0f} m²")
            x_start = bw / 2
        else:
            x_start = 0.0
            if H <= 0: ax.text(0, crest_ft + 25, "地面高于坝顶：此处无坝，天然山体挡水", ha="center", fontsize=9, color="0.25")
        # 水：从坝的上游面起连续到第一次开挖面高于水位处
        j = int(np.argmin(np.abs(off - x_start))); wz = znew / FT
        k2 = j
        while k2 < len(off) and (wz[k2] < nwl_ft or (H > 0 and off[k2] <= (d["toe_us"] or 0))): k2 += 1
        if k2 > j: ax.fill_between(off[j:k2], np.minimum(wz[j:k2], nwl_ft), nwl_ft, color="#2c7fb8", alpha=0.25, lw=0, label="水")
        ax.axhline(nwl_ft, color="b", ls="--", lw=0.7, label=f"水位 {nwl_ft:.0f} ft"); ax.axhline(crest_ft, color="r", lw=0.6, label=f"坝顶 {crest_ft:.0f} ft")
        out_ds = None
        for jj in range(len(off)):
            if off[jj] < 0 and not parcel.contains(Point(P[jj])): out_ds = off[jj]
        if out_ds is not None:
            ax.axvline(out_ds, color="k", ls="-.", lw=0.7); ax.text(out_ds, crest_ft + 25, "宗地界", fontsize=7, ha="center")
        cut_area = float(np.nansum(np.maximum(zsec - znew, 0)) * (off[1] - off[0])); depth = float(np.nanmax(zsec - znew))
        rec = {"scheme": key, "station": float(s[i]), "reach": lab, "E": float(pts[i][0]), "N": float(pts[i][1]), "z_ground_ft": float(zg[i] / FT), "H_m": float(H),
               "slope_ds_deg": float(slope_ds[i]), "slope_us_deg": float(slope_us[i]), "geology": geol(pts[i][0], pts[i][1], zg[i] / FT),
               "kind": kind if H > 0 else "无坝", "area_m2": d["area"] if H > 0 else 0.0, "base_m": d.get("base_width", np.nan) if H > 0 else np.nan,
               "toe_us_m": d["toe_us"] if H > 0 else np.nan, "toe_ds_m": d["toe_ds"] if H > 0 else np.nan, "z_toe_ds_ft": d.get("z_toe_ds", np.nan) / FT if H > 0 else np.nan,
               "us_wedge_m2": d["area_us_below_nwl"] if H > 0 else 0.0, "cut_area_m2": cut_area, "cut_depth_max_m": depth, "parcel_edge_ds_m": out_ds if out_ds is not None else np.nan}
        table.append(rec)
        lo = np.nanmin(np.where(np.isfinite(znew), znew, np.inf)) / FT
        ax.set_xlim(-900, 600); ax.set_ylim(min(lo, crest_ft - 650) - 30, crest_ft + 70)
        ax.set_title(f"S{lab}  桩号 {s[i]:.0f} m（段 {lab}）：地面 {zg[i]/FT:.0f} ft，H = {H:.0f} m，下游坡 {slope_ds[i]:.0f}°，库侧坡 {slope_us[i]:.0f}°，{rec['geology']}；断面开挖 {cut_area:.0f} m²，最深 {depth:.0f} m", fontsize=9)
        ax.set_xlabel("离坝轴偏移 m（+ 为库侧）"); ax.set_ylabel("ft"); ax.grid(alpha=0.3); ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, fname), dpi=105); plt.close(fig)
    return table

TAB_A = section_figure("A", PICKS_A, "D2v2_sections_A.png"); TAB_B = section_figure("B", PICKS_B, "D2v2_sections_B.png")
with open(os.path.join(OUT, "D2v2_sections_table.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(TAB_A[0].keys())); w.writeheader(); w.writerows(TAB_A + TAB_B)

# ---------------------------------------------------------------- HTML
def img64(path): return "data:image/png;base64," + base64.b64encode(open(path, "rb").read()).decode()
def fig_html(path, caption): return f'<figure><img src="{img64(path)}" alt="{caption}"><figcaption>{caption}</figcaption></figure>'
def esc(x): return str(x).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
def table_html(cols, rows):
    h = ['<div class="wrap"><table><thead><tr>' + "".join(f"<th>{esc(c)}</th>" for c in cols) + "</tr></thead><tbody>"]
    for r in rows:
        h.append("<tr>" + "".join(f"<td>{esc(r.get(c, '') if isinstance(r, dict) else r[j])}</td>" for j, c in enumerate(cols)) + "</tr>")
    h.append("</tbody></table></div>"); return "\n".join(h)
CSS = """
body{font-family:"Noto Sans CJK SC","WenQuanYi Zen Hei","PingFang SC","Microsoft YaHei","Helvetica Neue",Arial,sans-serif;max-width:1500px;margin:0 auto;padding:18px 26px 60px;line-height:1.55;color:#1f2328;background:#fff}
h1{font-size:22px;border-bottom:2px solid #333;padding-bottom:6px} h2{font-size:18px;margin-top:34px;border-left:5px solid #7b1fa2;padding-left:10px} h3{font-size:15px;margin-top:22px}
table{border-collapse:collapse;font-size:12.5px;margin:8px 0 14px;background:#fff} th,td{border:1px solid #bfc4cc;padding:4px 7px;vertical-align:top;text-align:left}
th{background:#eef1f5;position:sticky;top:0;z-index:1} tbody tr:nth-child(even) td{background:#f8f9fb} td:first-child{white-space:nowrap}
.wrap{overflow-x:auto;max-width:100%} figure{margin:14px 0} img{max-width:100%;height:auto;border:1px solid #d8dbe0;background:#fff}
figcaption{font-size:12.5px;color:#444;margin-top:4px} .note{background:#fff8e1;border-left:4px solid #f0b400;padding:8px 14px;margin:10px 0}
.key{background:#f3e5f5;border-left:4px solid #7b1fa2;padding:8px 14px;margin:10px 0} .fix{background:#fdecea;border-left:4px solid #d62728;padding:8px 14px;margin:10px 0}
.small{font-size:12px;color:#555} ul{margin-top:4px} li{margin:3px 0} code{background:#f1f3f5;padding:1px 4px;border-radius:3px;font-size:12px} nav a{margin-right:14px}
"""
def html_doc(title, body):
    return f'<!DOCTYPE html>\n<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>{esc(title)}</title><style>{CSS}</style></head><body>\n{body}\n</body></html>'
fmt = lambda v: (f"{v:.0f}" if isinstance(v, float) and np.isfinite(v) else ("—" if isinstance(v, float) else str(v)))
eA, eB, eC = EX[("A", "deep")], EX[("B", "deep")], EX[("C", "deep")]
bA, bB, bC, bR = (SCHEMES[i]["bal"] for i in range(4))
uA = bal_uni.get((round(CREST_A), "combo")); uC = bal_uni.get((round(CREST_C), "cfrd")); uB = bal_uni.get((4120, "rcc"))

# ---- D1
d1 = []
d1.append(f"<h1>交付物 1：坝型比较表（打孔为主框架）——Nevada County PSH 上库，坝线 {LINE}</h1>")
d1.append('<nav><a href="#s0">0 平面</a><a href="#s1">1.1 方案级比较</a><a href="#s2">1.2 料的供需</a><a href="#s3">1.3 准则对照</a><a href="#s4">1.4 逐段表</a><a href="#s5">1.5 三个环</a><a href="#s6">1.6 未回答的</a><a href="#s7">复算</a></nav>')
d1.append('<div class="key"><b>框架（笔记 14）。</b>上库是在天然地里打出来的孔，所有坝料来自这个孔。坝线固定为宗地内沿 3,700 ft 等高线闭合的环（4,589 m，环内 139 acre），<b>坝顶高程与库底高程是变量</b>。'
          '坑是碗形的：坑壁顶线在上游坝趾（无坝处为坝轴/宗地界）以内 20 m 的平台边，坑壁以 0.75:1 向内下降，底面平到库底；要求挖方 ≥ 坝体需料。'
          f'在此约束下比较：<b>A</b> 组合 III（第 1、4 段 RCC，其余 CFRD）坝顶 {CREST_A:.0f} ft——孔养得起的最高坝顶；<b>B</b> 全环 RCC 坝顶 4,120 ft；<b>C</b> 全环 CFRD 坝顶 {CREST_C:.0f} ft——全堆石环养得起的最高坝顶；参照 <b>R</b> 全环 CFRD 固定坝顶 4,120（笔记 12/13 里"料缺 86%"的那个）。'
          '每个方案给两个库底：<i>最少开挖平衡点</i>（挖方刚够坝料的最浅底面）与<i>深挖参考点</i>（净库容达到几何极限——坑壁在中间相遇——95% 的最浅底面）。</div>')
if uA and uC and uB:
    d1.append(f'<div class="fix"><b>对笔记 14 初版数字的修正。</b>初版把可挖区算成"坝趾多边形统一向内退 (最低坝基 − 库底) × 0.75 后挖平"，在地面高的段落（东界、南台地、CFRD 上游坝趾以内的坡）坑壁被算到了坝体下面和宗地外，挖方偏大。'
              f'改为碗形坑后，孔的极限挖方从 {fnum(uA["max_cut_Mm3"]):.0f} → {fnum(bA["max_cut_Mm3"]):.0f} Mm³（组合 III，坝顶 {CREST_A:.0f}）、{fnum(uC["max_cut_Mm3"]):.0f} → {fnum(bC["max_cut_Mm3"]):.0f}（全 CFRD，坝顶 {CREST_C:.0f}）、{fnum(uB["max_cut_Mm3"]):.0f} → {fnum(bB["max_cut_Mm3"]):.0f}（全 RCC，坝顶 4,120）；'
              f'堆石环养得起的坝顶相应下移（组合 III 3,950 → {CREST_A:,.0f}；全 CFRD 仍 {CREST_C:,.0f} 但余料由 ×{fnum(uC["maxnet_surplus_ratio"]):.1f} 降到 ×{fnum(bC["maxnet_surplus_ratio"]):.1f}），全 RCC 到 4,120 由宽松变成勉强（极限挖方÷需料 {fnum(bB["max_cut_Mm3"])/fnum(bB["need_bank_Mm3"]):.2f}，要挖到 {EX[("B","min")]["floor_ft"]:.0f} ft）。见图 1.4；本文件所有数字都是修正后的。</div>')
d1.append('<h2 id="s0">0. 方案在地面上是什么样</h2>')
d1.append(fig_html(os.path.join(OUT, "D1v2_plan.png"), "图 1.1 方案 A（左）与方案 B（右）在深挖参考点的平面：坝体足迹按坝型着色（坝趾到坝趾），黄–红色阶是碗形坑的开挖深度，红线是坑壁顶线，蓝色是开挖后的水面；底图坡度灰阶 + 100 ft 等高线。方案 A 的东界与南台地没有坝，水被天然山体挡在宗地内。段号对应 1.4 节；S# 是交付物 2 代表断面的位置。"))
d1.append('<h2 id="s1">1.1 方案级比较</h2>')
d1.append(table_html(scols, SROWS))
d1.append(f'<p class="small">网格法方量计入平面曲率（环外侧的扇形），平均断面法不计；需料按网格法取。三态换算：堆石 1 m³ 原岩 → 1.2 m³ 压实方（松方 = 原岩 ÷ 0.65）；RCC 每 m³ 需骨料原岩 0.9 m³。'
          f'"盖层"= 坑内 3,950 ft 以上的 Mehrten 泥流（能否作坝料待试验采场）。车队与拌合站只是把方量换成年的尺子（5,000 BCM/h ≈ 8 台 D10T2 松土 + 4 台大挖 + 40 台 740C；700 m³/h ≈ 一座 1,000 yd³/h 拌合站），不是配置建议。'
          f'GWh 按下库 2,400 ft（北界河谷）、η 0.85，只作尺度。复核：本脚本按同一规则重挖一遍与 balanced 扫描一致（{"; ".join(r["_check"] for r in SROWS)}）。</p>')
d1.append('<div class="note"><b>读法。</b>'
          f'A：坝顶 {CREST_A:.0f} 时孔养得起堆石环（深挖参考点挖÷需 {eA["cut"]/1e6/fnum(bA["need_bank_Mm3"]):.2f}），净库容 {eA["net"]/1e6:.1f} Mm³（{gwh(eA["net"], CREST_A-20):.0f} GWh）；RCC 只占坝体 {100*fnum(bA["grid_rcc_Mm3"])/fnum(bA["dam_grid_Mm3"]):.0f}%，却换掉了 CFRD 最费料的两段陡坡。'
          + (f'B：全环 RCC 到坝顶 4,120 只是勉强平衡——孔的极限挖方 {fnum(bB["max_cut_Mm3"]):.1f} Mm³ 对需料 {fnum(bB["need_bank_Mm3"]):.1f} Mm³（×{fnum(bB["max_cut_Mm3"])/fnum(bB["need_bank_Mm3"]):.2f}），要挖到 {EX[("B","min")]["floor_ft"]:.0f} ft 才够，净库容 {EX[("B","min")]["net"]/1e6:.1f} Mm³（{gwh(EX[("B","min")]["net"],4100):.0f} GWh）；深挖参考点 {eB["floor_ft"]:.0f} ft 的净库容 {eB["net"]/1e6:.1f} Mm³ 但料差 {100*(1-eB["cut"]/1e6/fnum(bB["need_bank_Mm3"])):.0f}%。代价还有 {fnum(bB["dam_grid_Mm3"]):.0f} Mm³ RCC 的胶材外运、拌合、温控与 129 m 直立面的抗滑；坝顶 4,050–4,100 之间 RCC 的余料比就宽松得多。' if ("B","min") in EX else f'B：全环 RCC 在 4,120 料不够（极限挖方 {fnum(bB["max_cut_Mm3"]):.1f} vs 需 {fnum(bB["need_bank_Mm3"]):.1f} Mm³）。')
          + f'C：全堆石环只能到坝顶 {CREST_C:.0f}，净库容 {eC["net"]/1e6:.1f} Mm³（{gwh(eC["net"], CREST_C-20):.0f} GWh），且第 4 段坝趾仍出界。R：固定坝顶 4,120 的全 CFRD 需料 {fnum(bR["need_bank_Mm3"]):.0f} Mm³，孔最多出 {fnum(bR["max_cut_Mm3"]):.0f} Mm³——这就是笔记 13 的"缺口"，它是框架的产物。</div>')
d1.append('<h2 id="s2">1.2 料的供需与坝顶高程</h2>')
d1.append(fig_html(os.path.join(OUT, "D1v2_balance_bars.png"), "图 1.2 每个方案两个库底下的需料（堆石原岩 + RCC 骨料）与碗形坑能出的料（基底 + 盖层）。柱上数字是挖方÷需料。"))
d1.append(fig_html(os.path.join(OUT, f"balanced_{LINE}.png"), f"图 1.3 坝顶高程扫描（scripts/balanced_design_scan.py，碗形坑）：三种断面策略在每个坝顶下的净库容电量。实心 = 孔养得起，空心 = 料不够。全 CFRD 上限 {CREST_C:.0f} ft，组合 III {CREST_A:.0f} ft，全 RCC 到 4,120 ft。"))
if bal_uni:
    d1.append(fig_html(os.path.join(OUT, "D1v2_need_vs_cut.png"), "图 1.4 坝料需求（实线）与孔的极限挖方（虚线 = 碗形坑；点线 = 初版的统一退让规则）随坝顶高程的变化。堆石需求随 H² 上升、孔的产量随坝顶上升而下降（上游坝趾更靠内、坑更小），两线交点就是各策略的坝顶上限。"))
d1.append('<h2 id="s3">1.3 准则对照（资料 1–5、7、8、10、11）</h2>')
HA = round(B["A"]["H_max"])
crit = [
    ["断面依据", "EM 2200：上游直立、下游 0.7–0.8:1；RCC 无模板 ≥0.8:1；坝顶 ≥20 ft（RCC 手册）", "小坝设计 / DS-13：上游 1.3–1.7:1，下游休止角 1.3–1.4:1；顶宽 ≥20 ft；拱度 1% H", "案例：上游 1.8–2:1，下游 1.6–2:1；心墙底宽 ≥25% 水头"],
    ["在本坝线上能否闭合", "0.8:1 面在 ≤51° 坡上闭合：全环闭合，下游趾 60–315 m", "1.4:1 面在 ≥35.5° 坡上永不闭合；第 1、4 段（30–42°）趾 500–930 m，第 4 段出宗地北界", "2.0/1.8 面 ≥28–29° 不闭合：第 3、4 段 12–15% 桩号不闭合，几乎全环不适用 → 不进方案"],
    ["稳定初筛（φ=40°、c=0、排水有效率 0.5）", f"H {H6} m FS {GS[H6][0]:.2f}、{HA} m {GS[HA][0]:.2f}、129 m {GS[129][0]:.2f}；FERC 2A（≥1.5）全过，EM 2200 usual（≥2.0）到 ~60 m，FERC 高危害（≥3.0）只到 30 m → 方案 B 的 129 m 段要靠层面粘聚力（RCC 手册 c 300 psi 级）或锚固；方案 A 的 {HA} m 段接近 2.0", "好岩 + 好基础可只做无限边坡（DS-13 2.3.5.3）；否则 Spencer FS 1.5 / 1.3", "Spencer 1.5 / 1.3；心墙开裂假定必有"],
    ["基础要求", "微–中风化微裂隙岩；基底水平或向上游倾；台阶接头应力集中；固结灌浆 30 ft；帷幕 H/3 + 50 ft", "趾板落在可灌浆基岩上；堆石区可容忍差些的基础", "心墙下坝肩 ≤0.5H:1V 整形"],
    ["斜坡基础（本地块 18–42°）", "台阶基底，坝趾开挖深 ≈1.55 H（笔记 4）；30–42° 段（第 1、4 段）台阶开挖量大", "趾板在 21–34° 坡上；下游堆石落在 30° 以上坡上时底宽 ≥ 8 H", "—"],
    ["坑与坝的关系（本框架特有）", "直立上游面：坑壁顶线离坝轴 ≥ 20 m 平台即可；坑壁 0.75:1 在坝趾附近的岩体里，需按坝基应力核", "上游 1.4:1 坡脚落在坑壁顶线外 20 m；坑壁与堆石体之间的平台是趾板与帷幕的位置，宽度要按坑壁稳定与趾板宽度核", "—"],
    ["地震（PGA 0.35 g / 2,475 年）", ">50 ft 必须动力分析；震后 FS 1.3；全层面垫层砂浆", f"超高 3–5% H（A：{0.03*HA:.1f}–{0.05*HA:.1f} m）；防御措施清单", "心墙塑性、宽反滤贯穿全高"],
    ["料（资料 6/7，本框架下全部来自孔）", "骨料加工厂在坑边；胶材山路外运；水源 14 天养护", "完好岩块：变质火山岩候选；千枚岩片状料的 1.2 折算与压实特性要试验填筑证明；泥流盖层耐久性待验", "心墙料可用风化岩；反滤要制砂"],
    ["压实与记录（资料 7/8/10/11）", "每层 Vebe / 核子密度、层面成熟度；试验段 100 ft", "层厚 0.6–1.0 m、10–20 t 振动碾、遍数由试验填筑定；回数 100% 轨迹记录；ICMV 找弱区", "心墙点测 + 回数；反滤专门验收"],
    ["优点（资料 1/3）", "可短时漫顶；坝体最小；陡坡上唯一能闭合的断面", "可检修、灌浆不在关键路径、坝体不饱和、易加高、雨季可施工", "材料适应性最广"],
    ["本框架下的主要开口", f"B：4.6 km 环上 {fnum(bB['dam_grid_Mm3']):.0f} Mm³ RCC 的热控、拌合与十余年拌合工期，余料只有 {fnum(bB['maxnet_surplus_ratio']):.2f} 倍；A：第 1、4 段 RCC 与 CFRD 的接头（应力集中、地震开裂位置）", f"A/C：坝顶被孔的产量封顶在 {CREST_A:.0f} / {CREST_C:.0f}；上游楔吃掉 3–6 Mm³ 库容；第 5、6、12 段坝趾出界要内移或改线", "不进方案"],
]
d1.append(table_html(["项", "RCC 重力（方案 A 的 1、4 段；方案 B 全环）", "面板堆石 CFRD（方案 A 其余段；方案 C 全环）", "心墙堆石 ECRD（为何不进方案）"], crit))
d1.append(f'<h2 id="s4">1.4 逐段表：方案 A（坝顶 {CREST_A:.0f}）与方案 B（坝顶 4,120）</h2>')
d1.append('<p class="small">分段沿用固定坝顶版交付物 1 的 12 段（按 H 档与 CFRD 坝趾距离档切分）。每格：H 范围；坝型；下游坝趾距轴线距离与底宽（段内最大）；平均断面法方量；下游坝趾出宗地的桩号比例 / 最大出界距离。"同坝顶若全段 CFRD / RCC"两列给出方案 A 的坝顶下换用另一种断面的几何，是每段选型的依据。地质是 1:250k 图上的"可能单元"。</p>')
d1.append(table_html(list(RROWS[0].keys()), RROWS))
d1.append('<div class="note"><b>三条贯穿各段的事实。</b>（1）第 4 段（东北坡 31–42°）在任何坝顶下堆石坝趾都出北界，所以这块地上没有一个坝顶高程能让整条坝线都用堆石而坝趾都在界内——第 4 段总要 RCC 或改线。'
          f'（2）坝顶降到 {CREST_A:.0f} 后，东界与南台地大部分地面高于坝顶，水靠天然山体挡住，坝只剩 {B["A"]["dam_len"]/1000:.1f} km；固定坝顶 4,120 时这些段坝趾 100% 出界的问题随之消失。'
          '（3）第 5、12 段（沿东界、西界）两种坝型的下游坝趾都在界外：要么坝轴内移（第 5 段 80–100 m），要么用地（第 12 段西侧 60–100 m 宽的带），要么在宗地内沿更高等高线回折闭合。</div>')
d1.append('<h2 id="s5">1.5 换一条坝线会怎样：三个环</h2>')
rings = []
for L, nm, ln, ac in (("A_parcel_ring_3700", "3,700 ft 环（本文）", "4,589 m / 139 acre", ""), ("A_parcel_ring_3800", "3,800 ft 环", "4,222 m / 92 acre", ""), ("A_parcel_ring_3900", "3,900 ft 环", "3,470 m / 56 acre", "")):
    fn = os.path.join(OUT, f"balanced_{L}.csv")
    if not os.path.exists(fn): continue
    rr = {(round(float(r["crest_ft"])), r["type"]): r for r in csv.DictReader(open(fn))}
    def best_of(ty):
        cs = [c for (c, tt), r in rr.items() if tt == ty and r["balanced"] == "True"]
        if not cs: return "无平衡坝顶"
        c = max(cs); r = rr[(c, ty)]; return f"{c} ft：净 {fnum(r['maxnet_Mm3']):.0f} Mm³ / {fnum(r['maxnet_GWh']):.0f} GWh（底 {fnum(r['maxnet_floor_ft']):.0f}，挖÷需 {fnum(r['maxnet_surplus_ratio']):.1f}）"
    rB = rr.get((4120, "rcc"))
    rings.append([nm, ln, best_of("cfrd"), best_of("combo") if any(tt == "combo" for (_, tt) in rr) else "（未定义）", f"净 {fnum(rB['maxnet_Mm3']):.0f} Mm³ / {fnum(rB['maxnet_GWh']):.0f} GWh（底 {fnum(rB['maxnet_floor_ft']):.0f}，挖÷需 {fnum(rB['maxnet_surplus_ratio']):.2f}，{'平衡' if rB['balanced']=='True' else '料不够'}）" if rB else "—"])
d1.append(table_html(["坝线", "长 / 环内面积", "全 CFRD 最高平衡坝顶与库容", "组合 III 最高平衡坝顶与库容", "全 RCC 坝顶 4,120 的库容"], rings))
d1.append('<p class="small">坝线越往上抬，堆石能平衡的坝顶越高，但环内面积缩得更快，库容反而下降；所以"打孔为主"把坝线留在 3,700 环。数据：outputs/balanced_A_parcel_ring_{3700,3800,3900}.csv（碗形坑规则）。</p>')
d1.append('<h2 id="s6">1.6 这张表没有回答的（留给交付物 3–6 与现场）</h2>')
d1.append("<ul><li>孔里的岩石够不够格：Calaveras 千枚岩片状颗粒的压实特性与 1.2 折算、变质火山岩的分布、Mehrten 盖层能否作坝料（交付物 3 与试验采场，资料 7）。盖层若只能弃，堆石上限约再降 50 ft。</li>"
          "<li>坑本身：碗形坑坑壁 0.75:1 是假设——完好岩可更陡（0.5:1 带马道）则孔能出更多料，千枚岩沿片理可能更缓；坑壁稳定、渗漏与排水、进出水口位置都在这个数字之外。</li>"
          "<li>每段坝基的真实岩性、风化深度、RQD、软弱面产状（决定 RCC 能不能落在第 1、4 段，CFRD 趾板怎么做）。</li>"
          "<li>A 的 RCC/CFRD 接头（第 1|2、3|4、4|5 段之间）的三维应力与地震开裂——DS-13 建议数值分析。</li>"
          "<li>下游坝趾出界段的用地或改线；净库容 → 电量还要下库位置与水位。</li></ul>")
d1.append('<h2 id="s7">复算</h2>')
d1.append("<pre><code>python3 scripts/terrain_model.py --build-lines --level-ft 3700                 # 坝线\npython3 scripts/balanced_design_scan.py A_parcel_ring_3700 --types=cfrd,rcc,combo # 坝顶扫描（碗形坑）\npython3 scripts/deliverables_1_2_html.py                                        # 本文件与交付物 2</code></pre>")
d1.append('<p class="small">数据：data/dem_3dep_1m_utm10.tif（USGS 3DEP 1 m，EPSG:32610）、data/parcel_014160001000_utm10.json、data/dam_line_candidates.json；中间结果 outputs/D1v2_schemes.csv、outputs/D1v2_reaches.csv。固定坝顶 4,120 的旧版见 deliverables/D1_dam_type_comparison.md。</p>')
open(os.path.join(DEL, "D1_dam_type_comparison.html"), "w").write(html_doc("交付物 1：坝型比较表（打孔为主框架）", "\n".join(d1)))

# ---- D2
fA, fB = SCHEMES[0]["floor_deep"], SCHEMES[1]["floor_deep"]
d2 = []
d2.append(f"<h1>交付物 2：整条坝线纵剖面与代表断面（打孔为主框架）——坝线 {LINE}</h1>")
d2.append('<nav><a href="#p1">2.1 纵剖面</a><a href="#p2">2.2 方案 A 断面</a><a href="#p3">2.3 方案 B 断面</a><a href="#p4">2.4 断面表</a><a href="#p5">2.5 图直接说明的事</a></nav>')
d2.append(f'<div class="key">坝线：宗地内沿 3,700 ft 等高线闭合的环，桩号从南界西端起{CW}。方案 A：组合 III（第 1、4 段 RCC，其余 CFRD），坝顶 {CREST_A:.0f} / 水位 {CREST_A-20:.0f} / 库底 {fA:.0f} ft（深挖参考点）；方案 B：全环 RCC，坝顶 4,120 / 水位 4,100 / 库底 {fB:.0f} ft。'
          '断面几何：RCC 顶宽 6.1 m、上游直立、下游 0.8:1；CFRD 顶宽 10 m、上下游 1.4:1。坑按笔记 14 修正后的规则：坑壁顶线在上游坝趾内 20 m，坑壁 0.75:1，底面平到库底；断面上画的开挖面就是这条规则给出的面。</div>')
d2.append('<h2 id="p1">2.1 纵剖面</h2>')
d2.append(fig_html(os.path.join(OUT, "D2v2_profile.png"), "图 2.1 上：轴线地面（黑）；方案 A 的坝体按坝型着色（橙 RCC、绿 CFRD），方案 B 在其上再加的 RCC 用点状阴影；两条坝顶/水位与两条库底；绿点 / 红叉是 A / B 的下游坝趾高程（越低说明坝要伸到多深的坡上）；底色是可能的地质单元。下：两个坝顶下的 H、下游与库侧地面坡角；35.5° 是 1.4:1 面的坡角。"))
ra = RR["A"]; rb = RR["B"]
d2.append(f'<ul><li>0–2,695 m（第 1–4 段，沿 3,700 线绕鼻尖）：A 的 H {ra[0]["H_min"]:.0f}–{max(r["H_max"] for r in ra[:4]):.0f} m，B 的 H 126–129 m；下游 100 m 内坡角 18–42°。A 在第 1、4 段（30–42°）用 RCC，下游趾 {-ra[0]["toe_ds_max"]:.0f}–{-ra[3]["toe_ds_max"]:.0f} m；第 2、3 段 CFRD 趾 {-ra[1]["toe_ds_max"]:.0f}–{-ra[2]["toe_ds_max"]:.0f} m。B 全 RCC 趾 213–311 m。</li>'
          f'<li>2,705–2,893 m（第 5 段，东界过渡）：地面从 3,700 升到 3,855 ft；两方案的坝趾都出东界。</li>'
          f'<li>2,903–4,339 m（第 6–11 段，东界与南台地）：地面 3,867–4,234 ft。A {"只在第 6 段前 %.0f m 有 H ≤ %.0f m 的低坝，其余" % (ra[5]["dam_len_m"], ra[5]["H_max"]) if ra[5]["dam_len_m"] > 0 else "全段"}地面高于坝顶，天然山体挡水；B 在第 6、7、8、10 段有 H 2–77 m 的坝，第 9、11 段无坝。</li>'
          f'<li>4,349–4,589 m（第 12 段，西南闭合）：沿西界顺坡而下，A 的 H 0–{ra[11]["H_max"]:.0f} m，B 0–128 m；坝趾都在界外。</li></ul>')
d2.append(f'<h2 id="p2">2.2 方案 A 的代表断面（坝顶 {CREST_A:.0f}，库底 {fA:.0f}）</h2>')
d2.append(fig_html(os.path.join(OUT, "D2v2_sections_A.png"), f"图 2.2 方案 A 七个代表断面（位置见交付物 1 图 1.1 的 S#）：黑线地面，彩色是坝体（橙 RCC / 绿 CFRD），橙红是碗形坑在该断面上的开挖（库底 {fA:.0f} ft），蓝是水；点划线是宗地界在断面上的位置（下游侧）。+ 为库侧；只画库侧的开挖。"))
d2.append(f'<h2 id="p3">2.3 方案 B 的代表断面（坝顶 4,120，库底 {fB:.0f}）</h2>')
d2.append(fig_html(os.path.join(OUT, "D2v2_sections_B.png"), f"图 2.3 方案 B 七个代表断面：全环 RCC（直立上游面、0.8:1 下游面），库底 {fB:.0f} ft。"))
d2.append('<h2 id="p4">2.4 断面表</h2>')
tcols = ["方案", "断面", "桩号 m", "段", "地面 ft", "H m", "下游坡°", "库侧坡°", "地质（可能）", "坝型", "面积 m²", "底宽 m", "上游趾 m", "下游趾 m", "下游趾高程 ft", "上游楔 m²", "断面开挖 m²", "最大挖深 m", "宗地界（下游）m"]
trows = [[rec["scheme"], f"S{rec['reach']}", f"{rec['station']:.0f}", rec["reach"], f"{rec['z_ground_ft']:.0f}", f"{rec['H_m']:.0f}", f"{rec['slope_ds_deg']:.0f}", f"{rec['slope_us_deg']:.0f}", rec["geology"],
          rec["kind"].upper() if rec["kind"] != "无坝" else "无坝", f"{rec['area_m2']:.0f}", fmt(rec["base_m"]), fmt(rec["toe_us_m"]), fmt(rec["toe_ds_m"]), fmt(rec["z_toe_ds_ft"]), f"{rec['us_wedge_m2']:.0f}", f"{rec['cut_area_m2']:.0f}", f"{rec['cut_depth_max_m']:.0f}", fmt(rec["parcel_edge_ds_m"])] for rec in TAB_A + TAB_B]
d2.append(table_html(tcols, trows))
d2.append('<p class="small">面积 = 坝面与地面之间的断面面积；上游楔 = 库侧水位以下坝面与地面之间的面积（算库容要扣）；断面开挖 = 该断面库侧地面与开挖面之间的面积；宗地界列给出下游侧宗地界离坝轴的距离（空 = 900 m 内没有出界）。数据：outputs/D2v2_sections_table.csv。</p>')
d2.append('<h2 id="p5">2.5 这些图直接说明的事</h2>')
sA = {r["reach"]: r for r in TAB_A}
d2.append(f'<ul><li><b>坝顶从 4,120 降到 {CREST_A:.0f}，坝线上一半的坝消失了。</b>A 的坝只剩 {B["A"]["dam_len"]/1000:.1f} km（第 1–5 段{"、第 6 段前 %.0f m" % ra[5]["dam_len_m"] if ra[5]["dam_len_m"] > 0 else ""}、第 12 段），东界与南台地是天然山体；B 的坝 {B["B"]["dam_len"]/1000:.1f} km 绕整个环，其中第 5–8、10、12 段的坝趾在界外。</li>'
          f'<li><b>第 1、4 段为什么是 RCC。</b>S1（{sA["1"]["slope_ds_deg"]:.0f}°）与 S4（{sA["4"]["slope_ds_deg"]:.0f}°）上，1.4:1 的堆石面几乎与地面平行，坝趾要伸到 {-RR["A_cfrd"][0]["toe_ds_max"]:.0f}–{-RR["A_cfrd"][3]["toe_ds_max"]:.0f} m 外、底宽 {RR["A_cfrd"][0]["base_max"]:.0f}–{RR["A_cfrd"][3]["base_max"]:.0f} m，方量是 RCC 的 {RR["A_cfrd"][0]["vol_ea_Mm3"]/RR["A_rcc"][0]["vol_ea_Mm3"]:.0f}–{RR["A_cfrd"][3]["vol_ea_Mm3"]/RR["A_rcc"][3]["vol_ea_Mm3"]:.0f} 倍；S4 的堆石趾还出宗地北界。RCC 在同一位置底宽 {sA["1"]["base_m"]:.0f}–{sA["4"]["base_m"]:.0f} m。</li>'
          f'<li><b>第 2、3 段为什么是 CFRD。</b>S2（{sA["2"]["slope_ds_deg"]:.0f}°）、S3（{sA["3"]["slope_ds_deg"]:.0f}°）地面离 35.5° 够远，CFRD 趾 {-ra[1]["toe_ds_max"]:.0f}–{-ra[2]["toe_ds_max"]:.0f} m 都在界内；这两段 1.7 km 是孔料的主要去处，坑壁顶线就在上游坝趾以内 20 m。</li>'
          f'<li><b>坑长什么样。</b>库底 {fA:.0f} ft 比最低坝基低 {(zmin_axis - fA*FT):.0f} m：坑壁从上游坝趾内 20 m 处以 0.75:1 下到库底，S2 断面最深 {sA["2"]["cut_depth_max_m"]:.0f} m、开挖 {sA["2"]["cut_area_m2"]/1e4:.1f} 万 m²；地面高的段落（S3、S5）坑壁占去的带更宽，这就是初版"统一退让"低估的部分。B 的库底 {fB:.0f} 更深。</li>'
          f'<li><b>上游楔。</b>CFRD 断面水位以下的坝体（S2 约 {sA["2"]["us_wedge_m2"]/1e4:.1f} 万 m²/m）合计 {B["A"]["wedge"]/1e6:.1f} Mm³ 要从库容里扣；RCC 直立面几乎没有。</li>'
          + (f'<li><b>低坝段。</b>S6（桩号 {sA["6"]["station"]:.0f}，H {sA["6"]["H_m"]:.0f} m）是 A 唯一的低坝段，题设"低段可用混凝土"在这里成立；它的下游坝趾出东界。</li>' if sA["6"]["H_m"] > 0 else f'<li><b>S6（桩号 {sA["6"]["station"]:.0f}）地面高于坝顶</b>：方案 A 在东界没有坝，坑壁顶线离宗地界 20 m。</li>') + '</ul>')
d2.append('<p class="small">复算：<code>python3 scripts/deliverables_1_2_html.py</code>；图 outputs/D2v2_profile.png、D2v2_sections_A.png、D2v2_sections_B.png；固定坝顶 4,120 的旧版见 deliverables/D2_profile_and_sections.md。</p>')
open(os.path.join(DEL, "D2_profile_and_sections.html"), "w").write(html_doc("交付物 2：纵剖面与代表断面（打孔为主框架）", "\n".join(d2)))
print("written:", [os.path.getsize(os.path.join(DEL, f)) // 1024 for f in ("D1_dam_type_comparison.html", "D2_profile_and_sections.html")], "KB; CREST_A", CREST_A, "CREST_C", CREST_C)
