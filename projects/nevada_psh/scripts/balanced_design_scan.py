#!/usr/bin/env python3
"""平衡设计扫描：以"打孔"为主的上库——所有坝料来自库内开挖，坝按开挖量来修。

对一条坝线，把坝顶高程当变量（正常水位 = 坝顶 − 20 ft），对每种断面（CFRD 1.4/1.4、RCC 直立/0.8）：
  1. 逐站放坝 → 坝体压实方（平均断面法 + 网格法）、上游坝趾 → 可挖区多边形；
  2. 库底高程 F 从坝顶以下向下扫：可挖区按 0.75:1 切坡退让；挖方 cut(F)；挖后毛库容；净库容 = 毛 − 上游楔；
  3. 需求：堆石 = 压实方 / SF_rock（1 m³ 原岩压实后约 1.2 m³ 堆石，取 1.1–1.3 的中值）；RCC 骨料 = 0.9 × RCC 方量；
  4. 找最浅的 F 使 cut(F) ≥ 需求 —— 这就是"料全部来自孔里"的平衡点；记录该点的库容、电量、坝体、孔深。
输出：outputs/balanced_<line>.csv 与 outputs/balanced_<line>.png。
坑的规则（--wall=bowl，默认）：碗形坑，坑壁顶线 = 上游坝趾多边形向内退 20 m 平台，坑壁 0.75:1 向内下降到底面 F；
  开挖面 = min(地面, max(F, 坑壁面))。旧规则（--wall=uniform）：底面平、可挖区 = 坝趾多边形统一退让 (最低坝基 − F)×0.75，
  在地面高的段落把坑壁算到了坝体下面和宗地外，挖方偏大 1.5–2.5 倍（见笔记 14 §2.0 的修正）。
"深挖参考点"（列名仍为 maxnet_*）= 净库容达到几何极限（坑壁相遇）95% 的最浅底面；"平衡点" = 挖方刚够坝料的最浅底面。

用法：python3 balanced_design_scan.py A_parcel_ring_3700 [3800 3850 ... 坝顶ft列表]
"""
from __future__ import annotations
import csv, json, os, sys
import numpy as np
import shapely
from shapely.geometry import LineString, Polygon
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from terrain_model import (Terrain, SECTIONS, FT, DATA, OUT, cross_section, dam_section, densify,
                           inward_normals, longitudinal_profile)  # noqa

ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
TYPES_ARG = [a for a in sys.argv[1:] if a.startswith("--types=")]
WALL_ARG = [a for a in sys.argv[1:] if a.startswith("--wall=")]
WALL = WALL_ARG[0].split("=", 1)[1] if WALL_ARG else "bowl"   # bowl（默认）| uniform（旧规则，见下）
WALL_M = 0.75      # 坑壁坡比 H:V（0.75:1 ≈ 53°）
BERM = 20.0        # 上游坝趾（或无坝处的坝轴/宗地界）与坑壁顶之间留的平台 m
CAP_FT = 3950.0    # Mehrten 泥流盖层下界：以上算"盖层"挖方
LINE = ARGS[0] if ARGS else "A_parcel_ring_3700"
CRESTS = [float(x) for x in ARGS[1:]] or [3800, 3850, 3900, 3950, 4000, 4050, 4120]
SF_ROCK = 1.2      # 压实堆石 m³ / 原岩 m³（松方 1/0.65 ≈ 1.54，压实后约 ×0.78）
AGG_RCC = 0.9      # 骨料原岩 m³ / RCC m³
LOWER_FT = 2400.0  # 下库若在北界河谷
TYPES = tuple(TYPES_ARG[0].split("=", 1)[1].split(",")) if TYPES_ARG else ("cfrd", "rcc")
# "combo"：陡坡段（交付物 1 的第 1、4 段，下游坡 30°+）用 RCC，其余用 CFRD——只对 3,700 环定义
COMBO_RCC_STA = {"A_parcel_ring_3700": [(0.0, 606.0), (2302.8, 2695.2)]}
def station_spec(ty, st):
    if ty != "combo": return SECTIONS[ty], ty
    for a, b in COMBO_RCC_STA.get(LINE, []):
        if a - 1e-6 <= st <= b + 1e-6: return SECTIONS["rcc"], "rcc"
    return SECTIONS["cfrd"], "cfrd"

t = Terrain(os.path.join(DATA, "dem_3dep_1m_utm10.tif"), os.path.join(DATA, "dem_3dep_1m_utm10.json"))
cand = json.load(open(os.path.join(DATA, "dam_line_candidates.json")))[LINE]
line = LineString(cand["coords"]); ring = Polygon(cand["coords"])
s, pts, zg = longitudinal_profile(t, line, 10.0); nrm = inward_normals(pts, ring, None)
secs = [cross_section(t, pts[i], nrm[i]) for i in range(len(pts))]
zmin_axis = float(np.nanmin(zg))
# 库内像元
win, m_ring = t.mask_inside(ring); z_ring = t.z[win]
X, Y = np.meshgrid(t.xc[win[1]], t.yc[win[0]])
# 网格法准备
from scipy.spatial import cKDTree
dense = densify(line, 2.0); dn = inward_normals(dense, ring, None); tree = cKDTree(dense)
st_dense = np.concatenate([[0], np.cumsum(np.hypot(np.diff(dense[:, 0]), np.diff(dense[:, 1])))])
bwin, bm = t.mask_inside(line.buffer(950.0))
BX, BY = np.meshgrid(t.xc[bwin[1]], t.yc[bwin[0]]); bz = t.z[bwin]; bsel = bm & np.isfinite(bz)
_, bidx = tree.query(np.column_stack([BX[bsel], BY[bsel]]))
bvec = np.column_stack([BX[bsel], BY[bsel]]) - dense[bidx]; boff = np.einsum("ij,ij->i", bvec, dn[bidx])

def grid_volume(spec, crest, toe_us, toe_ds, feas):
    a = np.abs(boff); b = spec["crest_w"]; mm = np.where(boff > 0, spec["m_us"], spec["m_ds"])
    with np.errstate(divide="ignore", invalid="ignore"):
        face = np.where(a <= b / 2, crest, np.where(mm > 0, crest - (a - b / 2) / np.where(mm > 0, mm, 1), -np.inf))
    thick = np.maximum(face - bz[bsel], 0.0)
    fd = np.interp(st_dense[bidx], s, feas.astype(float)) >= 0.999
    tu = np.interp(st_dense[bidx], s, np.nan_to_num(toe_us, nan=0.0)); td = np.interp(st_dense[bidx], s, np.nan_to_num(toe_ds, nan=0.0))
    return float(np.sum(thick[(boff <= tu) & (boff >= td) & fd]) * t.px ** 2)

def wall_field(pit0):
    """碗形坑的坑壁约束场：对坑壁顶线（上游坝趾多边形向内退 BERM）上每 5 m 一个点 b（地面 z_b），
    库内像元 p 的坑壁面高程 wmax(p) = max_b [z_b − d(p,b)/WALL_M]。返回 inside 掩膜与 wmax。"""
    top = pit0.buffer(-BERM)
    if top.is_empty: return {"inside": np.zeros_like(z_ring, bool), "wmax": np.full_like(z_ring, np.inf)}
    if top.geom_type != "Polygon": top = max(top.geoms, key=lambda g: g.area)
    bpts = np.array(shapely.segmentize(top.exterior, 5.0).coords); zb = t.elev(bpts[:, 0], bpts[:, 1])
    ok = np.isfinite(zb); bpts = bpts[ok]; zb = zb[ok]
    inside = shapely.contains_xy(top, X.ravel(), Y.ravel()).reshape(X.shape) & np.isfinite(z_ring)
    wmax = np.full(z_ring.shape, np.inf)
    idx = np.where(inside.ravel())[0]; px_ = X.ravel()[idx]; py_ = Y.ravel()[idx]; out = np.empty(len(idx))
    CH = 20000
    for a in range(0, len(idx), CH):
        dx = px_[a:a + CH, None] - bpts[None, :, 0]; dy = py_[a:a + CH, None] - bpts[None, :, 1]
        out[a:a + CH] = np.max(zb[None, :] - np.hypot(dx, dy) / WALL_M, axis=1)
    wmax.ravel()[idx] = out
    return {"inside": inside, "wmax": wmax}


rows = []
for crest_ft in CRESTS:
    crest = crest_ft * FT; nwl = crest - 20 * FT
    for ty in TYPES:
        if ty == "combo" and LINE not in COMBO_RCC_STA: continue
        A = []; W = []; TU = []; TD = []; KIND = []
        for i in range(len(pts)):
            spec, kind = station_spec(ty, s[i]); KIND.append(kind)
            d = dam_section(*secs[i], crest, spec, nwl)
            A.append(d["area"]); W.append(d["area_us_below_nwl"]); TU.append(d["toe_us"] if d["toe_us"] is not None else np.nan); TD.append(d["toe_ds"] if d["toe_ds"] is not None else np.nan)
        A = np.array(A); W = np.array(W); TU = np.array(TU); TD = np.array(TD); feas = np.isfinite(A); KIND = np.array(KIND)
        V_ea = float(sum(0.5 * (A[i] + A[i + 1]) * 10 for i in range(len(A) - 1) if feas[i] and feas[i + 1]))
        wedge = float(sum(0.5 * (W[i] + W[i + 1]) * 10 for i in range(len(A) - 1) if feas[i] and feas[i + 1]))
        V_grid_by = {k: grid_volume(SECTIONS[k], crest, TU, TD, feas & (KIND == k)) for k in set(KIND)}
        V_grid = float(sum(V_grid_by.values()))
        need = float(sum(v / SF_ROCK if k != "rcc" else v * AGG_RCC for k, v in V_grid_by.items()))
        # 可挖区
        pit_pts = [pts[i] + max(np.nan_to_num(TU[i], nan=0.0), 0.0) * nrm[i] for i in range(len(pts))]
        pit0 = Polygon(pit_pts).buffer(0)
        if pit0.geom_type != "Polygon": pit0 = max(pit0.geoms, key=lambda g: g.area)
        best = None; curve = []
        if WALL == "bowl":
            # 碗形坑：坑壁从（上游坝趾 − 平台）处的天然地面以 0.75:1 向内下降，底面平到 F；
            # 任一点的开挖面 = min(地面, max(F, max_b[z_b − d(p,b)/0.75]))，b 取坑壁顶线上的点。
            # 这样坑壁在地面高的段落自动占更宽的带，不会伸到坝体下面或宗地外。
            wall = wall_field(pit0)
        for F_ft in np.arange(crest_ft - 50, 2950, -50):
            Fm = F_ft * FT
            if WALL == "bowl":
                z_new = np.where(wall["inside"], np.minimum(z_ring, np.maximum(Fm, wall["wmax"])), z_ring)
            else:
                sb = max(0.0, zmin_axis - Fm) * WALL_M
                pit = pit0.buffer(-sb)
                if pit.is_empty: break
                mp = shapely.contains_xy(pit, X.ravel(), Y.ravel()).reshape(X.shape)
                z_new = np.where(mp, np.minimum(z_ring, Fm), z_ring)
            dz = np.maximum(z_ring - z_new, 0.0); cut = float(np.nansum(dz))
            cap = float(np.nansum(np.maximum(z_ring - np.maximum(z_new, CAP_FT * FT), 0.0)))
            gross = float(np.nansum(np.where(m_ring, np.maximum(nwl - z_new, 0.0), 0.0)))
            net = gross - wedge; area = float((dz > 0.01).sum() / 4046.86)
            curve.append((F_ft, cut, gross, net, area, cap))
            if best is None and cut >= need:
                best = curve[-1]
        head = (nwl / FT - LOWER_FT) * FT
        gwh = lambda v: v * 1000 * 9.81 * head * 0.85 / 3.6e12
        max_cut = max(c[1] for c in curve) if curve else 0.0
        if WALL == "bowl":
            # 深挖参考点：净库容达到几何极限（坑壁在中间相遇、底面消失）的 95% 的最浅底面
            net_inf = curve[-1][3] if curve else 0.0
            top = next((c for c in curve if c[3] >= 0.95 * net_inf), curve[-1]) if curve else None
        else:
            top = max(curve, key=lambda c: c[3]) if curve else None
        rec = {"line": LINE, "crest_ft": crest_ft, "nwl_ft": crest_ft - 20, "type": ty, "H_max_m": float(np.nanmax(crest - zg)),
               "dam_endarea_Mm3": V_ea / 1e6, "dam_grid_Mm3": V_grid / 1e6, "grid_cfrd_Mm3": V_grid_by.get("cfrd", 0.0) / 1e6, "grid_rcc_Mm3": V_grid_by.get("rcc", 0.0) / 1e6, "infeasible_m": float((~feas).sum() * 10),
               "need_bank_Mm3": need / 1e6, "wedge_Mm3": wedge / 1e6,
               "net_nocut_Mm3": (curve[0][2] - wedge) / 1e6 if curve else np.nan,
               "max_cut_Mm3": max_cut / 1e6,
               "balanced": best is not None,
               "floor_ft": best[0] if best else "", "cut_Mm3": best[1] / 1e6 if best else "", "cap_Mm3": best[5] / 1e6 if best else "", "pit_acre": best[4] if best else "",
               "net_Mm3": best[3] / 1e6 if best else "", "net_GWh": gwh(best[3]) if best else "",
               "maxnet_floor_ft": top[0] if top else "", "maxnet_cut_Mm3": top[1] / 1e6 if top else "", "maxnet_cap_Mm3": top[5] / 1e6 if top else "", "maxnet_pit_acre": top[4] if top else "",
               "maxnet_Mm3": top[3] / 1e6 if top else "", "maxnet_GWh": gwh(top[3]) if top else "",
               "maxnet_surplus_ratio": (top[1] / need) if (top and need) else "",
               "netinf_Mm3": (curve[-1][3] / 1e6) if curve else "", "wall": WALL, "berm_m": BERM, "wall_m": WALL_M,
               "head_m": head}
        rec["curve"] = ";".join(f"{c[0]:.0f}:{c[1]/1e6:.2f}:{c[3]/1e6:.2f}" for c in curve)  # 底面:挖方:净库容
        rows.append(rec)
        fm = lambda v: f"{v:.1f}" if isinstance(v, float) else str(v)
        print(f"crest {crest_ft:.0f} {ty:4s} Hmax {rec['H_max_m']:5.0f} dam {rec['dam_grid_Mm3']:6.1f}/{rec['dam_endarea_Mm3']:6.1f} need {rec['need_bank_Mm3']:6.1f} "
              f"maxcut {rec['max_cut_Mm3']:6.1f} {'BAL' if best else '---'} minfloor {rec['floor_ft']} net {fm(rec['net_Mm3'])} | maxnet floor {rec['maxnet_floor_ft']} "
              f"cut {fm(rec['maxnet_cut_Mm3'])} net {fm(rec['maxnet_Mm3'])} GWh {fm(rec['maxnet_GWh'])} surplus x{fm(rec['maxnet_surplus_ratio'])} (no-cut net {fm(rec['net_nocut_Mm3'])})")
cols = list(rows[0].keys())
main = os.path.join(OUT, f"balanced_{LINE}.csv")
if TYPES_ARG and os.path.exists(main):  # 只算了部分坝型：并入已有结果（同坝型同坝顶的旧行被替换）
    old = [r for r in csv.DictReader(open(main)) if not ((r["type"] in TYPES) and (float(r["crest_ft"]) in CRESTS))]
    for r in old:
        for c in cols: r.setdefault(c, "")
    rows = sorted(old + [{c: r.get(c, "") for c in cols} for r in rows], key=lambda r: (float(r["crest_ft"]), str(r["type"])))
with open(main, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
from balanced_plot import plot_balanced  # noqa
plot_balanced(LINE)
