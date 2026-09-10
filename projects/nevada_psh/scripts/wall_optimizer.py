#!/usr/bin/env python3
""""墙段自由组合"的最优库区轮廓：给定水位，每个网格单元划入库区的收益 = 蓄水量；库区边界穿过低于坝顶的地面要付坝的填方，
穿过高于坝顶的天然地面免费。对每个"填方价格" λ 解最小割（最大流），得到最优轮廓；λ 从大到小扫描，解是嵌套的，
等价于"每次加进增量库容÷增量填方最高的墙段"，得到累计填方–累计库容曲线与拐点。入围方案在 1 m 地形上逐站放坝精算，
并按笔记 14 的碗形坑规则算"围住后再挖"的库容与挖方。

阶段 1：4 m 网格、8 邻域几何割权重（Boykov–Kolmogorov），坝成本 = 断面积(H, 局部坡度) × 边长，H = 坝顶 − 边上地面。
阶段 2：outputs/wall_curves.csv（各水位、各 λ 的轮廓统计）、wall_eval.csv（精算）、wall_curves.png、wall_plan_*.png、wall_growth.png。
用法：python3 scripts/wall_optimizer.py [--crests=3900,3950,4000,4050] [--cost=rcc|cfrd] [--fast]"""
from __future__ import annotations
import csv, json, math, os, sys, time
import numpy as np, shapely, contourpy
from shapely.geometry import LineString, Point, Polygon, MultiPolygon
from shapely.ops import unary_union
from scipy import ndimage
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_flow, breadth_first_order
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
from terrain_model import Terrain, SECTIONS, FT, ACRE, DATA, OUT, Bowl, cross_section, dam_section, densify, inward_normals, load_parcel  # noqa
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
matplotlib.rcParams["font.family"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]; matplotlib.rcParams["axes.unicode_minus"] = False
import logging; logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

ARG = {a.split("=")[0]: a.split("=")[1] for a in sys.argv[1:] if "=" in a}
CRESTS = [float(x) for x in ARG.get("--crests", "3900,3950,4000,4050").split(",")]
COST_KIND = ARG.get("--cost", "rcc")
FAST = "--fast" in sys.argv
LAMBDAS = [40, 20, 10, 7, 5, 4, 3.5, 3, 2.7, 2.4, 2.1, 1.85, 1.6, 1.4, 1.2, 1.0, 0.85, 0.7, 0.55, 0.4] if not FAST else [20, 8, 4, 2.7, 2, 1.4, 1.0, 0.6]
HMAX_M = float(ARG.get("--hmax", "130"))   # 墙段最大坝高（m）：更低的地面不划入库区
K = 4; FREEBOARD_FT = 20.0; LOWER_FT = 2400.0; ETA = 0.85; EPS = 1
BERM = 20.0; WALL_M = 0.75; SF = 1.2; AGG = 0.9
B1 = {"fill": 29.7, "need": 26.7, "net_nocut": 26.5, "net": 61.4, "nwl": 4030.0, "dam_len_km": 3.9, "cut": 36.3, "gwh": 70.6}
GULLY = {"fill_rcc": 2.1, "fill_cfrd": 6.1, "net": 6.4, "nwl": 3980.0, "len": 838}
gwh = lambda v_m3, nwl_ft: v_m3 * 1000 * 9.81 * (nwl_ft - LOWER_FT) * FT * ETA / 3.6e12

t = Terrain(os.path.join(DATA, "dem_3dep_1m_utm10.tif"), os.path.join(DATA, "dem_3dep_1m_utm10.json"))
parcel = load_parcel(); pb = parcel.bounds
# ------------------------------------------------------------- 4 m 网格（山鼻区域）
ny4, nx4 = t.ny // K, t.nx // K
Z4 = t.z[:ny4 * K, :nx4 * K].reshape(ny4, K, nx4, K).mean(axis=(1, 3))
x4 = t.xmin + K * (0.5 + np.arange(nx4)); y4 = t.ymax - K * (0.5 + np.arange(ny4))
X4, Y4 = np.meshgrid(x4, y4); INP4 = shapely.contains_xy(parcel, X4.ravel(), Y4.ravel()).reshape(Z4.shape)
sel = INP4 & (Z4 / FT >= 3500)
rr, cc = np.where(sel); r0, r1 = max(rr.min() - 10, 0), min(rr.max() + 11, ny4); c0, c1 = max(cc.min() - 10, 0), min(cc.max() + 11, nx4)
z = Z4[r0:r1, c0:c1]; inp = INP4[r0:r1, c0:c1]; xg = x4[c0:c1]; yg = y4[r0:r1]; NY, NX = z.shape; N = NY * NX
gy, gx = np.gradient(np.nan_to_num(z, nan=np.nanmean(z)), K); tanb = np.clip(np.hypot(gx, gy), 0, 1.5)
print(f"网格 {NY}×{NX} = {N} 节点，{K} m")
idx = np.arange(N).reshape(NY, NX)

def section_area(H, tb, kind):
    sp = SECTIONS[kind]; m_ds, m_us, b = sp["m_ds"], sp["m_us"], sp["crest_w"]
    f_ds = 1.0 / np.maximum(1.0 - m_ds * tb, 0.3)
    return np.where(H > 0, 0.5 * H ** 2 * (m_ds * f_ds + m_us) + b * H, 0.0)

def build_edges(crest_m, kind):
    """8 邻域几何割：返回 (u, v, A_len) —— 每条无向边的坝成本(m³, λ=1) 与其长度贡献。"""
    U = []; V = []; C = []; L = []
    for (dr, dc, lf) in ((0, 1, K * math.pi / 8), (1, 0, K * math.pi / 8), (1, 1, K * math.pi / (8 * math.sqrt(2))), (1, -1, K * math.pi / (8 * math.sqrt(2)))):
        ra, rb = (0, NY - dr) if dr >= 0 else (-dr, NY); ca, cb = (max(0, -dc), NX - max(0, dc))
        a = idx[ra:rb, ca:cb]; b = idx[ra + dr:rb + dr, ca + dc:cb + dc]
        za = z[ra:rb, ca:cb]; zb = z[ra + dr:rb + dr, ca + dc:cb + dc]
        ze = 0.5 * (za + zb); H = crest_m - ze; tb = 0.5 * (tanb[ra:rb, ca:cb] + tanb[ra + dr:rb + dr, ca + dc:cb + dc])
        A = section_area(np.nan_to_num(H, nan=0.0), tb, kind)
        U.append(a.ravel()); V.append(b.ravel()); C.append((A * lf).ravel()); L.append(np.full(a.size, lf))
    return np.concatenate(U), np.concatenate(V), np.concatenate(C), np.concatenate(L)

def solve(crest_m, nwl_m, kind, lam, edges):
    U, V, C, Lf = edges
    allowed = inp & np.isfinite(z) & (z >= crest_m - HMAX_M)
    benefit = np.where(allowed, np.maximum(nwl_m - np.nan_to_num(z, nan=1e9), 0) * K * K, 0.0)
    benefit = np.where(allowed & (benefit <= 0), EPS, benefit)          # 高于水位的自由高地：象征性收益，方便后面挖
    src, snk = N, N + 1
    w = np.round(lam * C).astype(np.int64)
    rows = np.concatenate([np.full(N, src), U, V, np.where(~allowed)[0]])
    cols = np.concatenate([np.arange(N), V, U, np.full((~allowed).sum(), snk)])
    caps = np.concatenate([np.round(benefit.ravel()).astype(np.int64), w, w, np.full((~allowed).sum(), 2_000_000_000)])
    keep = caps > 0
    G = coo_matrix((caps[keep].astype(np.int32), (rows[keep], cols[keep])), shape=(N + 2, N + 2)).tocsr(); G.sum_duplicates()
    res = maximum_flow(G, src, snk, method="dinic")
    Rm = (G - res.flow.tocsr()).tocsr(); Rm.data = (Rm.data > 0).astype(np.int32); Rm.eliminate_zeros()
    reach = breadth_first_order(Rm, src, directed=True, return_predecessors=False)
    inR = np.zeros(N + 2, bool); inR[reach] = True; inR = inR[:N].reshape(NY, NX) & allowed
    # 统计
    cut = inR.ravel()[U] != inR.ravel()[V]
    fill_est = float(C[cut].sum()); dam_len = float(Lf[cut & (C > 0)].sum()); free_len = float(Lf[cut & (C == 0)].sum())
    stor = float(np.where(inR, np.maximum(nwl_m - np.nan_to_num(z, nan=1e9), 0) * K * K, 0).sum())
    return {"lam": lam, "mask": inR, "fill_est_Mm3": fill_est / 1e6, "dam_len_m": dam_len, "free_len_m": free_len, "storage_est_Mm3": stor / 1e6,
            "area_acre": float(inR.sum() * K * K / ACRE), "flow": res.flow_value}

PLAN = [(float(a.split(":")[0]), float(a.split(":")[1])) for a in ARG["--plan"].split(",")] if "--plan" in ARG else None
# ------------------------------------------------------------- 阶段 1：λ 扫描
curves = {}
for cf in ([] if PLAN else CRESTS):
    cm = cf * FT; nwl = cm - FREEBOARD_FT * FT; edges = build_edges(cm, COST_KIND); out = []
    for lam in LAMBDAS:
        t0 = time.time(); r = solve(cm, nwl, COST_KIND, lam, edges); out.append(r)
        print(f"坝顶 {cf:.0f} λ {lam:5.2f}: 填 {r['fill_est_Mm3']:6.2f} Mm³  库 {r['storage_est_Mm3']:6.2f} Mm³  坝长 {r['dam_len_m']:6.0f} m  自由岸 {r['free_len_m']:6.0f} m  面积 {r['area_acre']:5.0f} ac  ({time.time()-t0:.1f} s)")
    curves[cf] = out
if not PLAN:
  with open(os.path.join(OUT, f"wall_curves_{COST_KIND}.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["crest_ft", "nwl_ft", "lambda", "fill_est_Mm3", "storage_est_Mm3", "dam_len_m", "free_len_m", "area_acre", "marg_eff"])
    for cf, out in curves.items():
        prev = None
        for r in out:
            me = ((r["storage_est_Mm3"] - prev["storage_est_Mm3"]) / (r["fill_est_Mm3"] - prev["fill_est_Mm3"])) if prev and r["fill_est_Mm3"] > prev["fill_est_Mm3"] + 1e-6 else ""
            w.writerow([cf, cf - FREEBOARD_FT, r["lam"], f"{r['fill_est_Mm3']:.3f}", f"{r['storage_est_Mm3']:.3f}", f"{r['dam_len_m']:.0f}", f"{r['free_len_m']:.0f}", f"{r['area_acre']:.1f}", f"{me:.2f}" if me != "" else ""]); prev = r

def marg(out):
    """相邻 λ 解之间的边际效率 = Δ库容 / Δ填方。"""
    m = [np.nan]
    for a, b in zip(out[:-1], out[1:]):
        df = b["fill_est_Mm3"] - a["fill_est_Mm3"]; m.append((b["storage_est_Mm3"] - a["storage_est_Mm3"]) / df if df > 0.05 else np.nan)
    return m

def knee(out, thr=1.0):
    """停止点：最后一个"边际效率 ≥ thr"的解（再加墙，每 m³ 填方换到的库容 < thr m³；thr=1 时不如去挖：1 m³ 挖方 ≈ 1 m³ 库容还出料）。"""
    m = marg(out); best = None
    for r, e in zip(out, m):
        if r["fill_est_Mm3"] > 0.3 and (np.isnan(e) or e >= thr): best = r
        elif r["fill_est_Mm3"] > 0.3 and e < thr: break
    return best or next((r for r in out if r["fill_est_Mm3"] > 0.3), out[-1])

# ------------------------------------------------------------- 阶段 2：精算
def polygons_of(mask):
    gen = contourpy.contour_generator(x=xg, y=yg, z=mask.astype(float), name="serial")
    polys = [Polygon(l).buffer(0) for l in gen.lines(0.5) if len(l) > 3]
    polys = [p for p in polys if p.area > 0.5 * ACRE]
    # 去掉岛（被包含的多边形当作洞）
    outer = [p for p in polys if not any(q is not p and q.contains(p) for q in polys)]
    return outer

def evaluate(r, cf, label):
    cm = cf * FT; nwl = cm - FREEBOARD_FT * FT
    polys = polygons_of(r["mask"]); segs = []; tot = {k: {"V": 0.0, "wedge": 0.0, "len": 0.0, "infeas": 0.0, "toe_out": 0.0} for k in ("rcc", "cfrd")}
    free_len = 0.0; dam_len = 0.0; Hmax = 0.0; toe_polys = {}; walls = []
    for poly in polys:
        pts = densify(LineString(poly.exterior.coords), 10.0); nrm = inward_normals(pts, poly, None); zg = t.elev(pts[:, 0], pts[:, 1]); H = cm - zg
        n = len(pts); run = []
        for i in range(n):
            if H[i] > 0: run.append(i)
            if (H[i] <= 0 or i == n - 1) and run:
                seg = {"poly": poly, "idx": run, "len": len(run) * 10.0, "Hmax": float(H[run].max())}
                for k in ("rcc", "cfrd"):
                    A = []; W = []; TU = []; TD = []; out_n = 0
                    for j in run:
                        off, zs = cross_section(t, pts[j], nrm[j]); d = dam_section(off, zs, cm, SECTIONS[k], nwl)
                        A.append(d["area"]); W.append(d["area_us_below_nwl"] if np.isfinite(d["area_us_below_nwl"]) else 0)
                        TU.append(d["toe_us"] if d["toe_us"] is not None else np.nan); TD.append(d["toe_ds"] if d["toe_ds"] is not None else np.nan)
                        if d["toe_ds"] is not None and not parcel.contains(Point(pts[j] + d["toe_ds"] * nrm[j])): out_n += 1
                    A = np.array(A); ok = np.isfinite(A)
                    V = float(np.nansum([0.5 * (A[q] + A[q + 1]) * 10 for q in range(len(A) - 1) if ok[q] and ok[q + 1]]))
                    Wd = float(np.nansum([0.5 * (W[q] + W[q + 1]) * 10 for q in range(len(A) - 1) if ok[q] and ok[q + 1]]))
                    seg[k] = {"V": V, "wedge": Wd, "infeas_m": float((~ok).sum() * 10), "toe_us_max": float(np.nanmax(TU)) if np.isfinite(TU).any() else 0.0, "toe_out_m": out_n * 10.0}
                    tot[k]["V"] += V; tot[k]["wedge"] += Wd; tot[k]["infeas"] += seg[k]["infeas_m"]; tot[k]["toe_out"] += seg[k]["toe_out_m"]
                segs.append(seg); dam_len += seg["len"]; Hmax = max(Hmax, seg["Hmax"])
                walls.append(LineString(pts[run]) if len(run) > 1 else Point(pts[run[0]]).buffer(1))
                run = []
        free_len += float((H <= 0).sum() * 10)
        # 坑顶线：库区多边形减去坝趾带（按坝型分别）
        for k in ("rcc", "cfrd"):
            band = unary_union([LineString(pts[s["idx"]]).buffer(s[k]["toe_us_max"] + BERM) for s in segs if s["poly"] is poly and len(s["idx"]) > 1]) if any(s["poly"] is poly for s in segs) else None
            tp = poly.difference(band) if band is not None else poly
            if not tp.is_empty:
                if tp.geom_type != "Polygon": tp = max(tp.geoms, key=lambda g: g.area)
                toe_polys.setdefault(k, []).append(tp)
    # 1 m 洪水填充复核（墙 = 坝段栅格化）
    C0 = max(0, int(min(p.bounds[0] for p in polys) - 150 - t.xmin)); C1 = min(t.nx, int(max(p.bounds[2] for p in polys) + 150 - t.xmin))
    R0 = max(0, int(t.ymax - max(p.bounds[3] for p in polys) - 150)); R1 = min(t.ny, int(t.ymax - min(p.bounds[1] for p in polys) + 150))
    ZW = t.z[R0:R1, C0:C1]; XW, YW = np.meshgrid(t.xc[C0:C1], t.yc[R0:R1])
    wall = np.zeros(ZW.shape, bool)
    for wl in walls:
        for p in densify(wl, 1.0) if wl.geom_type == "LineString" else [np.array(wl.centroid.coords[0])]:
            cc_ = int(p[0] - t.xmin) - C0; rr_ = int(t.ymax - p[1]) - R0; wall[max(rr_ - 3, 0):rr_ + 4, max(cc_ - 3, 0):cc_ + 4] = True
    wet = (ZW < nwl) & ~wall & np.isfinite(ZW); lab, _ = ndimage.label(wet)
    comp = np.zeros(ZW.shape, bool); spill = False
    for poly in polys:
        m = shapely.contains_xy(poly, XW.ravel(), YW.ravel()).reshape(ZW.shape) & (lab > 0)
        for lb in np.unique(lab[m]):
            if lb == 0: continue
            cmp_ = lab == lb; comp |= cmp_
            if cmp_[0, :].any() or cmp_[-1, :].any() or cmp_[:, 0].any() or cmp_[:, -1].any(): spill = True
    inpoly = np.zeros(ZW.shape, bool)
    for poly in polys: inpoly |= shapely.contains_xy(poly, XW.ravel(), YW.ravel()).reshape(ZW.shape)
    gross = float(np.nansum(np.where(comp, nwl - ZW, 0))); leak_acre = float((comp & ~inpoly).sum() / ACRE)
    rec = {"label": label, "crest_ft": cf, "nwl_ft": cf - FREEBOARD_FT, "lambda": r["lam"], "n_regions": len(polys), "n_segments": len(segs), "dam_len_m": dam_len, "free_len_m": free_len, "H_max_m": Hmax,
           "area_acre": sum(p.area for p in polys) / ACRE, "fill_est_Mm3": r["fill_est_Mm3"], "storage_est_Mm3": r["storage_est_Mm3"],
           "fill_rcc_Mm3": tot["rcc"]["V"] / 1e6, "fill_cfrd_Mm3": tot["cfrd"]["V"] / 1e6, "need_rcc_Mm3": tot["rcc"]["V"] * AGG / 1e6, "need_cfrd_Mm3": tot["cfrd"]["V"] / SF / 1e6,
           "cfrd_infeasible_m": tot["cfrd"]["infeas"], "rcc_toe_out_m": tot["rcc"]["toe_out"], "cfrd_toe_out_m": tot["cfrd"]["toe_out"],
           "gross_Mm3": gross / 1e6, "net_rcc_Mm3": (gross - tot["rcc"]["wedge"]) / 1e6, "net_cfrd_Mm3": (gross - tot["cfrd"]["wedge"]) / 1e6, "spill": spill, "leak_acre": leak_acre,
           "net_rcc_GWh": gwh(gross - tot["rcc"]["wedge"], cf - FREEBOARD_FT), "head_m": (cf - FREEBOARD_FT - LOWER_FT) * FT,
           "segments": [{"len": s["len"], "Hmax": s["Hmax"], "fill_rcc": s["rcc"]["V"] / 1e6, "fill_cfrd": s["cfrd"]["V"] / 1e6, "E": float(pts_center(s)[0]), "N": float(pts_center(s)[1])} for s in segs]}
    # 再挖：碗形坑（按坝型分别：坑顶线离 RCC 直立面 24 m，离 CFRD 上游趾 20 m）
    if toe_polys and not spill:
        for k in ("rcc", "cfrd"):
            if k not in toe_polys: continue
            top = unary_union(toe_polys[k])
            if top.geom_type != "Polygon": top = max(top.geoms, key=lambda g: g.area)
            bowl = Bowl(t, top, BERM, WALL_M); inside, wmax = bowl.field(XW, YW, ZW); curve = []
            for F_ft in np.arange(cf - FREEBOARD_FT - 50, cf - FREEBOARD_FT - 750, -50):
                Fm = F_ft * FT; z_new = np.where(inside, np.minimum(ZW, np.maximum(Fm, wmax)), ZW)
                cut = float(np.nansum(np.maximum(ZW - z_new, 0))); g2 = float(np.nansum(np.where(comp | inside, np.maximum(nwl - z_new, 0), 0)))
                curve.append((F_ft, cut, g2))
            ginf = curve[-1][2]; deep = next((c for c in curve if c[2] >= 0.95 * ginf), curve[-1])
            need = tot[k]["V"] * AGG if k == "rcc" else tot[k]["V"] / SF; bal = next((c for c in curve if c[1] >= need), None)
            rec.update({f"pit_top_acre_{k}": bowl.top.area / ACRE if not bowl.top.is_empty else 0, f"deep_floor_ft_{k}": deep[0], f"deep_cut_Mm3_{k}": deep[1] / 1e6,
                        f"deep_net_Mm3_{k}": (deep[2] - tot[k]["wedge"]) / 1e6, f"deep_net_GWh_{k}": gwh(deep[2] - tot[k]["wedge"], cf - FREEBOARD_FT),
                        f"deep_cut_over_need_{k}": deep[1] / need if need else np.nan,
                        f"bal_floor_ft_{k}": bal[0] if bal else "", f"bal_cut_Mm3_{k}": bal[1] / 1e6 if bal else "", f"bal_net_Mm3_{k}": (bal[2] - tot[k]["wedge"]) / 1e6 if bal else "",
                        f"curve_{k}": ";".join(f"{c[0]:.0f}:{c[1]/1e6:.2f}:{c[2]/1e6:.2f}" for c in curve)})
            if k == "rcc": rec["_pit"] = (inside, wmax, deep[0], (R0, R1, C0, C1))
    rec["_polys"] = polys; rec["_walls"] = walls; rec["_comp"] = (comp, (R0, R1, C0, C1))
    return rec

def pts_center(s):
    pts = densify(LineString(s["poly"].exterior.coords), 10.0); return pts[s["idx"]].mean(axis=0)

EVAL = []
if PLAN:
    for cf, lam in PLAN:
        cm = cf * FT; nwl = cm - FREEBOARD_FT * FT; edges = build_edges(cm, COST_KIND); r = solve(cm, nwl, COST_KIND, lam, edges)
        rec = evaluate(r, cf, f"λ={lam:g}"); EVAL.append(rec)
        print(f"PLAN 坝顶 {cf:.0f} λ {lam}: 段 {rec['n_segments']} 坝长 {rec['dam_len_m']:.0f} 自由岸 {rec['free_len_m']:.0f} Hmax {rec['H_max_m']:.0f} 填 RCC {rec['fill_rcc_Mm3']:.1f} CFRD {rec['fill_cfrd_Mm3']:.1f} 净库 {rec['net_rcc_Mm3']:.1f} 再挖 {rec.get('deep_net_Mm3_rcc', float('nan')):.1f}（挖 {rec.get('deep_cut_Mm3_rcc', float('nan')):.1f}）")
    with open(os.path.join(OUT, "wall_plan_segments.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["crest_ft", "lambda", "seg", "E", "N", "len_m", "H_max_m", "fill_rcc_Mm3", "fill_cfrd_Mm3"])
        for rec in EVAL:
            for i, s in enumerate(sorted(rec["segments"], key=lambda s: -s["len"]), 1): w.writerow([rec["crest_ft"], rec["lambda"], i, f"{s['E']:.0f}", f"{s['N']:.0f}", f"{s['len']:.0f}", f"{s['Hmax']:.0f}", f"{s['fill_rcc']:.2f}", f"{s['fill_cfrd']:.2f}"])
    cols = sorted({k for r in EVAL for k in r if not k.startswith("_") and k != "segments" and not k.startswith("curve")}) + ["curve_rcc", "curve_cfrd"]
    with open(os.path.join(OUT, "wall_plan_eval.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols + ["segments"]); w.writeheader()
        for r in EVAL: w.writerow({**{k: r.get(k, "") for k in cols}, "segments": json.dumps(r["segments"])})
    sl = t.slope_deg(); w1, _ = t.mask_inside(parcel.buffer(60)); ext = [t.xc[w1[1]][0], t.xc[w1[1]][-1], t.yc[w1[0]][-1], t.yc[w1[0]][0]]
    def plan2(ax, rec, title, number=False):
        ax.imshow(sl[w1], cmap="Greys", extent=ext, vmin=0, vmax=55, alpha=0.65)
        ax.contour(t.xc, t.yc, t.z / FT, levels=np.arange(2600, 4400, 100), colors="#a08000", linewidths=0.3)
        px_, py_ = parcel.exterior.xy; ax.plot(px_, py_, "k--", lw=1)
        comp, (R0, R1, C0, C1) = rec["_comp"]; e2 = [t.xc[C0], t.xc[C1 - 1], t.yc[R1 - 1], t.yc[R0]]
        if "_pit" in rec:
            inside, wmax, F_ft, _ = rec["_pit"]; ZW = t.z[R0:R1, C0:C1]; z_new = np.where(inside, np.minimum(ZW, np.maximum(F_ft * FT, wmax)), ZW); dz = np.maximum(ZW - z_new, 0)
            ax.imshow(np.ma.masked_where(dz <= 0.01, dz), cmap="YlOrRd", alpha=0.8, extent=e2, vmin=0, vmax=max(60, float(np.nanmax(dz))))
        ax.imshow(np.ma.masked_where(~comp, np.ones_like(comp, float)), cmap=matplotlib.colors.ListedColormap(["#2c7fb8"]), alpha=0.4, extent=e2)
        for poly in rec["_polys"]: xx, yy = poly.exterior.xy; ax.plot(xx, yy, "-", color="#1b5e20", lw=1.3)
        for wl in rec["_walls"]:
            if wl.geom_type == "LineString": xx, yy = wl.xy; ax.plot(xx, yy, "-", color="#e65100", lw=3.2)
        for i, s in enumerate(sorted(rec["segments"], key=lambda s: -s["len"]), 1):
            if s["len"] >= 40: ax.text(s["E"], s["N"], (f"W{i}: " if number else "") + f"{s['len']:.0f} m / H {s['Hmax']:.0f} / {s['fill_rcc']:.1f} Mm³", fontsize=6.5, color="#bf360c", bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.5))
        ax.set_xlim(685900, ext[1]); ax.set_ylim(4355800, 4357100); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=8.5)
    fig, axs = plt.subplots(1, len(EVAL), figsize=(6.2 * len(EVAL), 7.2)); axs = np.atleast_1d(axs)
    for ax, rec in zip(axs, EVAL):
        plan2(ax, rec, f"坝顶 {rec['crest_ft']:.0f}（水位 {rec['nwl_ft']:.0f}）λ={rec['lambda']:g}：{rec['n_segments']} 段墙 {rec['dam_len_m']:.0f} m，自由岸 {rec['free_len_m']:.0f} m\n填 RCC {rec['fill_rcc_Mm3']:.1f} / CFRD {rec['fill_cfrd_Mm3']:.1f} Mm³；天然净库容 {rec['net_rcc_Mm3']:.1f} Mm³\n再挖到 {rec.get('deep_floor_ft_rcc', float('nan')):.0f} ft（挖 {rec.get('deep_cut_Mm3_rcc', float('nan')):.0f} Mm³）→ {rec.get('deep_net_Mm3_rcc', float('nan')):.0f} Mm³")
    axs[0].legend(handles=[plt.Line2D([], [], color="#e65100", lw=3, label="墙段（坝）：长 / 最大坝高 m / RCC 填方"), plt.Line2D([], [], color="#1b5e20", lw=1.3, label="库区边界，其余为高于坝顶的天然岸"), Patch(facecolor="#2c7fb8", alpha=0.45, label="天然蓄水面"), Patch(facecolor="#fd8d3c", label="再挖深度（碗形坑，深挖参考点）"), plt.Line2D([], [], color="k", ls="--", label="宗地界")], fontsize=6.5, loc="lower left")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "wall_plan_selected.png"), dpi=110); plt.close(fig)
    rec = EVAL[0]; fig, ax = plt.subplots(figsize=(9, 9.5)); plan2(ax, rec, f"推荐方案：坝顶 {rec['crest_ft']:.0f} / 水位 {rec['nwl_ft']:.0f} ft，λ={rec['lambda']:g}；墙段编号对应分段表", number=True)
    ax.legend(handles=[plt.Line2D([], [], color="#e65100", lw=3, label="墙段（坝）"), plt.Line2D([], [], color="#1b5e20", lw=1.3, label="库区边界，其余为高于坝顶的天然岸"), Patch(facecolor="#2c7fb8", alpha=0.45, label="天然蓄水面"), Patch(facecolor="#fd8d3c", label="再挖深度"), plt.Line2D([], [], color="k", ls="--", label="宗地界")], fontsize=7, loc="lower left")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "wall_plan_recommended.png"), dpi=120); plt.close(fig)
    print("plan done"); sys.exit(0)
for cf, out in curves.items():
    picks = {"拐点(边际效率≥2)": knee(out, 2.0), "拐点(边际效率≥1)": knee(out, 1.0)}
    for target in (9.0, 13.0):
        c = min(out, key=lambda r: abs(r["fill_est_Mm3"] - target)); picks[f"填方≈{target:.0f}"] = c
    picks["库容最大(λ最小)"] = out[-1]
    seen = set()
    for lab, r in picks.items():
        if r["lam"] in seen: continue
        seen.add(r["lam"]); t0 = time.time(); rec = evaluate(r, cf, lab); EVAL.append(rec)
        print(f"精算 坝顶 {cf:.0f} {lab} λ {r['lam']}: 区 {rec['n_regions']} 段 {rec['n_segments']} 坝长 {rec['dam_len_m']:.0f} m 自由岸 {rec['free_len_m']:.0f} m Hmax {rec['H_max_m']:.0f} | 填 RCC {rec['fill_rcc_Mm3']:.1f} CFRD {rec['fill_cfrd_Mm3']:.1f}（估 {rec['fill_est_Mm3']:.1f}）| 净库 {rec['net_rcc_Mm3']:.1f}（估 {rec['storage_est_Mm3']:.1f}）漫 {rec['spill']} 漏 {rec['leak_acre']:.0f} ac | 再挖(RCC) 底 {rec.get('deep_floor_ft_rcc', float('nan')):.0f} 净 {rec.get('deep_net_Mm3_rcc', float('nan')):.1f} 挖 {rec.get('deep_cut_Mm3_rcc', float('nan')):.1f} 挖÷需 {rec.get('deep_cut_over_need_rcc', float('nan')):.2f}  ({time.time()-t0:.0f} s)")
cols = sorted({k for r in EVAL for k in r if not k.startswith("_") and k != "segments" and not k.startswith("curve")}) + ["curve_rcc", "curve_cfrd"]
with open(os.path.join(OUT, f"wall_eval_{COST_KIND}.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols + ["segments"]); w.writeheader()
    for r in EVAL: w.writerow({**{k: r.get(k, "") for k in cols}, "segments": json.dumps(r["segments"])})

# ------------------------------------------------------------- 图：曲线
fig, ax = plt.subplots(figsize=(11, 7)); cmap = plt.get_cmap("viridis")
for i, (cf, out) in enumerate(curves.items()):
    col = cmap(i / max(len(curves) - 1, 1)); xs = [r["fill_est_Mm3"] for r in out]; ys = [r["storage_est_Mm3"] for r in out]
    ax.plot(xs, ys, "o-", color=col, ms=3.5, lw=1.2, label=f"坝顶 {cf:.0f}：天然库容（4 m 网格估算）")
    kn = knee(out, 1.0); ax.scatter([kn["fill_est_Mm3"]], [kn["storage_est_Mm3"]], s=90, facecolors="none", edgecolors=col, lw=2)
    kn2 = knee(out, 2.0); ax.scatter([kn2["fill_est_Mm3"]], [kn2["storage_est_Mm3"]], s=140, facecolors="none", edgecolors=col, lw=1, ls="--")
    ev = [r for r in EVAL if r["crest_ft"] == cf and not r["spill"]]
    ax.scatter([r["fill_rcc_Mm3"] for r in ev], [r["net_rcc_Mm3"] for r in ev], marker="s", color=col, s=30)
    ev2 = [r for r in ev if "deep_net_Mm3_rcc" in r]
    ax.scatter([r["fill_rcc_Mm3"] for r in ev2], [r["deep_net_Mm3_rcc"] for r in ev2], marker="^", color=col, s=45)
    for r in ev2: ax.annotate(f"{r['label']}\n底{r['deep_floor_ft_rcc']:.0f} 挖{r['deep_cut_Mm3_rcc']:.0f}", (r["fill_rcc_Mm3"], r["deep_net_Mm3_rcc"]), fontsize=6.5, xytext=(4, 2), textcoords="offset points", color=col)
    me = marg(out)
    for r, e in zip(out, me):
        if not np.isnan(e) and r["fill_est_Mm3"] > 0.3: ax.annotate(f"{e:.1f}", (r["fill_est_Mm3"], r["storage_est_Mm3"]), fontsize=5.5, color=col, xytext=(-2, -9), textcoords="offset points")
ax.scatter([B1["fill"]], [B1["net_nocut"]], marker="*", s=200, color="red", label="B1 全环 RCC 4,050：不挖 26.5 / 挖后 61.4 Mm³")
ax.scatter([B1["fill"]], [B1["net"]], marker="*", s=200, color="darkred")
ax.scatter([GULLY["fill_rcc"]], [GULLY["net"]], marker="D", s=60, color="k", label="纯沟谷坝 S29 4,000（RCC 2.1 Mm³，6.4 Mm³）")
ax.add_patch(Rectangle((9, 30), 4, 45, facecolor="green", alpha=0.08, edgecolor="green", ls="--")); ax.text(9.2, 31, "目标区：填 9–13，库 ≥ 30", fontsize=8, color="green")
ax.set_xlabel("累计填方 Mm³（RCC 断面；曲线为 4 m 网格估算，方块/三角为 1 m 精算）"); ax.set_ylabel("净库容 Mm³"); ax.grid(alpha=0.3)
ax.legend(fontsize=7.5, loc="upper left"); ax.set_title("墙段自由组合：累计填方 vs 累计库容（每条曲线一个坝顶；点旁小数字 = 该步边际效率 Δ库容/Δ填方；实心圈 = 边际效率降到 1 的停止点，虚圈 = 降到 2；方块 = 1 m 精算天然库容；三角 = 再挖后）", fontsize=8.5)
fig.tight_layout(); fig.savefig(os.path.join(OUT, f"wall_curves_{COST_KIND}.png"), dpi=120); plt.close(fig)

# ------------------------------------------------------------- 图：拐点方案平面 + 生长过程
sl = t.slope_deg(); w1, _ = t.mask_inside(parcel.buffer(60)); ext = [t.xc[w1[1]][0], t.xc[w1[1]][-1], t.yc[w1[0]][-1], t.yc[w1[0]][0]]
def plan(ax, rec, title):
    ax.imshow(sl[w1], cmap="Greys", extent=ext, vmin=0, vmax=55, alpha=0.65)
    ax.contour(t.xc, t.yc, t.z / FT, levels=np.arange(2600, 4400, 100), colors="#a08000", linewidths=0.3)
    px_, py_ = parcel.exterior.xy; ax.plot(px_, py_, "k--", lw=1)
    comp, (R0, R1, C0, C1) = rec["_comp"]; e2 = [t.xc[C0], t.xc[C1 - 1], t.yc[R1 - 1], t.yc[R0]]
    if "_pit" in rec:
        inside, wmax, F_ft, _ = rec["_pit"]; ZW = t.z[R0:R1, C0:C1]; z_new = np.where(inside, np.minimum(ZW, np.maximum(F_ft * FT, wmax)), ZW); dz = np.maximum(ZW - z_new, 0)
        ax.imshow(np.ma.masked_where(dz <= 0.01, dz), cmap="YlOrRd", alpha=0.85, extent=e2, vmin=0, vmax=max(60, float(np.nanmax(dz))))
    ax.imshow(np.ma.masked_where(~comp, np.ones_like(comp, float)), cmap=matplotlib.colors.ListedColormap(["#2c7fb8"]), alpha=0.35, extent=e2)
    for poly in rec["_polys"]: xx, yy = poly.exterior.xy; ax.plot(xx, yy, "-", color="#1b5e20", lw=1.2)
    for wl in rec["_walls"]:
        if wl.geom_type == "LineString": xx, yy = wl.xy; ax.plot(xx, yy, "-", color="#e65100", lw=3)
    for s in rec["segments"]:
        if s["len"] >= 60: ax.text(s["E"], s["N"], f"{s['len']:.0f} m/H{s['Hmax']:.0f}/{s['fill_rcc']:.1f}", fontsize=6, color="#bf360c", bbox=dict(fc="white", ec="none", alpha=0.7, pad=0.5))
    ax.set_xlim(685850, ext[1]); ax.set_ylim(ext[2], 4357100); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9)
    ax.legend(handles=[plt.Line2D([], [], color="#e65100", lw=3, label="墙段（坝，标注：长/最大坝高/RCC 填方 Mm³）"), plt.Line2D([], [], color="#1b5e20", lw=1.2, label="库区边界（其余为高于坝顶的天然岸）"),
                       Patch(facecolor="#2c7fb8", alpha=0.4, label="天然蓄水面"), Patch(facecolor="#fd8d3c", label="再挖深度（碗形坑，深挖参考点）"), plt.Line2D([], [], color="k", ls="--", label="宗地界")], fontsize=6.5, loc="lower left")
knees = [r for r in EVAL if r["label"] == "拐点(边际效率≥1)"]
n = len(knees); fig, axs = plt.subplots(1, n, figsize=(7.5 * n, 8)); axs = np.atleast_1d(axs)
for ax, rec in zip(axs, knees):
    plan(ax, rec, f"坝顶 {rec['crest_ft']:.0f} 拐点方案：{rec['n_segments']} 段墙共 {rec['dam_len_m']:.0f} m（自由岸 {rec['free_len_m']:.0f} m），填 RCC {rec['fill_rcc_Mm3']:.1f}/CFRD {rec['fill_cfrd_Mm3']:.1f} Mm³\n天然净库容 {rec['net_rcc_Mm3']:.1f} Mm³；再挖到 {rec.get('deep_floor_ft_rcc', float('nan')):.0f} ft：挖 {rec.get('deep_cut_Mm3_rcc', float('nan')):.1f}，净 {rec.get('deep_net_Mm3_rcc', float('nan')):.1f} Mm³")
fig.tight_layout(); fig.savefig(os.path.join(OUT, f"wall_plan_knee_{COST_KIND}.png"), dpi=110); plt.close(fig)
# 生长过程（最高坝顶）
cf = CRESTS[-1]; out = curves[cf]; steps = [r for r in out if r["fill_est_Mm3"] > 0.05]
pick = [steps[k] for k in sorted(set([0, len(steps) // 4, len(steps) // 2, (3 * len(steps)) // 4, len(steps) - 1]))]
fig, axs = plt.subplots(1, len(pick), figsize=(4.2 * len(pick), 5)); axs = np.atleast_1d(axs)
for ax, r in zip(axs, pick):
    ax.imshow(sl[w1], cmap="Greys", extent=ext, vmin=0, vmax=55, alpha=0.6); ax.plot(*parcel.exterior.xy, "k--", lw=0.8)
    m = np.ma.masked_where(~r["mask"], np.ones(r["mask"].shape)); ax.imshow(m, cmap=matplotlib.colors.ListedColormap(["#2c7fb8"]), alpha=0.5, extent=[xg[0] - K / 2, xg[-1] + K / 2, yg[-1] - K / 2, yg[0] + K / 2])
    ax.set_xlim(685850, ext[1]); ax.set_ylim(ext[2], 4357100); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"λ={r['lam']}：填 {r['fill_est_Mm3']:.1f} Mm³，库 {r['storage_est_Mm3']:.1f} Mm³\n坝 {r['dam_len_m']:.0f} m，自由岸 {r['free_len_m']:.0f} m", fontsize=8)
fig.suptitle(f"坝顶 {cf:.0f}：λ（填方价格）降低时库区的生长（蓝 = 库区；只留效率 ≥ 1/λ 的墙段）", fontsize=10); fig.tight_layout()
fig.savefig(os.path.join(OUT, f"wall_growth_{COST_KIND}.png"), dpi=110); plt.close(fig)
print("done")
