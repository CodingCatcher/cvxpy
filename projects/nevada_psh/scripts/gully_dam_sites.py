#!/usr/bin/env python3
"""跨沟谷短坝的坝址搜索（团队思路：不围整圈，找天然窄口/沟谷，用一条短坝围出上库）。

阶段 1（筛选，4 m DEM）：D8 汇流找沟谷线 → 沿沟谷每 30 m 一个候选坝址 → 坝轴垂直于沟谷（±30° 内取坝长最短的方向）→
  每个坝顶高程：两岸坝肩（地面 ≥ 坝顶处）、坝长、最大坝高、粗估填方（平地断面 × 坡地放大）、坝后淹没库容（洪水填充：
  水位以下且与坝上游相连、不越过坝轴的像元；若连到窗口边缘 = 漫过别处的鞍部，判为"关不住"）、淹没区出宗地的面积。
阶段 2（精算，1 m DEM）：对入围坝址逐站放坝（terrain_model.dam_section，CFRD 1.4/1.4 与 RCC 0/0.8）算填方与上游楔，
  1 m 洪水填充库容，并按笔记 14 的碗形坑规则算"坝后再挖"的库容与挖方。
输出：outputs/gully_stage1.csv、gully_stage2.csv、gully_sites_plan.png、gully_fill_vs_storage.png、gully_thalwegs.png。
用法：python3 scripts/gully_dam_sites.py [--stage1-only] [--top N]"""
from __future__ import annotations
import csv, json, math, os, sys
import numpy as np, shapely
from shapely.geometry import LineString, Point, Polygon
from scipy import ndimage
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
from terrain_model import Terrain, SECTIONS, FT, ACRE, DATA, OUT, Bowl, cross_section, dam_section, load_parcel  # noqa
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
matplotlib.rcParams["font.family"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]; matplotlib.rcParams["axes.unicode_minus"] = False
import logging; logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

STAGE1_ONLY = "--stage1-only" in sys.argv
TOP_N = int([a for a in sys.argv if a.startswith("--top")][0].split("=")[1]) if any(a.startswith("--top=") for a in sys.argv) else 8
LOWER_FT = 2400.0; ETA = 0.85; FREEBOARD_FT = 20.0
ZMIN_SITE_FT = float([a for a in sys.argv if a.startswith("--zmin=")][0].split("=")[1]) if any(a.startswith("--zmin=") for a in sys.argv) else 3650.0   # 坝址（沟底）高程下限
SUF = "" if ZMIN_SITE_FT >= 3650 else f"_zmin{ZMIN_SITE_FT:.0f}"
CREST_MAX_FT = 4120.0
ACC_MIN_CELLS = 600              # 4 m 像元：≈1 ha 汇水面积以上才算沟谷线
SITE_SPACING = 30.0              # 沿沟谷候选间距 m
MAX_HALF_LEN = 700.0             # 单侧坝肩搜索距离 m
BERM = 20.0; WALL_M = 0.75; SF = 1.2; AGG = 0.9
B1 = {"name": "B1 全环 RCC 4,050", "fill": 29.7, "need": 26.7, "net_nocut": 26.5, "net": 61.4, "nwl_ft": 4030.0, "cut": 36.3}

t = Terrain(os.path.join(DATA, "dem_3dep_1m_utm10.tif"), os.path.join(DATA, "dem_3dep_1m_utm10.json"))
parcel = load_parcel(); pb = parcel.bounds
gwh = lambda v_m3, nwl_ft: v_m3 * 1000 * 9.81 * (nwl_ft - LOWER_FT) * FT * ETA / 3.6e12

# ------------------------------------------------------------------ 4 m DEM 与 D8
K = 4
ny4, nx4 = t.ny // K, t.nx // K
z4 = t.z[:ny4 * K, :nx4 * K].reshape(ny4, K, nx4, K).mean(axis=(1, 3))
z4s = ndimage.uniform_filter(np.nan_to_num(z4, nan=np.nanmin(z4)), 3)
x4 = t.xmin + K * (0.5 + np.arange(nx4)); y4 = t.ymax - K * (0.5 + np.arange(ny4))
def d8(z):
    ny, nx = z.shape; best = np.full(z.shape, -np.inf); down = np.full(z.shape, -1, dtype=np.int64)
    idx = np.arange(z.size).reshape(z.shape)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dy == 0 and dx == 0: continue
            zs = np.full(z.shape, np.nan); js = np.full(z.shape, -1, dtype=np.int64)
            ys, ye = max(0, -dy), ny - max(0, dy); xs, xe = max(0, -dx), nx - max(0, dx)
            zs[ys:ye, xs:xe] = z[ys + dy:ye + dy, xs + dx:xe + dx]; js[ys:ye, xs:xe] = idx[ys + dy:ye + dy, xs + dx:xe + dx]
            slope = (z - zs) / math.hypot(dx, dy)
            m = np.isfinite(slope) & (slope > best) & (slope > 0)
            best[m] = slope[m]; down[m] = js[m]
    return down.ravel()
down4 = d8(z4s)
order = np.argsort(z4s.ravel())[::-1]; acc = np.ones(z4s.size)
for i in order:
    j = down4[i]
    if j >= 0: acc[j] += acc[i]
acc = acc.reshape(z4s.shape)
X4, Y4 = np.meshgrid(x4, y4)
in_parcel4 = shapely.contains_xy(parcel, X4.ravel(), Y4.ravel()).reshape(z4s.shape)
thalweg = (acc >= ACC_MIN_CELLS) & in_parcel4 & (z4s / FT >= ZMIN_SITE_FT - 30)
print("4 m DEM", z4s.shape, "沟谷像元", int(thalweg.sum()))

# 候选点：沿沟谷线按间距抽稀（按汇流量从大到小）
cand = []
cells = np.argwhere(thalweg); cells = cells[np.argsort(-acc[thalweg])]
for r, c in cells:
    p = np.array([x4[c], y4[r]])
    if all(np.hypot(*(p - q["p"])) >= SITE_SPACING for q in cand):
        cand.append({"r": r, "c": c, "p": p, "z_ft": float(z4s[r, c] / FT), "acc_ha": float(acc[r, c] * K * K / 1e4)})
cand = [q for q in cand if q["z_ft"] >= ZMIN_SITE_FT]
print("候选坝址", len(cand))

def valley_dir(q):
    """沟谷方向：沿 D8 向下游走 12 步（≈50 m）的位移方向；走不了就用坡度方向。"""
    i = q["r"] * nx4 + q["c"]; k = 0; j = i
    while k < 12 and down4[j] >= 0: j = down4[j]; k += 1
    if k >= 3:
        r2, c2 = divmod(j, nx4); v = np.array([x4[c2] - x4[q["c"]], y4[r2] - y4[q["r"]]])
    else:
        gy, gx = np.gradient(z4s, K); v = -np.array([gx[q["r"], q["c"]], gy[q["r"], q["c"]]])
    return v / max(np.hypot(*v), 1e-9)

def abutments(p, u, crest_m, step=2.0):
    """沿坝轴方向 ±u 找地面 ≥ 坝顶的坝肩；返回 (左距, 右距) 或 None。"""
    out = []
    for sgn in (+1, -1):
        d = np.arange(step, MAX_HALF_LEN + step, step); E = p[0] + sgn * d * u[0]; N = p[1] + sgn * d * u[1]
        z = t.elev(E, N); k = np.where(z >= crest_m)[0]
        if len(k) == 0 or not np.isfinite(z[:max(k[0], 1)]).all(): return None
        out.append(float(d[k[0]]))
    return out

# 洪水填充窗口（4 m）：宗地外扩 200 m
c0 = max(0, int((pb[0] - 200 - t.xmin) / K)); c1 = min(nx4, int((pb[2] + 200 - t.xmin) / K)); r0 = max(0, int((t.ymax - pb[3] - 200) / K)); r1 = min(ny4, int((t.ymax - pb[1] + 200) / K))
zw = z4[r0:r1, c0:c1]; inpw = in_parcel4[r0:r1, c0:c1]
def flood4(p, u, half, nwl_m, v_up):
    """4 m 洪水填充：水位以下、与坝上游种子连通、不越过坝轴（坝轴栅格化为墙）。返回 (毛库容 m³, 面积 m², 出界面积 m², 是否漫出窗口)。"""
    wall = np.zeros(zw.shape, bool)
    for d in np.arange(-half[1] - 8, half[0] + 8, 2.0):
        E = p[0] + d * u[0]; N = p[1] + d * u[1]; cc = int((E - t.xmin) / K) - c0; rr = int((t.ymax - N) / K) - r0
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                if 0 <= rr + dr < zw.shape[0] and 0 <= cc + dc < zw.shape[1]: wall[rr + dr, cc + dc] = True
    wet = (zw < nwl_m) & ~wall & np.isfinite(zw)
    lab, n = ndimage.label(wet)
    seed = p + 10.0 * v_up; sc = int((seed[0] - t.xmin) / K) - c0; sr = int((t.ymax - seed[1]) / K) - r0
    if not (0 <= sr < zw.shape[0] and 0 <= sc < zw.shape[1]) or lab[sr, sc] == 0:
        # 种子落在墙上或水位以上：向上游多走几步
        for k in (14, 18, 24):
            seed = p + k * v_up; sc = int((seed[0] - t.xmin) / K) - c0; sr = int((t.ymax - seed[1]) / K) - r0
            if 0 <= sr < zw.shape[0] and 0 <= sc < zw.shape[1] and lab[sr, sc] > 0: break
        else:
            return 0.0, 0.0, 0.0, False
    comp = lab == lab[sr, sc]
    spill = bool(comp[0, :].any() or comp[-1, :].any() or comp[:, 0].any() or comp[:, -1].any())
    gross = float(np.nansum(np.where(comp, nwl_m - zw, 0)) * K * K)
    return gross, float(comp.sum() * K * K), float((comp & ~inpw).sum() * K * K), spill

def rough_fill(p, u, half, crest_m, v_up, kind):
    """粗估填方：沿坝轴每 10 m 一站，H_i = 坝顶 − 地面；平地断面 ½H²(m_us+m_ds)+bH，下游侧按沟谷纵坡放大 1/(1−m·tanθ)。"""
    spec = SECTIONS[kind]; d = np.arange(-half[1], half[0] + 10, 10.0); E = p[0] + d * u[0]; N = p[1] + d * u[1]
    H = np.maximum(crest_m - t.elev(E, N), 0)
    # 沟谷纵坡（下游 60 m 内）
    zd = t.elev(p[0] - 60 * v_up[0], p[1] - 60 * v_up[1]); tanth = max(0.0, float((t.elev(*p) - zd) / 60.0))
    fd = 1 / max(1 - spec["m_ds"] * tanth, 0.25); fu = 1.0
    A = 0.5 * H ** 2 * (spec["m_ds"] * fd + spec["m_us"] * fu) + spec["crest_w"] * H
    return float(np.sum(A) * 10.0), tanth

# ------------------------------------------------------------------ 阶段 1
rows1 = []
for qi, q in enumerate(cand):
    v = valley_dir(q); v_up = -v; p = q["p"]; zp = q["z_ft"] * FT
    crests = np.arange(math.ceil((q["z_ft"] + 100) / 50) * 50, CREST_MAX_FT + 1, 50)
    # 坝轴方向：±30° 内取参考坝顶（沟底 + 60 m）下最短的
    best_u = None; best_L = np.inf
    for ang in (-30, -15, 0, 15, 30):
        a = math.radians(ang); u = np.array([v[0] * math.cos(a) - v[1] * math.sin(a), v[0] * math.sin(a) + v[1] * math.cos(a)]); u = np.array([-u[1], u[0]])
        ab = abutments(p, u, zp + 60)
        if ab and sum(ab) < best_L: best_L, best_u = sum(ab), u
    if best_u is None: continue
    u = best_u
    for cf in crests:
        cm = cf * FT; ab = abutments(p, u, cm)
        if ab is None: continue
        L = ab[0] + ab[1]; nwl = cm - FREEBOARD_FT * FT
        gross, area, out_area, spill = flood4(p, u, ab, nwl, v_up)
        fc, tanth = rough_fill(p, u, ab, cm, v_up, "cfrd"); fr, _ = rough_fill(p, u, ab, cm, v_up, "rcc")
        rows1.append({"site": qi, "E": p[0], "N": p[1], "z_floor_ft": q["z_ft"], "acc_ha": q["acc_ha"], "valley_grade_tan": tanth,
                      "ux": u[0], "uy": u[1], "vx": v_up[0], "vy": v_up[1], "crest_ft": cf, "nwl_ft": cf - FREEBOARD_FT, "H_max_m": cm - zp,
                      "abut_L": ab[0], "abut_R": ab[1], "dam_len_m": L, "fill_cfrd_Mm3": fc / 1e6, "fill_rcc_Mm3": fr / 1e6,
                      "gross_Mm3": gross / 1e6, "water_acre": area / ACRE, "water_out_acre": out_area / ACRE, "spill": spill,
                      "abut_in_parcel": bool(parcel.contains(Point(p + ab[0] * u)) and parcel.contains(Point(p - ab[1] * u))),
                      "eff_cfrd": (gross / fc) if fc > 0 else 0.0, "head_m": (cf - FREEBOARD_FT - LOWER_FT) * FT})
cols = list(rows1[0].keys())
with open(os.path.join(OUT, f"gully_stage1{SUF}.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows1)
ok = [r for r in rows1 if not r["spill"]]
print(f"阶段 1：{len(rows1)} 个坝址×坝顶，其中关得住 {len(ok)}")
top_by_storage = sorted(ok, key=lambda r: -r["gross_Mm3"])[:15]
for r in top_by_storage:
    print(f"  site {r['site']:3d} 沟底 {r['z_floor_ft']:.0f} 坝顶 {r['crest_ft']:.0f} H {r['H_max_m']:4.0f} L {r['dam_len_m']:4.0f} 填 CFRD {r['fill_cfrd_Mm3']:5.2f} RCC {r['fill_rcc_Mm3']:5.2f} 库 {r['gross_Mm3']:5.2f} Mm³ 水面 {r['water_acre']:4.0f} ac 出界 {r['water_out_acre']:3.0f} ac 坝肩在界内 {r['abut_in_parcel']}")

# 沟谷线图
fig, ax = plt.subplots(figsize=(11, 9))
sl = t.slope_deg(); w1, _ = t.mask_inside(parcel.buffer(150))
ext = [t.xc[w1[1]][0], t.xc[w1[1]][-1], t.yc[w1[0]][-1], t.yc[w1[0]][0]]
ax.imshow(sl[w1], cmap="Greys", extent=ext, vmin=0, vmax=55, alpha=0.7)
ax.contour(t.xc, t.yc, t.z / FT, levels=np.arange(2600, 4400, 100), colors="#a08000", linewidths=0.3)
ax.contour(t.xc, t.yc, t.z / FT, levels=[3700], colors="cyan", linewidths=0.8)
tm = np.ma.masked_where(~thalweg, np.log10(acc)); ax.imshow(tm, cmap="Blues", extent=[x4[0] - K / 2, x4[-1] + K / 2, y4[-1] - K / 2, y4[0] + K / 2], alpha=0.9, vmin=2, vmax=5)
px_, py_ = parcel.exterior.xy; ax.plot(px_, py_, "k--", lw=1)
best_site = {}
for r in ok:
    if r["site"] not in best_site or r["gross_Mm3"] > best_site[r["site"]]["gross_Mm3"]: best_site[r["site"]] = r
for r in best_site.values():
    p = np.array([r["E"], r["N"]]); u = np.array([r["ux"], r["uy"]])
    ax.plot([p[0] - r["abut_R"] * u[0], p[0] + r["abut_L"] * u[0]], [p[1] - r["abut_R"] * u[1], p[1] + r["abut_L"] * u[1]], "-", color="red", lw=0.8, alpha=0.7)
    ax.text(p[0], p[1], f"{r['gross_Mm3']:.1f}", fontsize=6, color="darkred")
ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect("equal")
ax.set_title("沟谷线（蓝，D8 汇流 ≥ 1 ha）与各候选坝址在其最大库容坝顶下的坝轴（红），数字 = 毛库容 Mm³；青线 3,700 ft", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, f"gully_thalwegs{SUF}.png"), dpi=120); plt.close(fig)
if STAGE1_ONLY: sys.exit(0)

# ------------------------------------------------------------------ 阶段 2：入围坝址精算（1 m DEM）
# 入围：按"关得住"的最大毛库容取前 TOP_N 个不同坝址（相距 ≥ 120 m），外加库容/填方比最高的 3 个
def pick_sites():
    chosen = []
    for r in sorted(best_site.values(), key=lambda r: -r["gross_Mm3"]):
        if all(math.hypot(r["E"] - c["E"], r["N"] - c["N"]) >= 120 for c in chosen): chosen.append(r)
        if len(chosen) >= TOP_N: break
    for r in sorted([r for r in ok if r["gross_Mm3"] > 0.3], key=lambda r: -r["eff_cfrd"])[:3]:
        if all(math.hypot(r["E"] - c["E"], r["N"] - c["N"]) >= 120 for c in chosen): chosen.append(r)
    return chosen
SITES = pick_sites(); print("入围坝址", [s["site"] for s in SITES])

# 1 m 洪水填充窗口
C0 = max(0, int(pb[0] - 200 - t.xmin)); C1 = min(t.nx, int(pb[2] + 200 - t.xmin)); R0 = max(0, int(t.ymax - pb[3] - 200)); R1 = min(t.ny, int(t.ymax - pb[1] + 200))
ZW = t.z[R0:R1, C0:C1]; XW, YW = np.meshgrid(t.xc[C0:C1], t.yc[R0:R1]); INP = shapely.contains_xy(parcel, XW.ravel(), YW.ravel()).reshape(ZW.shape)
import contourpy
def flood1(p, u, half, nwl_m, v_up, toe_us_max):
    wall = np.zeros(ZW.shape, bool)
    for d in np.arange(-half[1] - 6, half[0] + 6, 1.0):
        E = p[0] + d * u[0]; N = p[1] + d * u[1]; cc = int(E - t.xmin) - C0; rr = int(t.ymax - N) - R0
        wall[max(rr - 2, 0):rr + 3, max(cc - 2, 0):cc + 3] = True
    wet = (ZW < nwl_m) & ~wall & np.isfinite(ZW); lab, n = ndimage.label(wet)
    comp = None
    for k in (10, 16, 24, 40):
        seed = p + k * v_up; sc = int(seed[0] - t.xmin) - C0; sr = int(t.ymax - seed[1]) - R0
        if 0 <= sr < ZW.shape[0] and 0 <= sc < ZW.shape[1] and lab[sr, sc] > 0: comp = lab == lab[sr, sc]; break
    if comp is None: return None
    spill = bool(comp[0, :].any() or comp[-1, :].any() or comp[:, 0].any() or comp[:, -1].any())
    gross = float(np.nansum(np.where(comp, nwl_m - ZW, 0)))
    # 淹没区多边形
    gen = contourpy.contour_generator(x=t.xc[C0:C1], y=t.yc[R0:R1], z=comp.astype(float), name="serial")
    lines = [l for l in gen.lines(0.5) if len(l) > 3]; poly = max((Polygon(l).buffer(0) for l in lines), key=lambda g: g.area) if lines else None
    return {"comp": comp, "gross": gross, "area": float(comp.sum()), "out_area": float((comp & ~INP).sum()), "spill": spill, "poly": poly}

rows2 = []
for sinfo in SITES:
    p = np.array([sinfo["E"], sinfo["N"]]); u = np.array([sinfo["ux"], sinfo["uy"]]); v_up = np.array([sinfo["vx"], sinfo["vy"]]); zp = sinfo["z_floor_ft"] * FT
    crests = sorted({r["crest_ft"] for r in ok if r["site"] == sinfo["site"]})
    for cf in crests:
        cm = cf * FT; nwl = cm - FREEBOARD_FT * FT; ab = abutments(p, u, cm, step=1.0)
        if ab is None: continue
        # 逐站放坝：坝轴上每 10 m 一站，法向 = 上游方向
        d = np.arange(-ab[1], ab[0] + 1e-6, 10.0); res = {}
        for kind in ("cfrd", "rcc"):
            A = []; W = []; TU = []; TD = []
            for dd in d:
                sp = p + dd * u; off, zs = cross_section(t, sp, v_up); ds = dam_section(off, zs, cm, SECTIONS[kind], nwl)
                A.append(ds["area"]); W.append(ds["area_us_below_nwl"] if np.isfinite(ds["area_us_below_nwl"]) else 0.0)
                TU.append(ds["toe_us"] if ds["toe_us"] is not None else np.nan); TD.append(ds["toe_ds"] if ds["toe_ds"] is not None else np.nan)
            A = np.array(A); feas = np.isfinite(A)
            V = float(np.nansum([0.5 * (A[i] + A[i + 1]) * 10 for i in range(len(A) - 1) if feas[i] and feas[i + 1]]))
            Wd = float(np.nansum([0.5 * (W[i] + W[i + 1]) * 10 for i in range(len(A) - 1) if feas[i] and feas[i + 1]]))
            res[kind] = {"V": V, "wedge": Wd, "infeasible_m": float((~feas).sum() * 10), "toe_us_max": float(np.nanmax(TU)) if np.isfinite(TU).any() else 0.0, "toe_ds_max": float(np.nanmin(TD)) if np.isfinite(TD).any() else np.nan}
        fl = flood1(p, u, ab, nwl, v_up, res["cfrd"]["toe_us_max"])
        rec = {"site": sinfo["site"], "E": p[0], "N": p[1], "z_floor_ft": sinfo["z_floor_ft"], "crest_ft": cf, "nwl_ft": cf - FREEBOARD_FT, "H_max_m": cm - zp, "dam_len_m": ab[0] + ab[1],
               "abut_in_parcel": bool(parcel.contains(Point(p + ab[0] * u)) and parcel.contains(Point(p - ab[1] * u))), "valley_grade_tan": sinfo["valley_grade_tan"],
               "fill_cfrd_Mm3": res["cfrd"]["V"] / 1e6, "fill_rcc_Mm3": res["rcc"]["V"] / 1e6, "cfrd_infeasible_m": res["cfrd"]["infeasible_m"], "rcc_infeasible_m": res["rcc"]["infeasible_m"],
               "cfrd_toe_ds_m": res["cfrd"]["toe_ds_max"], "rcc_toe_ds_m": res["rcc"]["toe_ds_max"], "wedge_cfrd_Mm3": res["cfrd"]["wedge"] / 1e6, "wedge_rcc_Mm3": res["rcc"]["wedge"] / 1e6,
               "head_m": (cf - FREEBOARD_FT - LOWER_FT) * FT}
        if fl is None:
            rec.update({"gross_Mm3": 0.0, "spill": True}); rows2.append(rec); continue
        rec.update({"gross_Mm3": fl["gross"] / 1e6, "water_acre": fl["area"] / ACRE, "water_out_acre": fl["out_area"] / ACRE, "spill": fl["spill"],
                    "net_cfrd_Mm3": (fl["gross"] - res["cfrd"]["wedge"]) / 1e6, "net_rcc_Mm3": (fl["gross"] - res["rcc"]["wedge"]) / 1e6,
                    "net_cfrd_GWh": gwh(fl["gross"] - res["cfrd"]["wedge"], cf - FREEBOARD_FT), "net_rcc_GWh": gwh(fl["gross"] - res["rcc"]["wedge"], cf - FREEBOARD_FT)})
        # 坝后再挖（碗形坑）：坑顶线 = 淹没区多边形 − 坝轴两侧 (CFRD 上游趾 + 平台) 的带，再退 20 m；底面向下扫，取库容饱和 95% 点
        if fl["poly"] is not None and not fl["spill"]:
            axis = LineString([p - (ab[1] + 20) * u, p + (ab[0] + 20) * u])
            toe_poly = fl["poly"].difference(axis.buffer(max(res["cfrd"]["toe_us_max"], 10.0) + BERM))
            if toe_poly.geom_type != "Polygon": toe_poly = max(toe_poly.geoms, key=lambda g: g.area) if not toe_poly.is_empty else None
            if toe_poly is not None and toe_poly.area > 4 * ACRE:
                bowl = Bowl(t, toe_poly, BERM, WALL_M); inside, wmax = bowl.field(XW, YW, ZW)
                curve = []
                for F_ft in np.arange(cf - FREEBOARD_FT - 50, cf - FREEBOARD_FT - 700, -50):
                    Fm = F_ft * FT; z_new = np.where(inside, np.minimum(ZW, np.maximum(Fm, wmax)), ZW)
                    cut = float(np.nansum(np.maximum(ZW - z_new, 0))); g2 = float(np.nansum(np.where(fl["comp"] | inside, np.maximum(nwl - z_new, 0), 0)))
                    curve.append((F_ft, cut, g2))
                ginf = curve[-1][2]; deep = next((c for c in curve if c[2] >= 0.95 * ginf), curve[-1])
                need_c = res["cfrd"]["V"] / SF; need_r = res["rcc"]["V"] * AGG
                rec.update({"pit_top_acre": bowl.top.area / ACRE if not bowl.top.is_empty else 0.0, "deep_floor_ft": deep[0], "deep_cut_Mm3": deep[1] / 1e6,
                            "deep_gross_Mm3": deep[2] / 1e6, "deep_net_cfrd_Mm3": (deep[2] - res["cfrd"]["wedge"]) / 1e6, "deep_net_rcc_Mm3": (deep[2] - res["rcc"]["wedge"]) / 1e6,
                            "deep_net_rcc_GWh": gwh(deep[2] - res["rcc"]["wedge"], cf - FREEBOARD_FT), "cut_over_need_cfrd": deep[1] / need_c if need_c else np.nan, "cut_over_need_rcc": deep[1] / need_r if need_r else np.nan,
                            "curve": ";".join(f"{c[0]:.0f}:{c[1]/1e6:.2f}:{c[2]/1e6:.2f}" for c in curve)})
        rows2.append(rec)
        print(f"site {sinfo['site']:3d} 坝顶 {cf:.0f} H {rec['H_max_m']:4.0f} L {rec['dam_len_m']:4.0f} 填 CFRD {rec['fill_cfrd_Mm3']:5.2f} RCC {rec['fill_rcc_Mm3']:5.2f} 毛库 {rec['gross_Mm3']:5.2f} 净(RCC) {rec.get('net_rcc_Mm3', 0):5.2f} 漫 {rec['spill']} 出界 {rec.get('water_out_acre', 0):3.0f} ac | 再挖：底 {rec.get('deep_floor_ft', '-')} 挖 {rec.get('deep_cut_Mm3', float('nan')):.1f} 净 {rec.get('deep_net_rcc_Mm3', float('nan')):.1f}")
cols2 = sorted({k for r in rows2 for k in r}, key=lambda k: (k == "curve", k))
with open(os.path.join(OUT, f"gully_stage2{SUF}.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols2); w.writeheader(); w.writerows(rows2)

# ------------------------------------------------------------------ 图：填方 vs 库容；平面
fig, ax = plt.subplots(figsize=(10, 6.5))
okr = [r for r in rows2 if not r["spill"]]
sc = ax.scatter([r["fill_rcc_Mm3"] for r in okr], [r["net_rcc_Mm3"] for r in okr], c=[r["nwl_ft"] for r in okr], cmap="viridis", s=28, label="沟谷坝（RCC 断面），天然库容")
ax.scatter([r["fill_rcc_Mm3"] for r in okr if "deep_net_rcc_Mm3" in r], [r["deep_net_rcc_Mm3"] for r in okr if "deep_net_rcc_Mm3" in r], c=[r["nwl_ft"] for r in okr if "deep_net_rcc_Mm3" in r], cmap="viridis", s=28, marker="^", label="沟谷坝 + 坝后碗形坑（库容饱和 95% 点）")
for r in okr:
    if r["net_rcc_Mm3"] > 1.5: ax.annotate(f"S{r['site']}/{r['crest_ft']:.0f}", (r["fill_rcc_Mm3"], r["net_rcc_Mm3"]), fontsize=6, xytext=(3, 2), textcoords="offset points")
ax.scatter([B1["fill"]], [B1["net_nocut"]], marker="*", s=180, color="red", label="B1 全环 RCC 4,050 不挖（29.7 Mm³ 填，26.5 Mm³）")
ax.scatter([B1["fill"]], [B1["net"]], marker="*", s=180, color="darkred", label="B1 全环 RCC 4,050 + 碗形坑（61.4 Mm³）")
ax.axvline(B1["fill"] / 10, color="0.4", ls="--", lw=0.8); ax.text(B1["fill"] / 10, ax.get_ylim()[1] * 0.02 + 30, "填方 = B1 的 1/10", fontsize=8, rotation=90, va="bottom")
ax.axhline(B1["net"] / 2, color="0.4", ls="--", lw=0.8); ax.text(0.2, B1["net"] / 2 + 0.5, "库容 = B1 的 1/2", fontsize=8)
ax.set_xscale("log"); ax.set_xlabel("填方 Mm³（RCC 断面；CFRD 约 ×3–4）"); ax.set_ylabel("净库容 Mm³"); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=7.5, loc="upper left")
cb = plt.colorbar(sc, ax=ax); cb.set_label("水位 ft")
ax.set_title("沟谷坝候选：填方 vs 库容（每点一个坝址×坝顶），与全环 B1 对照", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, f"gully_fill_vs_storage{SUF}.png"), dpi=120); plt.close(fig)

fig, ax = plt.subplots(figsize=(11, 9))
ax.imshow(sl[w1], cmap="Greys", extent=ext, vmin=0, vmax=55, alpha=0.7)
ax.contour(t.xc, t.yc, t.z / FT, levels=np.arange(2600, 4400, 100), colors="#a08000", linewidths=0.3)
ax.contour(t.xc, t.yc, t.z / FT, levels=[3700], colors="cyan", linewidths=0.8); ax.plot(px_, py_, "k--", lw=1)
best2 = {}
for r in okr:
    if r["site"] not in best2 or r["net_rcc_Mm3"] > best2[r["site"]]["net_rcc_Mm3"]: best2[r["site"]] = r
cmap = plt.get_cmap("tab10")
for k, (sid, r) in enumerate(sorted(best2.items(), key=lambda kv: -kv[1]["net_rcc_Mm3"])):
    p = np.array([r["E"], r["N"]]); sinfo = next(s for s in SITES if s["site"] == sid); u = np.array([sinfo["ux"], sinfo["uy"]]); v_up = np.array([sinfo["vx"], sinfo["vy"]])
    ab = abutments(p, u, r["crest_ft"] * FT, step=1.0)
    fl = flood1(p, u, ab, r["nwl_ft"] * FT, v_up, 0)
    col = cmap(k % 10)
    if fl and fl["poly"] is not None:
        xx, yy = fl["poly"].exterior.xy; ax.fill(xx, yy, color=col, alpha=0.25, lw=0); ax.plot(xx, yy, "-", color=col, lw=0.6)
    ax.plot([p[0] - ab[1] * u[0], p[0] + ab[0] * u[0]], [p[1] - ab[1] * u[1], p[1] + ab[0] * u[1]], "-", color=col, lw=2.5)
    ax.text(p[0], p[1], f"S{sid}: 坝顶 {r['crest_ft']:.0f}, L {r['dam_len_m']:.0f} m, H {r['H_max_m']:.0f} m\n填 {r['fill_rcc_Mm3']:.1f}/{r['fill_cfrd_Mm3']:.1f} Mm³, 库 {r['net_rcc_Mm3']:.1f} Mm³", fontsize=6.5, color="k", bbox=dict(fc="white", ec=col, alpha=0.85, lw=0.8))
ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect("equal")
ax.set_title("入围沟谷坝址（各取净库容最大的坝顶）：坝轴（粗线）、淹没区（同色填充）；填方 = RCC/CFRD", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, f"gully_sites_plan{SUF}.png"), dpi=120); plt.close(fig)
print("done")
