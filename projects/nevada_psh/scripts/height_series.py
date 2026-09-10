#!/usr/bin/env python3
"""S29、S2 两个沟谷坝址在不同坝高下的坝长、填方、库容；与东半环（3,950 / 4,050）和全环 B1 放在同一张填方–库容图上。
输出：outputs/height_series.csv、hs_map.png、hs_curve.png，deliverables/D_height_series_brief.html。"""
from __future__ import annotations
import base64, csv, json, math, os, sys
import numpy as np, shapely, contourpy
from shapely.geometry import LineString, Point, Polygon
from shapely import wkt
from scipy import ndimage
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
from terrain_model import Terrain, SECTIONS, FT, ACRE, DATA, OUT, Bowl, cross_section, dam_section, densify, load_parcel  # noqa
DEL = os.path.join(ROOT, "deliverables")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
matplotlib.rcParams["font.family"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]; matplotlib.rcParams["axes.unicode_minus"] = False
import logging; logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

LOWER_FT = 2400.0; ETA = 0.85; FREEBOARD_FT = 20.0; BERM = 20.0; WALL_M = 0.75; AGG = 0.9
HEIGHTS = list(range(20, 131, 10)); MARK = [30, 50, 80, 100, 130]
gwh = lambda v_m3, nwl_ft: v_m3 * 1000 * 9.81 * (nwl_ft - LOWER_FT) * FT * ETA / 3.6e12
t = Terrain(os.path.join(DATA, "dem_3dep_1m_utm10.tif"), os.path.join(DATA, "dem_3dep_1m_utm10.json")); parcel = load_parcel(); pb = parcel.bounds
s1 = {int(r["site"]): r for r in csv.DictReader(open(os.path.join(OUT, "gully_stage1.csv")))}
SITES = {}
for sid in (29, 2):
    r = s1[sid]; SITES[sid] = {"p": np.array([float(r["E"]), float(r["N"])]), "u": np.array([float(r["ux"]), float(r["uy"])]), "v_up": np.array([float(r["vx"]), float(r["vy"])]), "floor_ft": float(r["z_floor_ft"])}
C0 = max(0, int(pb[0] - 200 - t.xmin)); C1 = min(t.nx, int(pb[2] + 200 - t.xmin)); R0 = max(0, int(t.ymax - pb[3] - 200)); R1 = min(t.ny, int(t.ymax - pb[1] + 200))
ZW = t.z[R0:R1, C0:C1]; XW, YW = np.meshgrid(t.xc[C0:C1], t.yc[R0:R1]); INP = shapely.contains_xy(parcel, XW.ravel(), YW.ravel()).reshape(ZW.shape)

def abutments(p, u, crest_m, maxlen=900.0):
    out = []
    for sgn in (+1, -1):
        d = np.arange(1.0, maxlen + 1, 1.0); z = t.elev(p[0] + sgn * d * u[0], p[1] + sgn * d * u[1]); k = np.where(z >= crest_m)[0]
        if len(k) == 0: return None
        out.append(float(d[k[0]]))
    return out

def flood(p, u, half, nwl_m, v_up):
    wall = np.zeros(ZW.shape, bool)
    for d in np.arange(-half[1] - 6, half[0] + 6, 1.0):
        cc = int(p[0] + d * u[0] - t.xmin) - C0; rr = int(t.ymax - (p[1] + d * u[1])) - R0; wall[max(rr - 2, 0):rr + 3, max(cc - 2, 0):cc + 3] = True
    wet = (ZW < nwl_m) & ~wall & np.isfinite(ZW); lab, _ = ndimage.label(wet)
    for k in (10, 16, 24, 40):
        seed = p + k * v_up; sc = int(seed[0] - t.xmin) - C0; sr = int(t.ymax - seed[1]) - R0
        if 0 <= sr < ZW.shape[0] and 0 <= sc < ZW.shape[1] and lab[sr, sc] > 0:
            comp = lab == lab[sr, sc]; spill = bool(comp[0, :].any() or comp[-1, :].any() or comp[:, 0].any() or comp[:, -1].any())
            gen = contourpy.contour_generator(x=t.xc[C0:C1], y=t.yc[R0:R1], z=comp.astype(float), name="serial")
            lines = [l for l in gen.lines(0.5) if len(l) > 3]; poly = max((Polygon(l).buffer(0) for l in lines), key=lambda g: g.area) if lines else None
            return {"comp": comp, "gross": float(np.nansum(np.where(comp, nwl_m - ZW, 0))), "area": float(comp.sum()), "out": float((comp & ~INP).sum()), "spill": spill, "poly": poly}
    return None

rows = []
for sid, S in SITES.items():
    p, u, v_up = S["p"], S["u"], S["v_up"]; floor_m = S["floor_ft"] * FT
    for H in HEIGHTS:
        cm = floor_m + H; cf = cm / FT; nwl = cm - FREEBOARD_FT * FT; ab = abutments(p, u, cm)
        rec = {"site": f"S{sid}", "H_m": H, "crest_ft": cf, "nwl_ft": cf - FREEBOARD_FT, "floor_ft": S["floor_ft"]}
        if ab is None:
            rec.update({"feasible": False}); rows.append(rec); print(f"S{sid} H {H}: 无坝肩"); continue
        d = np.arange(-ab[1], ab[0] + 1e-6, 10.0); A = []; W = []; TU = []; TD = []; out_n = 0
        for dd in d:
            sp = p + dd * u; off, zs = cross_section(t, sp, v_up); ds = dam_section(off, zs, cm, SECTIONS["rcc"], nwl)
            A.append(ds["area"]); W.append(ds["area_us_below_nwl"] if np.isfinite(ds["area_us_below_nwl"]) else 0.0)
            TU.append(ds["toe_us"] if ds["toe_us"] is not None else np.nan); TD.append(ds["toe_ds"] if ds["toe_ds"] is not None else np.nan)
            if ds["toe_ds"] is not None and not parcel.contains(Point(sp + ds["toe_ds"] * v_up)): out_n += 1
        A = np.array(A); ok = np.isfinite(A)
        V = float(np.nansum([0.5 * (A[i] + A[i + 1]) * 10 for i in range(len(A) - 1) if ok[i] and ok[i + 1]])); Wd = float(np.nansum([0.5 * (W[i] + W[i + 1]) * 10 for i in range(len(A) - 1) if ok[i] and ok[i + 1]]))
        fl = flood(p, u, ab, nwl, v_up)
        rec.update({"feasible": True, "dam_len_m": ab[0] + ab[1], "abut_in_parcel": bool(parcel.contains(Point(p + ab[0] * u)) and parcel.contains(Point(p - ab[1] * u))), "fill_rcc_Mm3": V / 1e6, "need_Mm3": V * AGG / 1e6,
                    "wedge_Mm3": Wd / 1e6, "toe_out_m": out_n * 10.0})
        if fl is None or fl["spill"]:
            rec.update({"net_Mm3": np.nan, "spill": True}); rows.append(rec); print(f"S{sid} H {H}: 关不住"); continue
        net = fl["gross"] - Wd
        rec.update({"gross_Mm3": fl["gross"] / 1e6, "net_Mm3": net / 1e6, "water_acre": fl["area"] / ACRE, "water_out_acre": fl["out"] / ACRE, "spill": False, "gwh_nat": gwh(net, cf - FREEBOARD_FT)})
        # 再挖：碗形坑，坑顶线 = 淹没区 − 坝轴两侧 (RCC 上游趾 4 m + 20 m 平台)
        axis = LineString([p - (ab[1] + 20) * u, p + (ab[0] + 20) * u]); tp = fl["poly"].difference(axis.buffer(max(np.nanmax(TU) if np.isfinite(TU).any() else 4.0, 4.0) + BERM))
        if tp.geom_type != "Polygon": tp = max(tp.geoms, key=lambda g: g.area) if not tp.is_empty else None
        if tp is not None and tp.area > 2 * ACRE:
            bowl = Bowl(t, tp, BERM, WALL_M); inside, wmax = bowl.field(XW, YW, ZW); curve = []
            for F_ft in np.arange(cf - FREEBOARD_FT - 50, cf - FREEBOARD_FT - 700, -50):
                Fm = F_ft * FT; z_new = np.where(inside, np.minimum(ZW, np.maximum(Fm, wmax)), ZW)
                curve.append((F_ft, float(np.nansum(np.maximum(ZW - z_new, 0))), float(np.nansum(np.where(fl["comp"] | inside, np.maximum(nwl - z_new, 0), 0)))))
            ginf = curve[-1][2]; deep = next((c for c in curve if c[2] >= 0.95 * ginf), curve[-1])
            rec.update({"deep_floor_ft": deep[0], "deep_cut_Mm3": deep[1] / 1e6, "deep_net_Mm3": (deep[2] - Wd) / 1e6, "gwh_deep": gwh(deep[2] - Wd, cf - FREEBOARD_FT), "cut_over_need": deep[1] / (V * AGG) if V else np.nan})
        else:
            rec.update({"deep_floor_ft": np.nan, "deep_cut_Mm3": 0.0, "deep_net_Mm3": net / 1e6, "gwh_deep": gwh(net, cf - FREEBOARD_FT), "cut_over_need": 0.0})
        rows.append(rec)
        print(f"S{sid} H {H:3d}: 坝顶 {cf:.0f} L {rec['dam_len_m']:4.0f} 填 {rec['fill_rcc_Mm3']:.2f} 净 {rec['net_Mm3']:.2f} 再挖 {rec['deep_net_Mm3']:.2f}（挖 {rec['deep_cut_Mm3']:.2f}）坝肩在界内 {rec['abut_in_parcel']}")
cols = sorted({k for r in rows for k in r}, key=lambda k: (k != "site", k != "H_m", k))
with open(os.path.join(OUT, "height_series.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)

# 东半环与 B1
ev = {(r["crest_ft"][:4], r["lambda"]): r for r in csv.DictReader(open(os.path.join(OUT, "wall_plan_eval.csv")))}
f_ = lambda v: float(v) if v not in ("", None) else float("nan")
def ring(label, r, H):
    nwl = f_(r["nwl_ft"]); return {"site": label, "H_m": H, "crest_ft": f_(r["crest_ft"]), "nwl_ft": nwl, "dam_len_m": f_(r["dam_len_m"]), "fill_rcc_Mm3": f_(r["fill_rcc_Mm3"]), "need_Mm3": f_(r["need_rcc_Mm3"]),
                                  "net_Mm3": f_(r["net_rcc_Mm3"]), "gwh_nat": gwh(f_(r["net_rcc_Mm3"]) * 1e6, nwl), "deep_floor_ft": f_(r["deep_floor_ft_rcc"]), "deep_cut_Mm3": f_(r["deep_cut_Mm3_rcc"]), "deep_net_Mm3": f_(r["deep_net_Mm3_rcc"]),
                                  "gwh_deep": gwh(f_(r["deep_net_Mm3_rcc"]) * 1e6, nwl), "cut_over_need": f_(r["deep_cut_over_need_rcc"]), "feasible": True, "spill": False}
HALF = [ring("东半环 3,950", ev[("3950", "1.6")], 131), ring("东半环 4,050（W1–W3）", ev[("4050", "2.1")], 131)]
B1 = {"site": "全环 B1 4,050", "H_m": 108, "crest_ft": 4050, "nwl_ft": 4030, "dam_len_m": 3900, "fill_rcc_Mm3": 29.7, "need_Mm3": 26.7, "net_Mm3": 26.5, "gwh_nat": gwh(26.5e6, 4030), "deep_floor_ft": 3450, "deep_cut_Mm3": 36.3, "deep_net_Mm3": 61.4, "gwh_deep": 70.6, "cut_over_need": 1.36, "feasible": True, "spill": False}

# ---------------------------------------------------------------- 图 1：地图
geoms = json.load(open(os.path.join(OUT, "wall_plan_geoms.json")))[0]
ringline = LineString(json.load(open(os.path.join(DATA, "dam_line_candidates.json")))["A_parcel_ring_3700"]["coords"])
sl = t.slope_deg(); w1, _ = t.mask_inside(parcel.buffer(60)); ext = [t.xc[w1[1]][0], t.xc[w1[1]][-1], t.yc[w1[0]][-1], t.yc[w1[0]][0]]
fig, ax = plt.subplots(figsize=(9.5, 10))
ax.imshow(sl[w1], cmap="Greys", extent=ext, vmin=0, vmax=55, alpha=0.6)
ax.contour(t.xc, t.yc, t.z / FT, levels=np.arange(2600, 4400, 100), colors="#a08000", linewidths=0.3)
ax.contour(t.xc, t.yc, t.z / FT, levels=[3700, 4030], colors=["cyan", "#7b1fa2"], linewidths=[0.6, 0.6], linestyles=["-", "--"])
ax.plot(*parcel.exterior.xy, "k--", lw=1)
ax.plot(*ringline.xy, "-", color="#d62728", lw=2.2, alpha=0.8, label="全环 B1：沿 3,700 ft 线的环形坝轴（3.9 km）")
names = {2449: "W2", 762: "W3", 4297: "W1"}
for wk in geoms["walls"]:
    w = wkt.loads(wk)
    if w.length < 50: continue
    ax.plot(*w.xy, "-", color="#e65100", lw=4, alpha=0.9); c = w.interpolate(0.5, normalized=True)
    nm = "W1" if w.length > 1000 else ("W2" if w.length > 500 else "W3"); ax.text(c.x + 25, c.y, nm, fontsize=10, fontweight="bold", color="#bf360c")
ax.plot([], [], "-", color="#e65100", lw=4, label="东半环 4,050：W1 1,370 m / W2 860 m / W3 280 m")
for sid, col in ((29, "#1565c0"), (2, "#2e7d32")):
    S = SITES[sid]; Hd = max(r["H_m"] for r in rows if r["site"] == f"S{sid}" and r.get("feasible") and not r.get("spill"))
    ab = abutments(S["p"], S["u"], S["floor_ft"] * FT + Hd); p, u = S["p"], S["u"]
    ax.plot([p[0] - ab[1] * u[0], p[0] + ab[0] * u[0]], [p[1] - ab[1] * u[1], p[1] + ab[0] * u[1]], "-", color=col, lw=3.5, label=f"S{sid} 沟谷坝轴（画的是最高可行坝高 {Hd} m 时的坝长 {ab[0]+ab[1]:.0f} m；沟底 {S['floor_ft']:,.0f} ft）")
    ax.plot(p[0], p[1], "o", color=col, ms=6); ax.text(p[0] + 30, p[1] - 40, f"S{sid}", fontsize=11, fontweight="bold", color=col)
ax.text(686250, 4356980, "3,700 ft 等高线（青）\n4,030 ft 水位（紫虚线）", fontsize=8, color="0.2")
ax.set_xlim(685900, 686900); ax.set_ylim(4355800, 4357100); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
ax.plot([685950, 686450], [4355850, 4355850], "k-", lw=2); ax.text(686200, 4355870, "500 m", ha="center", fontsize=8)
ax.legend(fontsize=8, loc="lower left"); ax.set_title("图 1 方案位置：S29、S2 沟谷坝轴，东半环的 W1/W2/W3，全环 B1（宗地界为黑虚线）", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "hs_map.png"), dpi=120); plt.close(fig)

# ---------------------------------------------------------------- 图 2：填方–库容连续谱
fig, ax = plt.subplots(figsize=(11, 7))
for sid, col in ((29, "#1565c0"), (2, "#2e7d32")):
    rr = [r for r in rows if r["site"] == f"S{sid}" and r.get("feasible") and not r.get("spill")]
    ax.plot([r["fill_rcc_Mm3"] for r in rr], [r["net_Mm3"] for r in rr], "-", color=col, lw=1.8, label=f"S{sid} 沟谷坝：天然库容（坝高 20→{max(r['H_m'] for r in rr)} m）")
    ax.plot([r["fill_rcc_Mm3"] for r in rr], [r["deep_net_Mm3"] for r in rr], ":", color=col, lw=1.4, label=f"S{sid} 沟谷坝：坝后再挖后")
    for r in rr:
        if r["H_m"] in MARK:
            ax.plot(r["fill_rcc_Mm3"], r["net_Mm3"], "o", color=col, ms=7); ax.plot(r["fill_rcc_Mm3"], r["deep_net_Mm3"], "o", mfc="none", mec=col, ms=7)
            ax.annotate(f"H{r['H_m']}", (r["fill_rcc_Mm3"], r["net_Mm3"]), fontsize=7.5, color=col, xytext=(3, -10), textcoords="offset points")
for h, mk, col in ((HALF[0], "s", "#6a1b9a"), (HALF[1], "s", "#e65100")):
    ax.plot(h["fill_rcc_Mm3"], h["net_Mm3"], mk, color=col, ms=10, label=f"{h['site']}：天然 {h['net_Mm3']:.1f} → 再挖后 {h['deep_net_Mm3']:.1f} Mm³（H {h['H_m']} m）")
    ax.plot(h["fill_rcc_Mm3"], h["deep_net_Mm3"], mk, mfc="none", mec=col, ms=10, mew=2); ax.plot([h["fill_rcc_Mm3"]] * 2, [h["net_Mm3"], h["deep_net_Mm3"]], "-", color=col, lw=0.8, alpha=0.6)
    ax.annotate(f"H{h['H_m']}", (h["fill_rcc_Mm3"], h["net_Mm3"]), fontsize=7.5, color=col, xytext=(4, -10), textcoords="offset points")
ax.plot(B1["fill_rcc_Mm3"], B1["net_Mm3"], "*", color="#b71c1c", ms=16, label=f"全环 B1 4,050：天然 26.5 → 再挖后 61.4 Mm³（H 108 m）"); ax.plot(B1["fill_rcc_Mm3"], B1["deep_net_Mm3"], "*", mfc="none", mec="#b71c1c", ms=16, mew=2)
ax.plot([B1["fill_rcc_Mm3"]] * 2, [B1["net_Mm3"], B1["deep_net_Mm3"]], "-", color="#b71c1c", lw=0.8, alpha=0.6); ax.annotate("H108", (B1["fill_rcc_Mm3"], B1["net_Mm3"]), fontsize=7.5, color="#b71c1c", xytext=(6, -10), textcoords="offset points")
# 连续谱：天然库容的上包络
pts = sorted([(r["fill_rcc_Mm3"], r["net_Mm3"]) for r in rows if r.get("feasible") and not r.get("spill")] + [(h["fill_rcc_Mm3"], h["net_Mm3"]) for h in HALF] + [(B1["fill_rcc_Mm3"], B1["net_Mm3"])])
env = []; best = -1
for x, y in pts:
    if y > best: env.append((x, y)); best = y
ax.plot([e[0] for e in env], [e[1] for e in env], "--", color="0.35", lw=1, label="连续谱：每个填方量能换到的最大天然库容（上包络）")
ax.set_xscale("log"); ax.set_xlim(0.05, 60); ax.set_ylim(0, 70); ax.set_xticks([0.1, 0.3, 1, 3, 10, 30]); ax.set_xticklabels(["0.1", "0.3", "1", "3", "10", "30"])
ax.set_xlabel("累计填方 Mm³（RCC 断面，对数轴）"); ax.set_ylabel("累计库容 Mm³（实心 = 天然，空心 = 坝后再挖后）"); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=7.5, loc="upper left")
ax.set_title("图 2 填方–库容连续谱：S29、S2 沟谷坝随坝高变化，东半环 3,950 / 4,050，全环 B1", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "hs_curve.png"), dpi=120); plt.close(fig)

# ---------------------------------------------------------------- HTML
img64 = lambda p_: "data:image/png;base64," + base64.b64encode(open(p_, "rb").read()).decode()
esc = lambda x: str(x).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
def table(cols_, rows_):
    h = ['<div class="wrap"><table><tr>' + "".join(f"<th>{esc(c)}</th>" for c in cols_) + "</tr>"]
    for r in rows_: h.append("<tr>" + "".join(f"<td>{esc(x)}</td>" for x in r) + "</tr>")
    return "\n".join(h) + "</table></div>"
CSS = """body{font-family:"Noto Sans CJK SC","WenQuanYi Zen Hei","PingFang SC","Microsoft YaHei",Arial,sans-serif;max-width:1100px;margin:0 auto;padding:18px 26px 50px;line-height:1.6;color:#1f2328;background:#fff}
h1{font-size:20px;border-bottom:2px solid #333;padding-bottom:6px} h2{font-size:16px;margin-top:26px;border-left:5px solid #7b1fa2;padding-left:10px}
table{border-collapse:collapse;font-size:13px;margin:8px 0 12px} th,td{border:1px solid #bfc4cc;padding:4px 8px;text-align:left} th{background:#eef1f5} tr:nth-child(even) td{background:#f8f9fb}
.wrap{overflow-x:auto} figure{margin:12px 0} img{max-width:100%;border:1px solid #d8dbe0} figcaption{font-size:12.5px;color:#444} .ans{background:#e8f5e9;border-left:4px solid #2e7d32;padding:10px 14px;margin:10px 0} .small{font-size:12px;color:#555}"""
def fmt_row(r):
    if not r.get("feasible"): return [r["site"], r["H_m"], f"{r['crest_ft']:,.0f} / {r['nwl_ft']:,.0f}", "无坝肩（两岸到不了坝顶）", "—", "—", "—", "—", "—"]
    if r.get("spill"): return [r["site"], r["H_m"], f"{r['crest_ft']:,.0f} / {r['nwl_ft']:,.0f}", f"{r['dam_len_m']:,.0f}", f"{r['fill_rcc_Mm3']:.2f}", "关不住", "—", "—", "—"]
    return [r["site"], r["H_m"], f"{r['crest_ft']:,.0f} / {r['nwl_ft']:,.0f}", f"{r['dam_len_m']:,.0f}" + ("" if r.get("abut_in_parcel", True) else "（坝肩出界）"), f"{r['fill_rcc_Mm3']:.2f}", f"{r['net_Mm3']:.2f}",
            f"{r['deep_net_Mm3']:.1f}（挖 {r['deep_cut_Mm3']:.1f}，底 {r['deep_floor_ft']:,.0f}）" if np.isfinite(r.get("deep_floor_ft", np.nan)) else f"{r['deep_net_Mm3']:.1f}（坑太小）", f"{r['gwh_nat']:.1f} / {r['gwh_deep']:.1f}", f"{r['cut_over_need']:.1f}"]
trows = [fmt_row(r) for r in rows if r["H_m"] in MARK] + [fmt_row(h) for h in HALF] + [fmt_row(B1)]
h = []
h.append("<h1>S29、S2 两个坝址：坝高 → 填方 → 库容（简版）</h1>")
h.append('<figure><img src="' + img64(os.path.join(OUT, "hs_map.png")) + '"><figcaption>图 1 位置。S29（蓝）在东北小沟，S2（绿）在向西北流的主沟；东半环的三段墙 W1/W2/W3（橙）；全环 B1 的坝轴（红）沿 3,700 ft 等高线绕山鼻一圈。</figcaption></figure>')
h.append('<figure><img src="' + img64(os.path.join(OUT, "hs_curve.png")) + '"><figcaption>图 2 连续谱。S29、S2 的线上每 10 m 坝高一个点，标出 30/50/80/100/130 m；实心 = 天然库容，空心 = 坝后再挖后；灰虚线是"每个填方量能换到的最大天然库容"。横轴对数。</figcaption></figure>')
h.append("<h2>表</h2>")
h.append(table(["方案", "坝高 m", "坝顶 / 水位 ft", "坝长 m", "RCC 填方 Mm³", "天然库容 Mm³", "挖深后总库容 Mm³", "GWh 天然 / 挖后", "挖 ÷ 需"], trows))
h.append('<p class="small">坝长 = 两岸地面到达坝顶处之间的距离。填方按 RCC 断面（直立 / 0.8:1）逐站放坝。天然库容 = 水位（坝顶 − 20 ft）以下坝后淹没区的体积，已扣坝体上游楔。挖深后 = 在淹没区内按碗形坑（坑壁 0.75:1、离坝 24 m）挖到库容饱和 95% 的库底。GWh 按下库 2,400 ft、η 0.85。挖 ÷ 需 = 挖方 ÷ RCC 骨料需求（RCC × 0.9）。东半环与 B1 的数字来自笔记 14、16。</p>')
# 结论数字
s29 = {r["H_m"]: r for r in rows if r["site"] == "S29" and r.get("feasible") and not r.get("spill")}; s2 = {r["H_m"]: r for r in rows if r["site"] == "S2" and r.get("feasible") and not r.get("spill")}
hmax29 = max(s29); hmax2 = max(s2); H4050 = HALF[1]; H3950 = HALF[0]
eff = lambda r: r["net_Mm3"] / r["fill_rcc_Mm3"]
h.append("<h2>结论</h2>")
h.append(f'<div class="ans"><b>效率拐点在东半环 4,050。</b>拐点以下，每 1 m³ 填方换 2–3 m³ 天然库容：S29 从坝高 30 到 {hmax29} m，填方 {s29[30]["fill_rcc_Mm3"]:.1f} → {s29[hmax29]["fill_rcc_Mm3"]:.1f} Mm³，库容 {s29[30]["net_Mm3"]:.1f} → {s29[hmax29]["net_Mm3"]:.1f} Mm³（效率 {eff(s29[30]):.1f}–{eff(s29[hmax29]):.1f}）；'
         f'S2 从 30 到 {hmax2} m，填方 {s2[30]["fill_rcc_Mm3"]:.1f} → {s2[hmax2]["fill_rcc_Mm3"]:.1f} Mm³，库容 {s2[30]["net_Mm3"]:.1f} → {s2[hmax2]["net_Mm3"]:.1f} Mm³（效率 {eff(s2[hmax2]):.1f}）；两个沟谷坝再高就没有坝肩或坝肩出界，库容封顶在 {max(s29[hmax29]["net_Mm3"], s2[hmax2]["net_Mm3"]):.0f} Mm³ 以下。'
         f'东半环 4,050 用 {H4050["fill_rcc_Mm3"]:.1f} Mm³ 换 {H4050["net_Mm3"]:.1f} Mm³ 天然库容（效率 {eff(H4050):.1f}），再挖 {H4050["deep_cut_Mm3"]:.0f} Mm³ 到 {H4050["deep_net_Mm3"]:.0f} Mm³——这是曲线上最后一个效率还在 1.5 以上的点。'
         f'拐点以上，从东半环到全环 B1 多花 {B1["fill_rcc_Mm3"]-H4050["fill_rcc_Mm3"]:.0f} Mm³ 填方只多 {B1["net_Mm3"]-H4050["net_Mm3"]:.0f} Mm³ 天然库容（效率 {(B1["net_Mm3"]-H4050["net_Mm3"])/(B1["fill_rcc_Mm3"]-H4050["fill_rcc_Mm3"]):.1f}），挖后多 {B1["deep_net_Mm3"]-H4050["deep_net_Mm3"]:.0f} Mm³ 也是靠多挖 {B1["deep_cut_Mm3"]-H4050["deep_cut_Mm3"]:.0f} Mm³ 换来的。'
         f'所以：要 5 Mm³ 以内的库容，S29 或 S2 一条 50–100 m 高的沟谷坝，填方 1–2 Mm³；要 20–40 Mm³，东半环 4,050，填方 10 Mm³ 加挖方 20 Mm³；再往上每 m³ 填方换不到 0.5 m³ 库容。</div>')
open(os.path.join(DEL, "D_height_series_brief.html"), "w").write(f'<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>坝高–填方–库容（简版）</title><style>{CSS}</style></head><body>{"".join(h)}</body></html>')
print("written")
