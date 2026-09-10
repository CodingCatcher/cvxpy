#!/usr/bin/env python3
"""交付物 1、2 的简版（团队对齐版）：只保留结论、每一步的手算公式与上下限、以及 2–3 张图。
输入：outputs/balanced_A_parcel_ring_3700.csv（碗形坑扫描）、outputs/types_A_parcel_ring_3700_stations.csv（逐站坡角）、
      outputs/D1_reaches.csv（分段）、outputs/D2v2_profile.png、DEM 与坝线。
输出：deliverables/D1_brief.html、deliverables/D2_brief.html、outputs/brief_sections.png、outputs/brief_checks.csv。
用法：python3 scripts/deliverables_1_2_brief.py"""
from __future__ import annotations
import base64, csv, json, math, os, sys
import numpy as np, shapely
from shapely.geometry import LineString, Point, Polygon
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
from terrain_model import Terrain, SECTIONS, FT, ACRE, DATA, OUT, Bowl, cross_section, dam_section, inward_normals, longitudinal_profile, load_parcel  # noqa
DEL = os.path.join(ROOT, "deliverables")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
matplotlib.rcParams["font.family"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]; matplotlib.rcParams["axes.unicode_minus"] = False
import logging; logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

LINE = "A_parcel_ring_3700"; LOWER_FT = 2400.0; ETA = 0.85; SF = 1.2; AGG = 0.9; BERM = 20.0; WALL_M = 0.75
COMBO_RCC = [(0.0, 606.0), (2302.8, 2695.2)]
COL = {"rcc": "#e69500", "cfrd": "#2e8b57"}
fnum = lambda v: float(v) if v not in ("", None) else float("nan")
kind_of = lambda ty, st: ("rcc" if any(a - 1e-6 <= st <= b + 1e-6 for a, b in COMBO_RCC) else "cfrd") if ty == "combo" else ty
gwh_per_Mm3 = lambda nwl_ft: 1e6 * 1000 * 9.81 * (nwl_ft - LOWER_FT) * FT * ETA / 3.6e12

t = Terrain(os.path.join(DATA, "dem_3dep_1m_utm10.tif"), os.path.join(DATA, "dem_3dep_1m_utm10.json"))
cand = json.load(open(os.path.join(DATA, "dam_line_candidates.json")))[LINE]
line = LineString(cand["coords"]); ring = Polygon(cand["coords"]); parcel = load_parcel()
s, pts, zg = longitudinal_profile(t, line, 10.0); nrm = inward_normals(pts, ring, None)
secs = [cross_section(t, pts[i], nrm[i]) for i in range(len(pts))]
win, m_ring = t.mask_inside(ring); z_ring = t.z[win]; X, Y = np.meshgrid(t.xc[win[1]], t.yc[win[0]])
srows = list(csv.DictReader(open(os.path.join(OUT, f"types_{LINE}_stations.csv"))))
b_ds = np.array([fnum(r["slope_ds_deg"]) for r in srows]); b_us = np.array([fnum(r["slope_us_deg"]) for r in srows])
reaches = list(csv.DictReader(open(os.path.join(OUT, "D1_reaches.csv"))))
bal = {(r["type"], round(float(r["crest_ft"]))): r for r in csv.DictReader(open(os.path.join(OUT, f"balanced_{LINE}.csv")))}

# ------------------------------------------------ 方案（设计点）
SCH = [  # key, 名称, 坝型, 坝顶, 取哪个库底
    ("A", "A 组合 III：1、4 段 RCC，其余 CFRD", "combo", 3900, "deep"),
    ("B1", "B1 全环 RCC（推荐设计点）", "rcc", 4050, "deep"),
    ("B2", "B2 全环 RCC（题设坝顶，极限）", "rcc", 4120, "bal"),
    ("C", "C 全环 CFRD", "cfrd", 3850, "deep"),
    ("R", "R 参照：全环 CFRD，固定坝顶 4,120", "cfrd", 4120, "deep"),
]

def build(ty, crest_ft):
    crest = crest_ft * FT; nwl = crest - 20 * FT; n = len(pts)
    A = np.full(n, np.nan); W = np.zeros(n); TU = np.full(n, np.nan); TD = np.full(n, np.nan); K = []; D = []
    for i in range(n):
        k = kind_of(ty, s[i]); K.append(k); d = dam_section(*secs[i], crest, SECTIONS[k], nwl); D.append(d)
        A[i] = d["area"]; W[i] = d["area_us_below_nwl"] if np.isfinite(d["area_us_below_nwl"]) else 0
        if d["toe_us"] is not None: TU[i] = d["toe_us"]
        if d["toe_ds"] is not None: TD[i] = d["toe_ds"]
    H = crest - zg
    toe_poly = Polygon([pts[i] + max(np.nan_to_num(TU[i], nan=0.0), 0.0) * nrm[i] for i in range(n)]).buffer(0)
    if toe_poly.geom_type != "Polygon": toe_poly = max(toe_poly.geoms, key=lambda g: g.area)
    return {"ty": ty, "crest_ft": crest_ft, "nwl": nwl, "A": A, "W": W, "TU": TU, "TD": TD, "K": np.array(K), "D": D, "H": H, "toe_poly": toe_poly, "bowl": Bowl(t, toe_poly, BERM, WALL_M)}

def excavate(b, floor_ft):
    inside, wmax = b["bowl"].field(X, Y, z_ring); Fm = floor_ft * FT
    z_new = np.where(inside, np.minimum(z_ring, np.maximum(Fm, wmax)), z_ring); dz = np.maximum(z_ring - z_new, 0)
    cut = float(np.nansum(dz)); cut_below = float(np.nansum(np.maximum(np.minimum(z_ring, b["nwl"]) - z_new, 0)))
    gross = float(np.nansum(np.where(m_ring, np.maximum(b["nwl"] - z_new, 0), 0))); water = m_ring & (z_new < b["nwl"])
    top = b["bowl"].top; zin = z_ring[inside]
    gross0 = float(np.nansum(np.where(m_ring, np.maximum(b["nwl"] - z_ring, 0), 0)))   # 真正不挖时的毛库容
    return {"cut": cut, "cut_below_nwl": cut_below, "gross": gross, "gross_nocut": gross0, "water_acre": float(water.sum() * t.px ** 2 / ACRE),
            "top_acre": top.area / ACRE, "top_perim_m": top.length, "mean_ground_top_ft": float(np.nanmean(zin)) / FT, "max_depth_m": float(np.nanmax(dz))}

B = {}; EXC = {}; ROWS = []
for key, name, ty, crest, which in SCH:
    r = bal[(ty, crest)]; b = build(ty, crest); B[key] = b
    floor = fnum(r["floor_ft"]) if which == "bal" else fnum(r["maxnet_floor_ft"])
    e = excavate(b, floor); EXC[key] = e
    grid = fnum(r["dam_grid_Mm3"]); gc = fnum(r.get("grid_cfrd_Mm3", "")); gr = fnum(r.get("grid_rcc_Mm3", ""))
    if not np.isfinite(gc): gc, gr = (grid, 0) if ty == "cfrd" else (0, grid)
    need = fnum(r["need_bank_Mm3"]); cut = fnum(r["cut_Mm3"]) if which == "bal" else fnum(r["maxnet_cut_Mm3"])
    net = fnum(r["net_Mm3"]) if which == "bal" else fnum(r["maxnet_Mm3"]); g = fnum(r["net_GWh"]) if which == "bal" else fnum(r["maxnet_GWh"])
    ROWS.append({"key": key, "name": name, "ty": ty, "crest": crest, "nwl": crest - 20, "H": fnum(r["H_max_m"]), "dam_len_km": float((b["H"] > 0).sum() * 10) / 1000,
                 "grid": grid, "gc": gc, "gr": gr, "need": need, "need_hand": gc / SF + gr * AGG, "maxcut": fnum(r["max_cut_Mm3"]), "floor": floor, "cut": cut, "ratio": cut / need if need else np.nan,
                 "nocut": e["gross_nocut"] / 1e6 - fnum(r["wedge_Mm3"]), "wedge": fnum(r["wedge_Mm3"]), "net": net, "gwh": g, "gwh_hand": net * gwh_per_Mm3(crest - 20),
                 "cut_below": e["cut_below_nwl"] / 1e6, "identity": e["gross_nocut"] / 1e6 - fnum(r["wedge_Mm3"]) + e["cut_below_nwl"] / 1e6,  # 不挖净库容（已扣楔）+ 水位以下挖方
                 "gross": e["gross"] / 1e6, "water_acre": e["water_acre"], "mean_depth_m": e["gross"] / (e["water_acre"] * ACRE) if e["water_acre"] else np.nan,
                 "top_acre": e["top_acre"], "perim_m": e["top_perim_m"], "mean_ground_top_ft": e["mean_ground_top_ft"], "max_depth_m": e["max_depth_m"],
                 "pit_upper": e["top_acre"] * ACRE * max(0, (e["mean_ground_top_ft"] - floor) * FT) / 1e6,
                 "pit_lower": max(0, e["top_acre"] * ACRE * max(0, (e["mean_ground_top_ft"] - floor) * FT) - e["top_perim_m"] * 0.5 * WALL_M * (max(0, (e["mean_ground_top_ft"] - floor) * FT)) ** 2) / 1e6,
                 "balanced": r["balanced"] == "True"})
    print(f"{key}: cut {e['cut']/1e6:.1f} (csv {cut:.1f}) below-nwl {e['cut_below_nwl']/1e6:.1f} gross {e['gross']/1e6:.1f} net-identity {ROWS[-1]['identity']:.1f} vs net {net:.1f} | pit bound {ROWS[-1]['pit_lower']:.0f}–{ROWS[-1]['pit_upper']:.0f} vs cut {cut:.1f}")

# ------------------------------------------------ 分段手算（方案 A 的坝顶）
bA = B["A"]; crestA = SCH[0][3]
def reach_check(b, ty_full=None):
    out = []
    for r in reaches:
        a, c = fnum(r["sta_from"]), fnum(r["sta_to"]); idx = np.where((s >= a - 1e-6) & (s <= c + 1e-6))[0]
        H = b["H"][idx]; dam = H > 0
        if not dam.any(): out.append({"reach": int(float(r["reach"])), "L": len(idx) * 10, "dam": False}); continue
        k = ty_full or b["K"][idx][dam][0]; spec = SECTIONS[k]; m = spec["m_ds"]; mu = spec["m_us"]; bw = spec["crest_w"]
        Hm = float(H[dam].mean()); bd = float(np.nanmean(b_ds[idx][dam])); bu = float(np.nanmean(b_us[idx][dam])); L = float(dam.sum() * 10)
        td = Hm / (1 / m - math.tan(math.radians(bd))) if math.tan(math.radians(bd)) < 1 / m else float("inf")
        tu = Hm / (1 / mu + math.tan(math.radians(bu))) if mu > 0 else 0.0
        A_est = 0.5 * Hm * (td + tu) + bw * Hm if np.isfinite(td) else float("inf")
        if ty_full:
            AA = np.array([dam_section(*secs[i], b["crest_ft"] * FT, spec, b["nwl"])["area"] for i in idx]); TDm = np.array([dam_section(*secs[i], b["crest_ft"] * FT, spec, b["nwl"])["toe_ds"] or np.nan for i in idx])
        else:
            AA = b["A"][idx]; TDm = b["TD"][idx]
        V_model = float(np.nansum([0.5 * (AA[j] + AA[j + 1]) * 10 for j in range(len(AA) - 1) if np.isfinite(AA[j]) and np.isfinite(AA[j + 1])])) / 1e6
        out.append({"reach": int(float(r["reach"])), "L": L, "dam": True, "kind": k, "Hm": Hm, "Hmax": float(H.max()), "bd": bd, "bu": bu, "td_est": td, "td_model": float(-np.nanmin(TDm)) if np.isfinite(TDm).any() else np.nan,
                    "tu_est": tu, "A_est": A_est, "A_model": float(np.nanmean(AA[dam])), "V_est": A_est * L / 1e6 if np.isfinite(A_est) else float("inf"), "V_model": V_model})
    return out
RC_A = reach_check(bA); RC_Acfrd = reach_check(bA, "cfrd"); RC_Arcc = reach_check(bA, "rcc")

# ------------------------------------------------ 图：4 个关键断面
def draw_section(ax, b, i, kind, floor_ft, title):
    off, zsec = secs[i]; zft = zsec / FT; crest_ft = b["crest_ft"]; nwl_ft = crest_ft - 20
    d = dam_section(off, zsec, crest_ft * FT, SECTIONS[kind], b["nwl"]); H = b["H"][i]
    P = pts[i][None, :] + off[:, None] * nrm[i][None, :]; znew = b["bowl"].surface(P, zsec, floor_ft * FT); znew = np.where(off >= 0, znew, zsec)
    ax.fill_between(off, znew / FT, zft, where=znew < zsec - 1e-6, color="#f4a582", alpha=0.75, lw=0, label=f"开挖到库底 {floor_ft:.0f} ft")
    ax.plot(off, zft, "k-", lw=1.2, label="地面"); ax.plot(off, znew / FT, "-", color="#b2182b", lw=0.7)
    bw = SECTIONS[kind]["crest_w"]
    if d["toe_us"] is not None and d["toe_ds"] is not None:
        xs = [d["toe_ds"], -bw / 2, bw / 2, d["toe_us"]]; zs = [d["z_toe_ds"] / FT, crest_ft, crest_ft, d["z_toe_us"] / FT]
        m = (off >= d["toe_ds"]) & (off <= d["toe_us"])
        ax.fill(list(xs) + list(off[m][::-1]), list(zs) + list(zft[m][::-1]), color=COL[kind], alpha=0.65, lw=0.8, edgecolor="k", label=f"{kind.upper()}：底宽 {d['base_width']:.0f} m，面积 {d['area']:.0f} m²")
        j = int(np.argmin(np.abs(off - bw / 2))); k2 = j; wz = znew / FT
        while k2 < len(off) and (wz[k2] < nwl_ft or off[k2] <= d["toe_us"]): k2 += 1
        ax.fill_between(off[j:k2], np.minimum(wz[j:k2], nwl_ft), nwl_ft, color="#2c7fb8", alpha=0.25, lw=0, label="水")
    else:
        mds = SECTIONS[kind]["m_ds"]; xx = np.linspace(-900, -bw / 2, 50); ax.plot(xx, crest_ft - (-xx - bw / 2) / mds / FT, "--", color=COL[kind], lw=1.5, label=f"{kind.upper()} 下游面：900 m 内碰不到地面（不闭合）")
    out_ds = None
    for jj in range(len(off)):
        if off[jj] < 0 and not parcel.contains(Point(P[jj])): out_ds = off[jj]
    if out_ds is not None: ax.axvline(out_ds, color="k", ls="-.", lw=0.8); ax.text(out_ds, crest_ft + 30, "宗地界", fontsize=8, ha="center")
    ax.axhline(nwl_ft, color="b", ls="--", lw=0.7); ax.axhline(crest_ft, color="r", lw=0.6)
    lo = np.nanmin(np.where(np.isfinite(znew), znew, np.inf)) / FT
    ax.set_xlim(-900, 600); ax.set_ylim(min(lo, crest_ft - 650) - 30, crest_ft + 80); ax.grid(alpha=0.3); ax.legend(loc="lower right", fontsize=7); ax.set_title(title, fontsize=9)
    ax.set_xlabel("离坝轴偏移 m（+ 为库侧）"); ax.set_ylabel("ft")
    return d

iS2 = int(np.argmin(np.abs(s - 1298))); iS4 = int(np.argmin(np.abs(s - 2432)))
fig, axs = plt.subplots(2, 2, figsize=(16, 8.4))
tb = lambda i: math.tan(math.radians(b_ds[i]))
d1 = draw_section(axs[0, 0], bA, iS2, "cfrd", EXC["A"] and fnum(bal[("combo", 3900)]["maxnet_floor_ft"]), f"S2 桩号 1,298（第 2 段）方案 A：CFRD，H {bA['H'][iS2]:.0f} m，下游坡 {b_ds[iS2]:.0f}°  ——  tanβ={tb(iS2):.2f} < 1/1.4=0.71，闭合；手算坝趾 H/(0.71−tanβ)={bA['H'][iS2]/(1/1.4-tb(iS2)):.0f} m")
d2 = draw_section(axs[0, 1], bA, iS4, "cfrd", fnum(bal[("combo", 3900)]["maxnet_floor_ft"]), f"S4 桩号 2,432（第 4 段）若用 CFRD：下游坡 {b_ds[iS4]:.0f}°  ——  tanβ={tb(iS4):.2f} > 0.71，公式说落不了地；模型：600 m 外坡度变缓才落地，已出宗地界 200 m")
d3 = draw_section(axs[1, 0], bA, iS4, "rcc", fnum(bal[("combo", 3900)]["maxnet_floor_ft"]), f"S4 同一位置方案 A：RCC，H {bA['H'][iS4]:.0f} m  ——  1/0.8=1.25 > tanβ={tb(iS4):.2f}，闭合；手算坝趾 {bA['H'][iS4]/(1/0.8-tb(iS4)):.0f} m")
d4 = draw_section(axs[1, 1], B["B1"], iS2, "rcc", fnum(bal[("rcc", 4050)]["maxnet_floor_ft"]), f"S2 方案 B1：全环 RCC，坝顶 4,050，H {B['B1']['H'][iS2]:.0f} m，库底 {fnum(bal[('rcc', 4050)]['maxnet_floor_ft']):.0f} ft")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "brief_sections.png"), dpi=105); plt.close(fig)
SEC4 = [("S2 A CFRD", bA["H"][iS2], b_ds[iS2], 1.4, d1), ("S4 A 若 CFRD", bA["H"][iS4], b_ds[iS4], 1.4, d2), ("S4 A RCC", bA["H"][iS4], b_ds[iS4], 0.8, d3), ("S2 B1 RCC", B["B1"]["H"][iS2], b_ds[iS2], 0.8, d4)]

# ------------------------------------------------ 输出 CSV
with open(os.path.join(OUT, "brief_checks.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(ROWS[0].keys())); w.writeheader(); w.writerows(ROWS)

# ------------------------------------------------ HTML
img64 = lambda p: "data:image/png;base64," + base64.b64encode(open(p, "rb").read()).decode()
esc = lambda x: str(x).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
def fig_html(p, cap): return f'<figure><img src="{img64(p)}"><figcaption>{cap}</figcaption></figure>'
def table(cols, rows):
    h = ['<div class="wrap"><table><tr>' + "".join(f"<th>{esc(c)}</th>" for c in cols) + "</tr>"]
    for r in rows: h.append("<tr>" + "".join(f"<td>{esc(x)}</td>" for x in r) + "</tr>")
    return "\n".join(h) + "</table></div>"
CSS = """body{font-family:"Noto Sans CJK SC","WenQuanYi Zen Hei","PingFang SC","Microsoft YaHei",Arial,sans-serif;max-width:1150px;margin:0 auto;padding:18px 26px 50px;line-height:1.6;color:#1f2328;background:#fff}
h1{font-size:21px;border-bottom:2px solid #333;padding-bottom:6px} h2{font-size:17px;margin-top:28px;border-left:5px solid #7b1fa2;padding-left:10px}
table{border-collapse:collapse;font-size:13px;margin:8px 0 12px} th,td{border:1px solid #bfc4cc;padding:4px 8px;vertical-align:top;text-align:left} th{background:#eef1f5}
tr:nth-child(even) td{background:#f8f9fb} .wrap{overflow-x:auto} figure{margin:12px 0} img{max-width:100%;border:1px solid #d8dbe0} figcaption{font-size:12.5px;color:#444}
.box{background:#f3e5f5;border-left:4px solid #7b1fa2;padding:8px 14px;margin:10px 0} .warn{background:#fff8e1;border-left:4px solid #f0b400;padding:8px 14px;margin:10px 0}
.f{font-family:Menlo,Consolas,monospace;background:#f1f3f5;padding:1px 5px;border-radius:3px;font-size:12.5px} .small{font-size:12px;color:#555} li{margin:4px 0}"""
doc = lambda title, body: f'<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{esc(title)}</title><style>{CSS}</style></head><body>{body}</body></html>'
R = {r["key"]: r for r in ROWS}; A, B1, B2, C, RR = R["A"], R["B1"], R["B2"], R["C"], R["R"]
f0 = lambda v: f"{v:,.0f}"; f1 = lambda v: f"{v:.1f}"
ring_acre = ring.area / ACRE; ring_km = line.length / 1000; zmin = float(np.nanmin(zg)) / FT

# ---------- D1 简版
h = []
h.append("<h1>交付物 1（简版）：坝型比较——先对齐框架、公式和结论</h1>")
h.append('<div class="box"><b>一句话框架。</b>上库是在天然地里打出来的孔，所有坝料来自这个孔；坝线固定在宗地内沿 3,700 ft 等高线闭合的环，<b>坝顶高程和库底高程是变量</b>，问"这个孔能养多高的坝、库容多大"。'
         '坝型只在两种能在陡坡上落地的断面里选：面板堆石 CFRD（上下游 1.4:1，顶宽 10 m）与碾压混凝土 RCC（上游直立、下游 0.8:1，顶宽 6.1 m）。</div>')
h.append("<h2>1. 地形事实（三条，可在图上核对）</h2>")
h.append(f"<ul><li>环长 {ring_km:.1f} km，环内 {ring_acre:.0f} acre（{ring_acre*ACRE/1e6:.2f} km²）。西、北、东北三面沿 3,700 线（最低坝基 {zmin:,.0f} ft），东界与南界的地面本来就有 3,870–4,230 ft。</li>"
         f"<li>下游山坡：西坡约 30°，北鼻尖 18–21°，东北坡 31–42°，东界与南台地 12–22°。</li>"
         f"<li>坝顶 3,900 时东界与南台地不需要坝（地面高于坝顶），坝长 {A['dam_len_km']:.1f} km；坝顶 4,120 时坝绕整个环，{B2['dam_len_km']:.1f} km。</li></ul>")
h.append("<h2>2. 五个手算公式（每个数字都能这样核）</h2>")
h.append('<ol>'
         '<li><b>坝脚落地。</b>坝面坡 m:1、地面坡角 β，下游坝脚离坝轴 <span class="f">d = H / (1/m − tanβ)</span>；tanβ ≥ 1/m 时永不落地。CFRD 1.4:1 的极限是 35.5°，RCC 0.8:1 是 51°。库侧地面向上，<span class="f">d_上 = H / (1/m + tanβ_上)</span>。</li>'
         '<li><b>断面与方量。</b><span class="f">A ≈ ½·H·(d_下 + d_上) + 顶宽·H</span>，<span class="f">V ≈ A × 坝长</span>。平地上 CFRD 的 A ≈ 1.4 H²，RCC ≈ 0.4 H²，所以同高度 RCC 约为堆石的 1/4–1/3；坡地上按 d 放大。</li>'
         '<li><b>需料。</b><span class="f">原岩 = 堆石压实方 ÷ 1.2 + RCC 方 × 0.9</span>（1 m³ 原岩压实成 1.2 m³ 堆石；每 m³ RCC 需 0.9 m³ 骨料原岩）。</li>'
         '<li><b>孔能出多少。</b>上限 <span class="f">坑顶面积 × (坑内平均地面 − 库底)</span>；坑壁 0.75:1 占去的楔 ≈ <span class="f">坑顶周长 × 0.375 × 深²</span>，扣掉得下限。坑顶线在上游坝脚以内 20 m。</li>'
         '<li><b>库容与电量。</b><span class="f">净库容 = 不挖时的库容 + 水位以下的挖方 − 堆石上游坡占掉的库容</span>；也可按 <span class="f">水面面积 × 平均水深</span> 核。<span class="f">GWh ≈ 净库容(Mm³) × 水头(m) × 0.00232</span>（η 0.85，下库 2,400 ft）。</li></ol>')
h.append(fig_html(os.path.join(OUT, "D1v2_plan.png"), "图 1 方案 A（左，坝顶 3,900）与全环 RCC 4,120（右）的平面：橙/绿是坝体足迹（坝脚到坝脚），黄–红是坑的开挖深度，红线是坑顶线，蓝是水面。左图东界与南台地没有坝。"))
h.append("<h2>3. 方案比较（每个方案一个设计点）</h2>")
rows = []
for r in ROWS:
    rows.append([r["name"], f"{r['crest']:,} / {r['nwl']:,}", f0(r["H"]), f1(r["dam_len_km"]), f"{r['grid']:.1f}" + (f"（堆石 {r['gc']:.1f} + RCC {r['gr']:.1f}）" if r["gc"] > 0 and r["gr"] > 0 else ""),
                 f"{r['need']:.1f}（手算 {r['need_hand']:.1f}）", f"{r['maxcut']:.1f}", f"{r['floor']:,.0f}", f"{r['cut']:.1f}", f"{r['ratio']:.2f}" + ("" if r["balanced"] else "（料不够）"), f"{r['net']:.1f}", f"{r['gwh']:.0f}（手算 {r['gwh_hand']:.0f}）"])
h.append(table(["方案", "坝顶 / 水位 ft", "H 最大 m", "坝长 km", "坝体 Mm³（网格法）", "需原岩 Mm³", "孔的极限 Mm³", "库底 ft", "挖方 Mm³", "挖 ÷ 需", "净库容 Mm³", "GWh"], rows))
h.append('<div class="warn"><b>怎么读。</b>"孔的极限"是坑壁 0.75:1 挖到底面消失时的挖方；坝越高、堆石上游坡越长，坑顶线越靠里，极限越小，而需料按 H² 涨——两条线交叉处就是各坝型的坝顶上限：'
         f'全 CFRD {C["crest"]:,}，组合 III {A["crest"]:,}，全 RCC 到 4,120 只剩 {B2["ratio"]:.2f} 倍，4,050 还有 {B1["ratio"]:.2f} 倍。参照 R 就是笔记 13 说的"缺口"：固定坝顶 4,120 的全堆石要 {RR["need"]:.0f} Mm³，孔只有 {RR["maxcut"]:.0f} Mm³。</div>')
h.append(fig_html(os.path.join(OUT, "D1v2_need_vs_cut.png"), "图 2 需料（实线，随 H² 上升）与孔的极限挖方（虚线，随坝顶上升而下降）：实线在虚线之下才养得起。点线是初版把坑算大了的结果，仅作对照。"))
h.append("<h2>4. 三个数字的手算核对（方案 A 与 B1）</h2>")
chk = []
for r in (A, B1):
    chk.append([r["name"], f"{r['nocut']:.1f} + {r['cut_below']:.1f} = {r['identity']:.1f}，模型 {r['net']:.1f}", f"{r['water_acre']:.0f} acre × {r['mean_depth_m']:.0f} m = {r['gross']:.1f} 毛库容，减上游坡 {r['wedge']:.1f} → {r['gross']-r['wedge']:.1f}",
                f"坑顶 {r['top_acre']:.0f} acre，平均地面 {r['mean_ground_top_ft']:,.0f} ft，库底 {r['floor']:,.0f}：上限 {r['pit_upper']:.0f}，扣坑壁楔 → {r['pit_lower']:.0f}；模型 {r['cut']:.1f}（最深 {r['max_depth_m']:.0f} m）"])
h.append(table(["方案", "净库容 = 不挖库容 + 水位以下挖方（Mm³）", "面积 × 平均深", "挖方上下限（Mm³）"], chk))
h.append('<p class="small">"不挖库容"已扣上游坡；恒等式核法与模型只差舍入（水位以上的挖方——台地脚下那部分——不产生库容，所以"水位以下挖方"小于总挖方）。挖方的上下限是把坑当平地、坑壁当直墙的粗估，模型按真实地面逐米算。</p>')
h.append("<h2>5. 每一段用什么（坝顶 3,900，方案 A）</h2>")
seg = []
why = {1: "西坡 30°：堆石做得成但底宽 570 m、一段吃掉全环堆石料的 27%；RCC 底宽 118 m", 2: "北鼻尖 21°，堆石坝脚 270 m 在界内；这一段是孔料的主要去处", 3: "18° 但坡脚有台阶，堆石坝脚 390 m；坝顶想再高 50 ft 就先把这段改 RCC",
       4: "东北坡 31–42°：堆石坝脚伸出宗地北界 200 m，任何坝顶都如此；RCC 坝脚 164 m 在界内", 5: "沿东界的过渡段，两种坝型的坝脚都出界 40–65 m，坝轴要内移", 6: "开头 60 m 地面比坝顶低不到 6 m，一段小坝；其余天然山体",
       12: "沿西界顺坡而下的闭合段，坝脚必然出界 44–73 m：买一条带或改线"}
kinds = {1: "RCC", 2: "CFRD", 3: "CFRD", 4: "RCC", 5: "CFRD（坝轴内移）", 6: "小坝 60 m + 天然山体", 12: "RCC 或改线"}
for rc, rcf, rcr in zip(RC_A, RC_Acfrd, RC_Arcc):
    k = rc["reach"]
    if k in (7, 8, 9, 10, 11): continue
    seg.append([k, f"{rc['L']:.0f}", f"{rc['Hm']:.0f}", f"{rc['bd']:.0f}", kinds[k],
                (f"{rcf['td_est']:.0f}" if np.isfinite(rcf["td_est"]) else "不落地") + f" / {rcf['td_model']:.0f}", f"{rcr['td_est']:.0f} / {rcr['td_model']:.0f}",
                (f"{rcf['V_est']:.1f}" if np.isfinite(rcf["V_est"]) else "—") + f" / {rcf['V_model']:.1f}", f"{rcr['V_est']:.1f} / {rcr['V_model']:.1f}", why[k]])
seg.append(["7–11", f"{sum(rc['L'] for rc in RC_A if rc['reach'] in (7, 8, 9, 10, 11)):,.0f}（无坝）", "地面高于坝顶", "—", "无坝", "—", "—", "—", "—", "东界后段与南台地：地面 3,950–4,230 ft，天然挡水"])
h.append(table(["段", "有坝长 m", "平均 H m", "下游坡°", "选型", "CFRD 坝脚：公式 / 模型 m", "RCC 坝脚：公式 / 模型 m", "CFRD 方量：公式 / 模型 Mm³", "RCC 方量：公式 / 模型 Mm³", "一句理由"], seg))
h.append('<p class="small">公式用段内平均 H 与 100 m 窗口的平均坡角，模型逐站按真实断面找坝脚；两者差 20–40% 是局部地形（坡脚台阶、冲沟）造成的，趋势一致。段号与桩号见交付物 2 的纵剖面。</p>')
h.append("<h2>6. 结论与推荐</h2>")
h.append(f'<ul><li><b>目标"上面尽可能高"</b> → 全环 RCC，坝顶定 4,000–4,050（B1：净 {B1["net"]:.0f} Mm³ / {B1["gwh"]:.0f} GWh，料 {B1["ratio"]:.2f} 倍，坑最深 {B1["max_depth_m"]:.0f} m）。4,120 是极限（B2：料 {B2["ratio"]:.2f} 倍，要挖到 {B2["floor"]:,.0f} ft，坑深 {B2["max_depth_m"]:.0f} m）。代价：约 30–42 Mm³ RCC 的胶材外运、拌合站十年以上、129 m 直立面的抗滑。</li>'
         f'<li><b>目标"挖出来的石料堆成坝、少用水泥"</b> → 组合 III，坝顶 3,900（A：净 {A["net"]:.0f} Mm³ / {A["gwh"]:.0f} GWh，料 {A["ratio"]:.2f} 倍）。库容是 B1 的四成，但 RCC 只有 5 Mm³。</li>'
         '<li><b>两个方案都成立的事</b>：西坡与东北坡（第 1、4 段）必须 RCC；东界过渡段坝轴内移；西南闭合段要用地或改线。</li>'
         '<li><b>先要确认的三个假设</b>：千枚岩碎料压实 1.2 倍与台地泥流盖层能否作坝料（试验采场）；坑壁 0.75:1（完好岩 0.5:1 可多挖一到两成，组合 III 或可回到 3,950）；下库位置与水头（电量按 2,400 ft 河谷算）。</li></ul>')
h.append('<p class="small">细节版：deliverables/D1_dam_type_comparison.html；复算：python3 scripts/deliverables_1_2_brief.py。</p>')
open(os.path.join(DEL, "D1_brief.html"), "w").write(doc("交付物 1（简版）", "\n".join(h)))

# ---------- D2 简版
g = []
g.append("<h1>交付物 2（简版）：整条坝线的纵剖面与四个关键断面</h1>")
g.append('<div class="box">桩号从南界西端起顺时针：西坡（第 1 段）→ 北鼻尖（2、3）→ 东北坡（4）→ 东界（5、6）→ 南台地（7–11）→ 西南闭合（12）。看图只需要记三件事：<b>地面线在哪</b>（黑线）、<b>坝顶线在哪</b>（两条水平线），<b>两者之间就是坝</b>；坝顶线以上的地面段不需要坝。</div>')
g.append("<h2>1. 纵剖面</h2>")
g.append(fig_html(os.path.join(OUT, "D2v2_profile.png"), "上图：黑线是坝轴上的地面；紫线坝顶 3,900（方案 A），红线 4,120（方案 B2）；橙/绿是方案 A 的 RCC/CFRD 坝体，点状阴影是坝顶抬到 4,120 再加的 RCC；绿点/红叉是各段下游坝脚落地的高程，点越低说明坝脚伸到越深的坡上；下图是 H 与地面坡角，虚线 35.5° 是 CFRD 能落地的极限坡角。"))
g.append("<h2>2. 纵剖面上能直接手算核对的事</h2>")
g.append('<ul><li><b>H 就是坝顶减地面。</b>第 1–4 段地面 3,695–3,705 ft，坝顶 3,900 → H = 195 ft ≈ 60 m；坝顶 4,120 → 127 m。第 6–11 段地面 3,870–4,230，坝顶 3,900 时 H ≤ 6 m 或为负（无坝）。</li>'
         '<li><b>坝脚多远看坡角。</b>第 4 段下游坡角 31–42°，tanβ = 0.6–0.9，超过 CFRD 的 0.71，所以图上第 4 段绿点消失、只有红叉（RCC 的 1.25 还够）。第 2 段 21°，tanβ = 0.38，堆石坝脚 d = 61/(0.71 − 0.38) ≈ 185 m（模型 272 m，坡脚有台阶）。</li>'
         '<li><b>库底线。</b>方案 A 库底 3,500 ft 比最低坝基低 60 m；B1 3,450；B2 3,370。坑壁 0.75:1，从坝脚内 20 m 下到库底要 45–75 m 的水平距离。</li></ul>')
g.append("<h2>3. 四个关键断面</h2>")
g.append(fig_html(os.path.join(OUT, "brief_sections.png"), "左上：第 2 段用 CFRD（方案 A），坝脚在界内，上游坡脚以内是坑。右上：第 4 段若用 CFRD，1.4:1 的坝面比 42° 的地面还缓，900 m 内碰不到地面（虚线），实际在 100 m 窗口外坡度变缓才落地，但已出宗地北界。左下：同一位置用 RCC，坝脚 164 m 落地、在界内。右下：第 2 段在方案 B1（全 RCC 4,050）下的样子，坑更深。橙红是开挖，蓝是水，点划线是宗地界。"))
rows4 = []
for name, H, beta, m, d in SEC4:
    tb_ = math.tan(math.radians(beta)); td = H / (1 / m - tb_) if tb_ < 1 / m else float("inf")
    rows4.append([name, f"{H:.0f}", f"{beta:.0f}（tanβ {tb_:.2f}）", f"{m}:1（1/m {1/m:.2f}）", "不落地" if not np.isfinite(td) else f"{td:.0f}", "不落地" if d["toe_ds"] is None else f"{-d['toe_ds']:.0f}",
                  "—" if d["toe_ds"] is None else f"{d['base_width']:.0f}", "—" if d["toe_ds"] is None else f"{d['area']:.0f}", f"{0.5*H*H*(1/(1/m-tb_) if tb_<1/m else float('nan')):.0f}" if tb_ < 1 / m else "—"])
g.append(table(["断面", "H m", "下游坡角°", "坝面坡", "坝脚：公式 H/(1/m−tanβ) m", "坝脚：模型 m", "底宽 m", "断面面积 m²", "下游三角粗估 ½H·d m²"], rows4))
g.append('<p class="small">"下游三角粗估"只算下游侧的三角形，没算库侧和顶宽，所以比模型面积小；它说明的是量级和"落不落地"。第 2 段公式给 185 m、模型 272 m，是因为坡脚 200 m 外有一段台阶。</p>')
g.append("<h2>4. 图直接说明的三件事</h2>")
g.append(f'<ul><li>坝顶从 4,120 降到 3,900，坝从 {B2["dam_len_km"]:.1f} km 缩到 {A["dam_len_km"]:.1f} km：东界后段与南台地本来就高于水位，是天然的库岸。</li>'
         '<li>坝型由坡角决定：第 1、4 段（30° 与 31–42°）只有 RCC 能在宗地内落地；第 2、3、5 段（12–21°）堆石能落地，坑就在坝脚以内。</li>'
         '<li>坑在断面上是一个 0.75:1 的碗，深 60–170 m；库侧地面越高（第 3、5 段），坑壁占去的带越宽，这也是初版把坑算大了的原因。</li></ul>')
g.append('<p class="small">细节版（14 个断面、断面表）：deliverables/D2_profile_and_sections.html。</p>')
open(os.path.join(DEL, "D2_brief.html"), "w").write(doc("交付物 2（简版）", "\n".join(g)))
print("written", [os.path.getsize(os.path.join(DEL, f)) // 1024 for f in ("D1_brief.html", "D2_brief.html")], "KB")
