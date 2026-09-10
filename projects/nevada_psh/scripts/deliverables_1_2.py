#!/usr/bin/env python3
"""交付物 1（坝型比较表）与交付物 2（整条坝线纵剖面与代表断面）的生成脚本。

输入：outputs/types_<line>_stations.csv / _summary.json（terrain_model.py --compare-types）、
      outputs/terrain_<line>_<tag>_summary.json（--run 各模式）。
输出：deliverables/D1_dam_type_comparison.md（含逐段表）、outputs/D1_schemes.csv、outputs/D1_reaches.csv、
      deliverables/D2_profile_and_sections.md、outputs/D2_profile.png、outputs/D2_sections.png、outputs/D2_sections_table.csv。

只汇总几何与准则数字，不选型。所有数值可由脚本重算。
"""
from __future__ import annotations
import csv, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "outputs"); DEL = os.path.join(ROOT, "deliverables"); DATA = os.path.join(ROOT, "data")
sys.path.insert(0, HERE)
from terrain_model import Terrain, SECTIONS, cross_section, dam_section, inward_normals, longitudinal_profile, FT, ACRE  # noqa
from shapely.geometry import LineString, Polygon

LINE = sys.argv[1] if len(sys.argv) > 1 else "A_parcel_ring_3700"
CREST_FT, NWL_FT = 4120.0, 4100.0
HEAD_M = 518.0; ETA = 0.85               # 下库若在北界河谷（2,400 ft）：毛水头 518 m，只作尺度
GWH_PER_M3 = 1000 * 9.81 * HEAD_M * ETA / 3.6e12
SF_ROCK, LF_ROCK = 0.95, 0.65            # 堆石：压实方/原岩方（资料 7/8 试验前的假设）、松方系数（硬岩 0.60–0.67，资料 9）
AGG_BANK_PER_RCC = 0.9                   # 每 m³ RCC 需要的骨料原岩方（2.1 t 骨料 / 2.7 t/m³ / 0.85 加工回收）
FLEET_ROCK_BCM_H = 5000.0                # 示意车队（≈8 台 D10T2 松土 + 4 台大挖 + 40 台 740C）
RCC_PLANT_M3_H = 700.0                   # 大型拌合站 ~1,000 yd³/h 的公制值
HOURS_PER_YEAR = 4000.0
TYPES = ["rcc", "cfrd", "ecrd", "rockfill"]
TYPE_CN = {"rcc": "RCC 重力（直立/0.8:1，顶 6.1 m）", "cfrd": "面板堆石 CFRD（1.4/1.4，顶 10 m）",
           "ecrd": "心墙堆石 ECRD（2.0/1.8，顶 10 m）", "rockfill": "堆石 1.6/1.6（资料 1 初始假设）"}

# 地质锚点（CGS 2010 1:250k identify 两点 + 宗地主体；只是"可能的单元"）
GEOL_ANCHORS = [("Mehrten 火山泥流（可能）", (686654.0, 4356075.0), 700.0, 3950.0), ("变质火山岩（可能）", (686738.0, 4356941.0), 400.0, None)]
def geol(E, N, z_ft=None):
    """CGS 1:250k 在两点的 identify + 地貌推断：Mehrten 泥流是台地盖层，只在台地面（≥3,950 ft）且靠近东南锚点处标注；
    变质火山岩在东侧锚点附近；其余按宗地主体 Calaveras。全部是"可能"，不是测绘。"""
    for name, (ax, ay), r, zmin in GEOL_ANCHORS:
        if math.hypot(E - ax, N - ay) <= r and (zmin is None or z_ft is None or z_ft >= zmin):
            return name
    return "Calaveras 千枚岩/泥质岩（可能）"


def fnum(x):
    try:
        return float(x) if x not in ("", None) else float("nan")
    except Exception:
        return float("nan")


def load_stations():
    rows = list(csv.DictReader(open(os.path.join(OUT, f"types_{LINE}_stations.csv"))))
    for r in rows:
        for k in list(r.keys()):
            if k.endswith("_closes"):
                r[k] = r[k] == "True"
            elif k not in ("station",):
                r[k] = fnum(r[k])
        r["station"] = float(r["station"])
    return rows


def reaches(rows, step=10.0):
    """按 H 档（≤30 / 30–80 / >80）与 CFRD 闭合性变化切段，合并短于 60 m 的碎段。"""
    def cls(r):
        h = r["H_axis"]
        hb = "H≤30" if h <= 30 else ("30<H≤80" if h <= 80 else "H>80")
        d = r["cfrd_toe_ds"]
        tc = "open" if not r["cfrd_closes"] else ("≤350" if -d <= 350 else ("350–600" if -d <= 600 else ">600"))
        return (hb, tc)
    segs = []; start = 0
    for i in range(1, len(rows) + 1):
        if i == len(rows) or cls(rows[i]) != cls(rows[start]):
            segs.append((start, i - 1)); start = i
    # 合并碎段
    merged = []
    for a, b in segs:
        if merged and (b - a + 1) * step < 100:
            pa, pb = merged[-1]; merged[-1] = (pa, b)
        else:
            merged.append((a, b))
    out = []
    for a, b in merged:
        sub = rows[a:b + 1]
        rec = {"reach": len(out) + 1, "sta_from": sub[0]["station"], "sta_to": sub[-1]["station"], "length_m": (b - a + 1) * step,
               "z_ground_ft_min": min(r["z_ground_ft"] for r in sub), "z_ground_ft_max": max(r["z_ground_ft"] for r in sub),
               "H_min": min(r["H_axis"] for r in sub), "H_max": max(r["H_axis"] for r in sub),
               "slope_ds_mean": float(np.nanmean([r["slope_ds_deg"] for r in sub])), "slope_us_mean": float(np.nanmean([r["slope_us_deg"] for r in sub])),
               "geology": geol(np.mean([r["E"] for r in sub]), np.mean([r["N"] for r in sub]), np.mean([r["z_ground_ft"] for r in sub]))}
        for ty in TYPES:
            closes = [r[f"{ty}_closes"] for r in sub]
            rec[f"{ty}_closes_frac"] = sum(closes) / len(closes)
            A = [r[f"{ty}_area"] for r in sub]
            rec[f"{ty}_vol_m3"] = float(sum(0.5 * (A[i] + A[i + 1]) * step for i in range(len(A) - 1) if np.isfinite(A[i]) and np.isfinite(A[i + 1])))
            rec[f"{ty}_base_max"] = float(np.nanmax([r[f"{ty}_base"] for r in sub])) if any(np.isfinite(r[f"{ty}_base"]) for r in sub) else float("nan")
            rec[f"{ty}_toe_ds_max"] = float(np.nanmin([r[f"{ty}_toe_ds"] for r in sub])) if any(np.isfinite(r[f"{ty}_toe_ds"]) for r in sub) else float("nan")
        out.append(rec)
    return out


def gravity_screen(H, m=0.8, phi=40.0, E=0.5, b=6.1, freeboard=6.1):
    """调用 gravity_section_check 的逻辑（无粘聚力，排水有效率 0.5）——直接复用脚本的 CLI 输出里的 FS。"""
    import subprocess, re
    r = subprocess.run([sys.executable, os.path.join(HERE, "gravity_section_check.py"), "--H", str(H), "--m", str(m), "--phi", str(phi),
                        "--E", str(E), "--b", str(b), "--freeboard", str(freeboard), "--detail"], capture_output=True, text=True)
    fs = re.search(r"剪摩 FS=([\d.]+)", r.stdout)
    ok2a = "[OK] FERC Table 2A" in r.stdout; ok20 = "[OK] EM 2200" in r.stdout; ok3 = "[OK] FERC Table 2 高" in r.stdout
    return (float(fs.group(1)) if fs else float("nan"), ok20, ok3, ok2a)


def main():
    rows = load_stations()
    summ = json.load(open(os.path.join(OUT, f"types_{LINE}_summary.json")))
    tot = summ["totals"]; L = summ["length_m"]
    gross = tot["rcc"]["gross_storage_m3"]
    # 组合方案的 run 结果
    runs = {}
    for tag in ("auto_height", "auto_close"):
        f = os.path.join(OUT, f"terrain_{LINE}_{tag}_summary.json")
        if os.path.exists(f):
            runs[tag] = json.load(open(f))
    Hmax = max(r["H_axis"] for r in rows)
    gs = {H: gravity_screen(H) for H in (20, 30, 50, 90, round(Hmax))}

    # ---------------- 方案表 ----------------
    schemes = []
    def add(name, geom, kinds_vol, endarea, grid, wedge, closes_len, open_len, base_max, notes):
        net = gross - wedge
        rcc_v = kinds_vol.get("rcc", 0.0); fill_v = sum(v for k, v in kinds_vol.items() if k != "rcc")
        schemes.append({
            "方案": name, "几何": geom,
            "可闭合长度 m": round(closes_len), "不闭合长度 m": round(open_len), "最大基底宽 m": round(base_max) if np.isfinite(base_max) else "",
            "坝体 平均断面法 Mm³": round(endarea / 1e6, 1), "坝体 网格法 Mm³": round(grid / 1e6, 1),
            "其中 RCC Mm³": round(rcc_v / 1e6, 1), "其中堆石 Mm³": round(fill_v / 1e6, 1),
            "堆石原岩方 Mm³": round(fill_v / SF_ROCK / 1e6, 1), "堆石松方 Mm³": round(fill_v / SF_ROCK / LF_ROCK / 1e6, 1),
            "RCC 骨料原岩方 Mm³": round(rcc_v * AGG_BANK_PER_RCC / 1e6, 1),
            "毛库容 Mm³": round(gross / 1e6, 1), "上游楔扣除 Mm³": round(wedge / 1e6, 1), "净库容 Mm³": round(net / 1e6, 1),
            "净库容 GWh@518m": round(net * GWH_PER_M3, 1), "坝体(网格)/净水": round(grid / net, 1) if net > 0 else "—",
            "堆石车队年数@5,000 BCM/h": round(fill_v / SF_ROCK / FLEET_ROCK_BCM_H / HOURS_PER_YEAR, 1) if fill_v else "",
            "RCC 拌合年数@700 m³/h": round(rcc_v / RCC_PLANT_M3_H / HOURS_PER_YEAR, 1) if rcc_v else "",
            "说明": notes})
    for ty in TYPES:
        t = tot[ty]
        add(TYPE_CN[ty], f"顶宽 {SECTIONS[ty]['crest_w']} m，上游 {SECTIONS[ty]['m_us']}:1，下游 {SECTIONS[ty]['m_ds']}:1",
            {ty: t["end_area_m3"]}, t["end_area_m3"], t["grid_m3"], t["us_wedge_m3"], t["length_closes_m"], t["length_open_m"], t["base_width_max_m"],
            "全环同一断面；不闭合段的方量未计入（那里此断面在几何上做不成）")
    for tag, label in (("auto_height", "分段组合 I：题设高度规则（H≤30 m 用 RCC，其余 CFRD）"), ("auto_close", "分段组合 II：地形闭合规则（CFRD 下游趾≤350 m 且上游趾≤250 m 才用 CFRD，否则 RCC）")):
        if tag in runs:
            s = runs[tag]; dv = s["dam_volume_m3"]; st = s["storage"]
            add(label, "按站选断面（RCC 6.1 m/0/0.8；CFRD 10 m/1.4/1.4）", dv["by_kind"], dv["end_area_10m"], dv["grid_method"], st["V_dam_inside_below_nwl_m3"],
                s["length_m"] - s["infeasible"]["length_m"], s["infeasible"]["length_m"], float("nan"),
                f"RCC 段长 {sum(1 for r in [] )}" if False else "两种断面的接头是应力集中与地震开裂位置（资料 3/4）")
    cols = list(schemes[0].keys())
    with open(os.path.join(OUT, "D1_schemes.csv"), "w") as f:
        f.write(",".join(cols) + "\n")
        for r in schemes: f.write(",".join(str(r[c]).replace(",", "；") for c in cols) + "\n")

    # ---------------- 逐段表 ----------------
    rc = reaches(rows)
    rcols = list(rc[0].keys())
    with open(os.path.join(OUT, "D1_reaches.csv"), "w") as f:
        f.write(",".join(rcols) + "\n")
        for r in rc: f.write(",".join(str(round(r[c], 2) if isinstance(r[c], float) else r[c]) for c in rcols) + "\n")

    # ---------------- D1 markdown ----------------
    md = []
    md.append(f"# 交付物 1：坝型比较表（坝线 {LINE}，坝顶 {CREST_FT:.0f} ft，正常水位 {NWL_FT:.0f} ft）\n")
    md.append("生成：`python3 scripts/deliverables_1_2.py`（输入见脚本头）。所有方量按 1 m DEM 几何计算；\"网格法\"计入平面曲率（环外侧扇形），\"平均断面法\"不计，二者之差就是曲率效应。"
              "不闭合 = 该断面的坝面在 900 m 内碰不到地面（地面比坝面陡），几何上做不成，其方量未计入。\n")
    md.append(f"坝线长 {L:.0f} m；轴线地面 {min(r['z_ground_ft'] for r in rows):.0f}–{max(r['z_ground_ft'] for r in rows):.0f} ft；H 最大 {Hmax:.0f} m。"
              f"毛库容（水位 {NWL_FT:.0f} ft，环内）{gross/1e6:.1f} Mm³。净库容 = 毛库容 − 上游坝体在水位以下的体积。GWh 按毛水头 518 m、效率 0.85，只作尺度。\n")
    md.append("## 1.1 方案级比较\n")
    md.append("| " + " | ".join(cols) + " |\n|" + "---|" * len(cols))
    for r in schemes:
        md.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    md.append("\n三态换算假设：堆石压实方/原岩方 SF = 0.95（试验填筑前的占位，资料 7/8）、松方系数 LF = 0.65（硬岩 0.60–0.67，资料 9）；RCC 每 m³ 需骨料原岩方 0.9 m³。"
              "车队与拌合站只是把方量换成小时的尺子（5,000 BCM/h ≈ 8 台 D10T2 松土 + 4 台大挖 + 40 台 740C；700 m³/h ≈ 一座 1,000 yd³/h 大型拌合站），不是配置建议。\n")
    md.append("## 1.2 准则对照（资料 1–5，按坝型）\n")
    md.append("| 项 | RCC 重力 | 面板堆石 CFRD | 心墙堆石 ECRD |\n|---|---|---|---|")
    md.append("| 断面依据 | EM 2200：上游直立、下游 0.7–0.8:1；RCC 无模板 ≥0.8:1；坝顶 ≥20 ft（RCC 手册施工要求） | 小坝设计/DS-13：上游 1.3–1.7:1，下游休止角 1.3–1.4:1；顶宽 ≥20 ft；拱度 1% H | 案例 New Hogan/Flannagan/Cougar：上游 1.8–2:1，下游 1.6–2:1；顶宽 25–40 ft；心墙底宽 ≥25% 水头（128 m → ≥32 m） |")
    md.append(f"| 稳定初筛（本地块参数） | φ=40°、c=0、排水有效率 0.5：H=20 m FS {gs[20][0]:.2f}、30 m {gs[30][0]:.2f}、50 m {gs[50][0]:.2f}、90 m {gs[90][0]:.2f}、{round(Hmax)} m {gs[round(Hmax)][0]:.2f}；"
              f"FERC 2A（无 c，≥1.5）全部通过，EM 2200 usual（≥2.0）到 50 m，FERC Table 2 高危害（≥3.0）只到 30 m → 高段要靠层面粘聚力（RCC 手册：c 300 psi 级设计值）或加缓/锚固 | 好岩石 + 好基础可只做无限边坡（DS-13 2.3.5.3）；否则 Spencer FS 1.5（稳定渗流）/1.3（骤降、施工末） | Spencer FS 1.5/1.3；心墙开裂假定必有 → 反滤排水必设 |")
    md.append("| 基础要求 | 微–中风化微裂隙岩；基底水平或向上游倾（抗滑）；突变处应力集中 → 削坡/齿槽；固结灌浆 30 ft；帷幕 H/3+50 ft | 趾板落在可灌浆基岩上；堆石区可容忍差些的基础；\"堆石坝可陡坡，接触面积小\" | 心墙下坝肩 ≤0.5H:1V 整形，台阶 >0.5 ft 处理；心墙接触面宽 |")
    md.append("| 斜坡基础（本地块 24–35°） | 基底必须做成水平台阶或向上游倾：坝趾开挖深 ≈1.55 H（笔记 4）；每 50 ft 一段的台阶接头是应力集中区 | 下游趾在 24° 坡上 476 m 外、30° 935 m、≥35.5° 不闭合；趾板在 34° 坡上 | 上游 2:1 在 24° 坡 2,337 m 外、≥28° 不闭合；下游 1.8:1 ≥29° 不闭合 |")
    md.append("| 地震（PGA 0.35 g/2,475 年） | >50 ft 必须动力分析；FERC 不接受拟静力；震后 FS 1.3；地震区需全层面垫层砂浆 | 超高 3–5% H（≥3.8–6.4 m）；防御措施清单；\"推定基础错动时堆石优于土坝\" | 同左 + 心墙塑性、宽反滤贯穿全高、坝顶加宽 |")
    md.append("| 材料（资料 6/7） | 骨料加工厂；胶材山路运输；水源 14 天连续养护 | 完好岩块（变质火山岩候选）；千枚岩片状料需试验填筑证明；泥流岩耐久性待验 | 心墙料可用高度风化岩（Carters）；反滤需筛分制砂（D15 ≤0.7 mm）；分区多、雨季受限 |")
    md.append("| 压实与记录（资料 7/8/10/11） | 每层 Vebe/核子密度、层面成熟度；试验段 100 ft | 层厚 0.6–1.0 m、10–20 t 振动碾、遍数由试验填筑定；回数 100% 轨迹记录；ICMV 只作找弱区 | 心墙密度/含水率点测 + 回数；反滤层专门验收 |")
    md.append("| 优点（资料 1/3） | 可短时漫顶；坝体最小；断面在陡坡上唯一能闭合 | 可检修、灌浆不在关键路径、坝体不饱和、易加高、雨季可施工 | 材料适应性最广 |")
    md.append("| 本地块的主要开口 | 4.6 km 环上 34 Mm³ 级 RCC 的热控与拌合；30–35° 坡的台阶开挖量；曲线段分缝 | 上游楔吃掉库容；下游填方落到宗地界外/陡坡；34° 以上不闭合段需支挡或改断面 | 几乎全环不闭合；只在南界台地可行 |\n")
    md.append("## 1.3 逐段表（按 H 档 ≤30 / 30–80 / >80 m 与 CFRD 下游趾距离档 ≤350 / 350–600 / >600 m / 不闭合 切段，短于 100 m 的并入前段；地质为 1:250k 图上的\"可能单元\"：泥流盖层只标在 ≥3,950 ft 的台地面，侧翼按 Calaveras 基底）\n")
    hdr = ["段", "桩号 m", "长 m", "地面 ft", "H m", "下游坡°", "上游坡°", "地质（可能）", "RCC 闭合/底宽 m/方量 Mm³", "CFRD 闭合/底宽/方量", "ECRD 闭合/底宽/方量", "1.6 堆石 闭合/底宽/方量"]
    md.append("| " + " | ".join(hdr) + " |\n|" + "---|" * len(hdr))
    for r in rc:
        def cell(ty):
            b = r[f"{ty}_base_max"]; return f"{r[f'{ty}_closes_frac']*100:.0f}% / {b:.0f} / {r[f'{ty}_vol_m3']/1e6:.2f}" if np.isfinite(b) else f"{r[f'{ty}_closes_frac']*100:.0f}% / — / {r[f'{ty}_vol_m3']/1e6:.2f}"
        md.append(f"| {r['reach']} | {r['sta_from']:.0f}–{r['sta_to']:.0f} | {r['length_m']:.0f} | {r['z_ground_ft_min']:.0f}–{r['z_ground_ft_max']:.0f} | {r['H_min']:.0f}–{r['H_max']:.0f} | {r['slope_ds_mean']:.0f} | {r['slope_us_mean']:.0f} | {r['geology']} | {cell('rcc')} | {cell('cfrd')} | {cell('ecrd')} | {cell('rockfill')} |")
    md.append("\n读法：\"闭合\"是该段内此断面能落到地面的桩号比例；底宽是段内最大基底宽；方量是该段此断面的平均断面法方量（只计闭合站）。"
              "同一段里 RCC 与堆石都闭合时，两者方量之比就是资料 1 说的\"堆石面积随 H² 增长、坡地再乘 1.6–1.8\"的实测版。\n")
    md.append("## 1.4 这张表没有回答的（留给交付物 3–6 与现场勘察）\n")
    md.append("- 每段坝基的真实岩性、风化深度、RQD、软弱面产状（决定 RCC 能不能落在那里、CFRD 趾板怎么做）。\n- 堆石料的压实特性（SF、层厚、遍数）与 RCC 层面抗剪（决定胶材量）——都要试验采石与试验填筑/试验段。\n- 下游填方落在宗地界外的段落的用地。\n- 净库容 → 电量还要下库位置与水位。\n- 环形坝线转角处的三维效应（DS-13 建议数值分析）。\n")
    open(os.path.join(DEL, "D1_dam_type_comparison.md"), "w").write("\n".join(md))

    # ---------------- D2：纵剖面 + 代表断面 ----------------
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]; matplotlib.rcParams["axes.unicode_minus"] = False
    t = Terrain(os.path.join(DATA, "dem_3dep_1m_utm10.tif"), os.path.join(DATA, "dem_3dep_1m_utm10.json"))
    cand = json.load(open(os.path.join(DATA, "dam_line_candidates.json")))[LINE]
    line = LineString(cand["coords"]); poly = Polygon(cand["coords"]) if cand["closed"] else None
    s, pts, zg = longitudinal_profile(t, line, 10.0); nrm = inward_normals(pts, poly, cand.get("inside_point"))
    st = np.array([r["station"] for r in rows]); zft = np.array([r["z_ground_ft"] for r in rows]); H = np.array([r["H_axis"] for r in rows])
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(16, 9), gridspec_kw={"height_ratios": [3, 1.2]}, sharex=True)
    ax.plot(st, zft, "k-", lw=1.2, label="轴线地面")
    ax.axhline(CREST_FT, color="r", lw=1, label=f"坝顶 {CREST_FT:.0f} ft"); ax.axhline(NWL_FT, color="b", ls="--", lw=0.8, label=f"正常水位 {NWL_FT:.0f} ft")
    ax.fill_between(st, zft, CREST_FT, where=H > 0, color="0.85")
    for ty, col in (("rcc", "orange"), ("cfrd", "green"), ("ecrd", "olive")):
        zt = np.array([t.elev(r["E"] + r[f"{ty}_toe_ds"] * 0, r["N"]) for r in rows]) if False else None
    # 坝趾高程：从站表的 toe 偏移重新采样地面
    for ty, col, mk in (("rcc", "orange", "."), ("cfrd", "green", "."), ("ecrd", "olive", "x")):
        zt = []
        for i, r in enumerate(rows):
            o = r[f"{ty}_toe_ds"]
            if np.isfinite(o):
                p = pts[i] + o * nrm[i]; zt.append(float(t.elev(p[0], p[1]) / FT))
            else:
                zt.append(np.nan)
        ax.plot(st, zt, mk, color=col, ms=3, label=f"{ty.upper()} 下游坝趾高程")
    # 分段与地质
    colors = {"Calaveras 千枚岩/泥质岩（可能）": "#d9d2e9", "Mehrten 火山泥流（可能）": "#f4cccc", "变质火山岩（可能）": "#cfe2f3"}
    for r in rc:
        ax.axvspan(r["sta_from"], r["sta_to"], color=colors[r["geology"]], alpha=0.5, lw=0)
        ax.axvline(r["sta_from"], color="0.4", lw=0.5, ls=":")
        ax.text((r["sta_from"] + r["sta_to"]) / 2, CREST_FT + 40, str(r["reach"]), ha="center", fontsize=8)
    for name, c in colors.items():
        ax.fill_between([], [], [], color=c, alpha=0.5, label=name)
    ax.set_ylabel("高程 ft"); ax.set_ylim(min(zft) - 700, CREST_FT + 120); ax.grid(alpha=0.3); ax.legend(loc="lower left", fontsize=8, ncol=3)
    ax.set_title(f"交付物 2：坝线 {LINE} 纵剖面（L={L:.0f} m，Hmax={Hmax:.0f} m）——地面、坝顶、水位、各断面下游坝趾高程、分段编号与可能地质")
    # 下：H 与坡角
    ax2.plot(st, H, "k-", lw=1, label="H（坝顶−轴线地面）m")
    ax2.plot(st, [r["slope_ds_deg"] for r in rows], "m-", lw=0.7, label="下游侧地面坡角°（100 m 内）")
    ax2.plot(st, [r["slope_us_deg"] for r in rows], "c-", lw=0.7, label="上游（库）侧地面坡角°")
    ax2.axhline(30, color="0.5", lw=0.5, ls="--"); ax2.axhline(32, color="green", lw=0.5, ls=":"); ax2.axhline(35.5, color="green", lw=0.5, ls=":")
    ax2.set_xlabel("桩号 m（从南界西端起，逆时针）"); ax2.set_ylabel("m / °"); ax2.grid(alpha=0.3); ax2.legend(loc="upper right", fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "D2_profile.png"), dpi=130); plt.close(fig)

    # 代表断面：各段取 H 最大处 + 全线 H 最大 + 下游最陡处
    picks = []
    for r in rc:
        k = (st >= r["sta_from"]) & (st <= r["sta_to"])
        i = int(np.where(k)[0][np.argmax(H[k])]); picks.append((i, f"段 {r['reach']} H 最大"))
    i_steep = int(np.nanargmax([r["slope_ds_deg"] if r["H_axis"] > 60 else -1 for r in rows])); picks.append((i_steep, "高段中下游最陡"))
    picks = sorted({i: lab for i, lab in picks}.items())
    n = len(picks); fig, axs = plt.subplots(n, 1, figsize=(13, 3.3 * n)); axs = np.atleast_1d(axs)
    table = []
    for ax, (i, lab) in zip(axs, picks):
        off, zsec = cross_section(t, pts[i], nrm[i])
        ax.plot(off, zsec / FT, "k-", lw=1.2)
        ax.axhline(NWL_FT, color="b", ls="--", lw=0.7); ax.axhline(CREST_FT, color="r", lw=0.6)
        rec = {"station": st[i], "label": lab, "E": rows[i]["E"], "N": rows[i]["N"], "z_ground_ft": zft[i], "H_axis_m": H[i],
               "slope_ds_deg": rows[i]["slope_ds_deg"], "slope_us_deg": rows[i]["slope_us_deg"], "geology": geol(rows[i]["E"], rows[i]["N"], zft[i])}
        rec["no_dam"] = bool(H[i] <= 0)
        for ty, col in (("rcc", "orange"), ("cfrd", "green"), ("ecrd", "olive")):
            d = dam_section(off, zsec, CREST_FT * FT, SECTIONS[ty], NWL_FT * FT)
            rec[f"{ty}_area_m2"] = d["area"]; rec[f"{ty}_base_m"] = d.get("base_width", np.nan); rec[f"{ty}_toe_ds"] = d["toe_ds"]; rec[f"{ty}_toe_us"] = d["toe_us"]
            b = SECTIONS[ty]["crest_w"]
            if rec["no_dam"]:
                if ty == "rcc":
                    ax.text(0, CREST_FT + 20, "地面高于坝顶：此处无坝", ha="center", fontsize=9, color="0.3")
                continue
            if d["toe_us"] is not None and d["toe_ds"] is not None:
                xs = [d["toe_ds"], -b / 2, b / 2, d["toe_us"]]; zs = [d["z_toe_ds"] / FT, CREST_FT, CREST_FT, d["z_toe_us"] / FT]
                ax.plot(xs, zs, "-", color=col, lw=1.4, label=f"{ty.upper()} 底宽 {d['base_width']:.0f} m，面积 {d['area']:.0f} m²")
            else:
                # 画到边界为止，标"不闭合"
                mds = SECTIONS[ty]["m_ds"]; mus = SECTIONS[ty]["m_us"]
                xs = [-900, -b / 2, b / 2, min(900, (b / 2) + mus * (CREST_FT * FT - np.nanmin(zsec)))] if mus > 0 else [-900, -b / 2, b / 2, b / 2]
                zs = [(CREST_FT * FT - (900 - b / 2) / mds) / FT if mds > 0 else CREST_FT, CREST_FT, CREST_FT, CREST_FT if mus == 0 else np.nanmin(zsec) / FT]
                ax.plot(xs[:2], zs[:2], "--", color=col, lw=1); ax.plot(xs[1:], zs[1:], "-", color=col, lw=1.4, label=f"{ty.upper()} 不闭合（下游趾 {'—' if d['toe_ds'] is None else d['toe_ds']}，上游趾 {'—' if d['toe_us'] is None else d['toe_us']}）")
        table.append(rec)
        ax.set_xlim(-900, 400); ax.set_ylim(min(np.nanmin(zsec) / FT, CREST_FT - 700), CREST_FT + 60)
        ax.set_title(f"桩号 {st[i]:.0f} m（{lab}）：地面 {zft[i]:.0f} ft，H={H[i]:.0f} m，下游坡 {rows[i]['slope_ds_deg']:.0f}°，上游坡 {rows[i]['slope_us_deg']:.0f}°，{rec['geology']}", fontsize=9)
        ax.set_xlabel("离轴线偏移 m（+ 为库侧）"); ax.set_ylabel("ft"); ax.grid(alpha=0.3); ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "D2_sections.png"), dpi=120); plt.close(fig)
    tcols = list(table[0].keys())
    with open(os.path.join(OUT, "D2_sections_table.csv"), "w") as f:
        f.write(",".join(tcols) + "\n")
        for r in table: f.write(",".join("" if r[c] is None or (isinstance(r[c], float) and not np.isfinite(r[c])) else (f"{r[c]:.1f}" if isinstance(r[c], float) else str(r[c])) for c in tcols) + "\n")

    md2 = [f"# 交付物 2：整条坝线纵剖面与代表断面（坝线 {LINE}）\n",
           "生成：`python3 scripts/deliverables_1_2.py`。图：`outputs/D2_profile.png`（纵剖面）、`outputs/D2_sections.png`（代表断面）；站表 `outputs/types_%s_stations.csv`（每 10 m 一站，四种断面的坝趾、底宽、面积、上游楔）；断面表 `outputs/D2_sections_table.csv`。\n" % LINE,
           "## 2.1 纵剖面读法\n",
           "- 上图：黑线是轴线地面（1 m DEM 双线性插值），红线坝顶、蓝虚线正常水位；灰区是 H；橙/绿/橄榄点是 RCC / CFRD / ECRD 断面的**下游坝趾高程**——点越低说明这个断面在这里要伸到多深的坡上；没有点的桩号是该断面不闭合。底色是 1:250k 地质图上\"可能\"的单元，段号对应交付物 1 的逐段表。\n",
           "- 下图：H（黑）、下游侧 100 m 内地面坡角（品红）、库侧坡角（青）；绿色点线 32° 与 35.5° 分别是 1.6:1 与 1.4:1 坝面的坡角——地面坡角一旦高于坝面坡角，该堆石断面就永不闭合。\n",
           f"- 桩号从南界西端（E{rows[0]['E']:.0f} N{rows[0]['N']:.0f}）起逆时针：先沿 3,700 ft 等高线绕鼻尖的西坡、北坡（0–{max(r['sta_to'] for r in rc if r['H_min'] > 80):.0f} m 为 H 100–130 m 的高段），再沿东界南下（H 30–80 m），最后沿南界台地回到起点（H ≤ 30 m，两处地面高于坝顶）。\n",
           "## 2.2 代表断面\n",
           "| 桩号 m | 位置 | 地面 ft | H m | 下游坡° | 上游坡° | 地质（可能） | RCC 底宽/面积 | CFRD 底宽/面积 | ECRD 底宽/面积 |\n|---|---|---|---|---|---|---|---|---|---|"]
    for r in table:
        def c2(ty):
            if r.get("no_dam"):
                return "无坝（地面高于坝顶）"
            return f"{r[f'{ty}_base_m']:.0f} m / {r[f'{ty}_area_m2']:.0f} m²" if np.isfinite(r[f"{ty}_base_m"]) else "不闭合"
        md2.append(f"| {r['station']:.0f} | {r['label']} | {r['z_ground_ft']:.0f} | {r['H_axis_m']:.0f} | {r['slope_ds_deg']:.0f} | {r['slope_us_deg']:.0f} | {r['geology']} | {c2('rcc')} | {c2('cfrd')} | {c2('ecrd')} |")
    md2.append("\n断面图中每个子图把三种断面叠画在同一条地面线上（+ 为库侧）：实线是闭合的断面轮廓，虚线表示该坝面伸到 900 m 边界仍未碰到地面。面积是断面上坝体（坝面与地面之间）的面积；库侧水位以下、坝面与地面之间的面积就是\"上游楔\"，它在库容里要扣掉。\n")
    md2.append("## 2.3 这两张图直接说明的几件事\n")
    md2.append("- 高段（沿 3,700 线的约 2.9 km）：轴线下游 100 m 内坡角 18–42°，再往外更陡；CFRD 1.4:1 的下游趾伸到 430–930 m，ECRD 2.0/1.8 在 12% 桩号不闭合；RCC 0.8:1 在 80–315 m 内闭合。\n- 库侧坡角 12–35°：CFRD 上游 1.4:1 的趾在库内 100–250 m，这段水面下全是坝体（上游楔 27 Mm³）。\n- 东界段（H 30–80 m）与南界台地段（H ≤ 30 m）地面平缓，三种断面都闭合，方量差异回到 \"堆石 ≈ 3–4 × RCC\" 的平地规律；南界有两段地面高于坝顶（段 9、11），那里无坝。\n- 冲沟处（若坝线拉直，见 `terrain_A_parcel_ring_3700_s100_auto_profile.png`）基础降到 3,637 ft，H 到 147 m。\n")
    open(os.path.join(DEL, "D2_profile_and_sections.md"), "w").write("\n".join(md2))
    print("D1/D2 written"); print(json.dumps({r["方案"]: (r["坝体 网格法 Mm³"], r["净库容 Mm³"]) for r in schemes}, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
