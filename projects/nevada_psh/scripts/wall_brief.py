#!/usr/bin/env python3
""""墙段自由组合"结果的简版交付（HTML）：读 outputs/wall_curves_rcc.csv、wall_plan_eval.csv、wall_plan_segments.csv 与图，回答
"能否用 B1 三分之一到一半的填方（9–13 Mm³）换到 B1 一半以上的库容（≥ 30 Mm³）"。"""
import base64, csv, json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); OUT = os.path.join(ROOT, "outputs"); DEL = os.path.join(ROOT, "deliverables")
f = lambda v: float(v) if v not in ("", None) else float("nan")
FT = 0.3048; ACRE = 4046.8564224
B1 = {"len_km": 3.9, "fill": 29.7, "need": 26.7, "nat": 26.5, "cut": 36.3, "net": 61.4, "nwl": 4030.0, "gwh": 70.6}
GUL = {"len_km": 0.84, "fill": 2.1, "nat": 6.4, "cut": 0.6, "net": 7.0, "nwl": 3980.0}
gwh = lambda v_Mm3, nwl_ft: v_Mm3 * 1e6 * 1000 * 9.81 * (nwl_ft - 2400) * FT * 0.85 / 3.6e12
curves = list(csv.DictReader(open(os.path.join(OUT, "wall_curves_rcc.csv"))))
ev = list(csv.DictReader(open(os.path.join(OUT, "wall_plan_eval.csv"))))
segs = list(csv.DictReader(open(os.path.join(OUT, "wall_plan_segments.csv"))))
cf_ev = {r["crest_ft"]: r for r in ev}  # 注意：4050 有两条（λ 2.1 与 2.4）
img64 = lambda p: "data:image/png;base64," + base64.b64encode(open(p, "rb").read()).decode()
esc = lambda x: str(x).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
fig = lambda p, cap: f'<figure><img src="{img64(p)}"><figcaption>{cap}</figcaption></figure>'
def table(cols, rows):
    h = ['<div class="wrap"><table><tr>' + "".join(f"<th>{esc(c)}</th>" for c in cols) + "</tr>"]
    for r in rows: h.append("<tr>" + "".join(f"<td>{esc(x)}</td>" for x in r) + "</tr>")
    return "\n".join(h) + "</table></div>"
CSS = """body{font-family:"Noto Sans CJK SC","WenQuanYi Zen Hei","PingFang SC","Microsoft YaHei",Arial,sans-serif;max-width:1150px;margin:0 auto;padding:18px 26px 50px;line-height:1.6;color:#1f2328;background:#fff}
h1{font-size:21px;border-bottom:2px solid #333;padding-bottom:6px} h2{font-size:17px;margin-top:28px;border-left:5px solid #7b1fa2;padding-left:10px}
table{border-collapse:collapse;font-size:13px;margin:8px 0 12px} th,td{border:1px solid #bfc4cc;padding:4px 8px;vertical-align:top;text-align:left} th{background:#eef1f5}
tr:nth-child(even) td{background:#f8f9fb} .wrap{overflow-x:auto} figure{margin:12px 0} img{max-width:100%;border:1px solid #d8dbe0} figcaption{font-size:12.5px;color:#444}
.box{background:#f3e5f5;border-left:4px solid #7b1fa2;padding:8px 14px;margin:10px 0} .ans{background:#e8f5e9;border-left:4px solid #2e7d32;padding:10px 14px;margin:10px 0}
.warn{background:#fff8e1;border-left:4px solid #f0b400;padding:8px 14px;margin:10px 0} .f{font-family:Menlo,Consolas,monospace;background:#f1f3f5;padding:1px 5px;border-radius:3px;font-size:12.5px} .small{font-size:12px;color:#555} li{margin:4px 0}"""
doc = lambda title, body: f'<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{esc(title)}</title><style>{CSS}</style></head><body>{body}</body></html>'

def curve_pts(r):
    cur = [tuple(map(float, c.split(":"))) for c in r["curve_rcc"].split(";")] if r.get("curve_rcc") else []
    wedge = f(r["gross_Mm3"]) - f(r["net_rcc_Mm3"]); cuts = np.array([c[1] for c in cur]); gs = np.array([c[2] for c in cur])
    at = lambda k: (float(np.interp(k, cuts, gs)) - wedge) if len(cuts) and cuts.max() >= k else float("nan")
    return cur, wedge, at
REC = next(r for r in ev if r["crest_ft"].startswith("4050") and abs(f(r["lambda"]) - 2.1) < 1e-6)
cur, wedge, at = curve_pts(REC)
rec_segs = [s for s in segs if s["crest_ft"].startswith("4050") and abs(f(s["lambda"]) - 2.1) < 1e-6]

h = []
h.append("<h1>墙段自由组合：用最少的填方围出最多的库容（简版）</h1>")
h.append('<div class="box"><b>问题。</b>不预设闭合环，把坝当作可自由组合的墙段：南面与东面的天然高地（3,950–4,230 ft）免费，只在必要处修墙。从零开始每次加进"增量库容 ÷ 增量填方"最高的墙段，画累计填方–累计库容曲线，找拐点；再看围住的区域里挖深能加多少库容、出多少料。'
         '<b>核心问题：</b>有没有一个方案，只用 B1 三分之一到一半的填方（9–13 Mm³），换到 B1 一半以上的库容（≥ 30 Mm³）？</div>')
h.append(f'<div class="ans"><b>答案：有。</b>坝顶 4,050 / 水位 4,030 ft（与 B1 相同水头）下，只修 3 段墙共 {f(REC["dam_len_m"]):,.0f} m（其余 {f(REC["free_len_m"]):,.0f} m 库岸是高于坝顶的天然地面），RCC 填方 {f(REC["fill_rcc_Mm3"]):.1f} Mm³（B1 的 {100*f(REC["fill_rcc_Mm3"])/B1["fill"]:.0f}%），'
         f'天然库容 {f(REC["net_rcc_Mm3"]):.1f} Mm³（B1 不挖 26.5 的 {100*f(REC["net_rcc_Mm3"])/B1["nat"]:.0f}%）；围住的区域再挖 20 Mm³（B1 挖 36）得 {at(20):.0f} Mm³（B1 61.4 的 {100*at(20)/B1["net"]:.0f}%），挖到底 {f(REC["deep_floor_ft_rcc"]):.0f} ft（挖 {f(REC["deep_cut_Mm3_rcc"]):.0f}）得 {f(REC["deep_net_Mm3_rcc"]):.0f} Mm³。'
         f'挖方 {f(REC["deep_cut_Mm3_rcc"]):.0f} Mm³ 对 RCC 骨料需求 {f(REC["need_rcc_Mm3"]):.1f} Mm³ 是 {f(REC["deep_cut_over_need_rcc"]):.1f} 倍，料绰绰有余。这个方案就是 B1 的"东半环"：围住山鼻东侧那条主沟的流域和东半块台地，放弃西半块（B1 里 1.3 km 沿 3,700 线的高坝换来的库容最少）。</div>')
h.append('<h2>1. 我做了什么（等价于"每次加最高效的墙段"）</h2>')
h.append('<ul><li>把地块上 4 m 一格的每个格子看成"划入库区"或"不划入"。划入的收益是它在水位以下的蓄水量（水位 − 地面 × 面积）；库区边界穿过低于坝顶的地面时要付一段坝：填方 = 断面积(H, 局部坡度) × 长度，穿过高于坝顶的天然地面免费。</li>'
         '<li>给填方定一个价格 λ（每 m³ 填方要换到 λ m³ 库容才值得），求"总库容 − λ × 总填方"最大的库区轮廓。这是一个最小割问题，用最大流精确求解（和露天矿的最优坑界是同一个数学问题）。</li>'
         '<li>λ 从 40 降到 0.4：λ 大时只留效率最高的墙段，λ 小时墙越加越多，解是嵌套的——这就是"每次加进增量库容÷增量填方最高的墙段"的精确版。相邻两个解之间的 Δ库容/Δ填方 就是边际效率。</li>'
         '<li>坝顶试 3,900 / 3,950 / 4,000 / 4,050（水位 = 坝顶 − 20 ft）；墙高上限 130 m（与 B1 的 129 m 一致）；搜索用 RCC 断面成本，入围方案在 1 m 地形上逐站放坝精算 RCC 与 CFRD 两种填方，并按笔记 14 的碗形坑规则算"围住后再挖"。</li></ul>')
h.append("<h2>2. 累计填方–累计库容曲线</h2>")
h.append(fig(os.path.join(OUT, "wall_curves_rcc.png"), "图 1 每条曲线一个坝顶，从原点往右 = 填方价格 λ 降低、墙段增加；点旁小数字是该步的边际效率（Δ库容 ÷ Δ填方）。曲线本身是 4 m 网格的估算（填方偏高约 30%），方块是 1 m 精算的天然库容，三角是再挖后的库容。红星 B1，黑菱形纯沟谷坝。绿色是目标区。"))
rows = []
c4050 = [r for r in curves if r["crest_ft"].startswith("4050") and f(r["fill_est_Mm3"]) > 0.05]
for r in c4050: rows.append([f(r["lambda"]), f"{f(r['fill_est_Mm3']):.1f}", f"{f(r['storage_est_Mm3']):.1f}", f"{f(r['dam_len_m']):,.0f}", f"{f(r['free_len_m']):,.0f}", f"{f(r['area_acre']):.0f}", r["marg_eff"] or "—"])
h.append("<p>坝顶 4,050 的曲线逐步读（4 m 网格估算值）：</p>")
h.append(table(["λ（填方价格）", "累计填方 Mm³", "累计天然库容 Mm³", "墙总长 m", "天然岸长 m", "库区面积 acre", "这一步的边际效率"], rows))
h.append('<div class="warn"><b>拐点在哪。</b>λ 从 2.4 降到 1.85 时加进的是东半环（主沟流域 + 东半块台地），边际效率 1.0–1.6；再往下（λ 1.6–1.2）只是把墙往外挪几十米，效率 0.75–1.2；到 λ 0.85–0.7 加进西半环（沿 3,700 线绕鼻尖的 1.3 km 高坝，就是 B1 的第 1–3 段），边际效率只有 0.5–0.8——每 m³ 填方换不到 1 m³ 库容，还不如去挖（1 m³ 挖方 ≈ 1 m³ 库容，而且出料）。'
         '所以停止点是 λ ≈ 2.1：东半环修完就停。四个坝顶的曲线形状一样，只是水头不同。</div>')
h.append(fig(os.path.join(OUT, "wall_growth_rcc.png"), "图 2 坝顶 4,050 下库区随 λ 降低的生长：先是东半环（左三幅），λ ≤ 0.7 才把西半块台地围进来（右两幅）。"))
h.append("<h2>3. 各坝顶在 9–13 Mm³ 填方附近的方案</h2>")
rows = []
for r in ev:
    c, w, a = curve_pts(r); nwl = f(r["nwl_ft"])
    rows.append([f"{f(r['crest_ft']):,.0f} / {nwl:,.0f}", f"λ={f(r['lambda']):g}", f"{int(f(r['n_segments']))} 段 / {f(r['dam_len_m']):,.0f} m（天然岸 {f(r['free_len_m']):,.0f} m）", f"{f(r['H_max_m']):.0f}",
                 f"{f(r['fill_rcc_Mm3']):.1f} / {f(r['fill_cfrd_Mm3']):.1f}", f"{f(r['net_rcc_Mm3']):.1f}（{f(r['area_acre']):.0f} acre）", f"{f(r['bal_floor_ft_rcc']):,.0f} / {f(r['bal_cut_Mm3_rcc']):.0f} / {f(r['bal_net_Mm3_rcc']):.1f}",
                 f"{a(20):.1f}" if np.isfinite(a(20)) else "—", f"{a(36):.1f}" if np.isfinite(a(36)) else f"饱和 {f(r['deep_net_Mm3_rcc']):.1f}（挖 {f(r['deep_cut_Mm3_rcc']):.0f}）",
                 f"{100*(nwl-2400)/(B1['nwl']-2400):.0f}%", f"{gwh(a(20), nwl):.0f}" if np.isfinite(a(20)) else "—", f"{f(r['rcc_toe_out_m']):.0f}"])
h.append(table(["坝顶 / 水位 ft", "解", "墙段 / 总长", "最大坝高 m", "填方 RCC / CFRD Mm³", "天然净库容 Mm³", "料平衡点：底 ft / 挖 / 库容 Mm³", "挖 20 Mm³ 后库容", "挖 36 Mm³ 后库容", "水头占 B1", "GWh（挖 20）", "下游坝趾出界 m"], rows))
h.append('<p class="small">"料平衡点"= 挖方刚够 RCC 骨料（RCC × 0.9）的库底；"挖 20 / 36 Mm³ 后"按各方案的碗形坑曲线插值（36 是 B1 的挖方），"饱和"= 坑壁相遇、再挖库容不增。GWh 按下库 2,400 ft、η 0.85。坝顶越低，免费的天然岸越长、墙越短，但水头越小；坝顶 4,050 与 B1 同水头。</p>')
h.append(fig(os.path.join(OUT, "wall_plan_selected.png"), "图 3 五个入围方案的平面：橙 = 墙段（长 / 最大坝高 / RCC 填方），绿 = 库区边界中不需要坝的天然岸，蓝 = 天然蓄水面，黄–红 = 再挖深度（深挖参考点）。所有方案都是同一个形状：围主沟流域和东半块台地。"))
h.append("<h2>4. 推荐方案：坝顶 4,050，东半环</h2>")
h.append(fig(os.path.join(OUT, "wall_plan_recommended.png"), "图 4 推荐方案的墙段：W1 沿主沟西侧山脊、横跨沟口（最高 131 m）；W2 沿东北面；W3 东界上一小段低坝。南面与东面其余库岸是 3,950–4,230 ft 的天然高地。"))
rows = []
for s in rec_segs:
    if f(s["len_m"]) < 40: continue
    rows.append([f"W{s['seg']}", f"E {f(s['E']):,.0f} / N {f(s['N']):,.0f}", f"{f(s['len_m']):,.0f}", f"{f(s['H_max_m']):.0f}", f"{f(s['fill_rcc_Mm3']):.2f}", f"{f(s['fill_cfrd_Mm3']):.2f}",
                 {"1": "跨主沟口 + 沿沟西侧山脊向北，是整个方案的主体；沟底 3,620 ft，坝高 131 m", "2": "东北面，把东半块台地的北缘封上；下游坝趾出宗地东界，要邻地或内移", "3": "东界过渡段的低坝（B1 的第 5–6 段），坝趾出界"}.get(s["seg"], "")])
h.append(table(["墙段", "位置（UTM 10N，段中点）", "长 m", "最大坝高 m", "RCC 填方 Mm³", "CFRD 填方 Mm³", "说明"], rows))
h.append(f'<p><b>手算核对。</b>W1：坝高从两端 0 到沟底 131 m，按三角形分布平均 H² ≈ 131²/3 = 5,700 m²，RCC 断面 <span class="f">A ≈ 0.4·H²·1.3（坡地放大）+ 6.1·H</span> ≈ 0.4×5,700×1.3 + 6.1×65 ≈ 3,360 m²，乘 1,370 m ≈ 4.6 Mm³；模型逐站算是 8.6 Mm³，差别来自沟口那 300 m 坝高集中在 100 m 以上（H² 权重）。'
         f'天然库容：水面 {f(REC["area_acre"]):.0f} acre × 平均水深 {f(REC["net_rcc_Mm3"])*1e6/(f(REC["area_acre"])*ACRE):.0f} m = {f(REC["net_rcc_Mm3"]):.1f} Mm³。再挖：碗形坑曲线（库底 : 挖方 / 库容）'
         + "，".join(f"{c[0]:.0f}: {c[1]:.0f}/{c[2]-wedge:.0f}" for c in cur[3:12:2]) + f'——挖到 3,530 ft 之前每 m³ 挖方换约 1 m³ 库容，再深就饱和（坑壁相遇）。</p>')
h.append("<h2>5. 与 B1、纯沟谷坝对照</h2>")
rows = [["B1 全环 RCC 4,050", "3.9", "29.7", "26.5", "36.3", "61.4", "4,030", "70.6"],
        [f"东半环 4,050（推荐，λ 2.1）", f"{f(REC['dam_len_m'])/1000:.1f}", f"{f(REC['fill_rcc_Mm3']):.1f}（{100*f(REC['fill_rcc_Mm3'])/B1['fill']:.0f}%）", f"{f(REC['net_rcc_Mm3']):.1f}（{100*f(REC['net_rcc_Mm3'])/B1['nat']:.0f}%）", "20", f"{at(20):.1f}（{100*at(20)/B1['net']:.0f}%）", "4,030（100%）", f"{gwh(at(20), 4030):.0f}"],
        ["东半环 4,050（λ 2.4，更省）", "2.3", "8.0（27%）", "16.8（63%）", "20", "36.0（59%）", "4,030", f"{gwh(36.0, 4030):.0f}"],
        ["东半环 4,000（λ 2.1）", "2.1", "9.2（31%）", "18.3（69%）", "30", "47.0（77%）", "3,980（97%）", f"{gwh(47.0, 3980):.0f}"],
        ["东半环 3,950（λ 1.6）", "1.9", "8.1（27%）", "16.3（62%）", "36", "49.0（80%）", "3,930（94%）", f"{gwh(49.0, 3930):.0f}"],
        ["纯沟谷坝 S29 4,000", "0.84", "2.1（7%）", "6.4（24%）", "0.6", "7.0（11%）", "3,980（97%）", f"{gwh(7.0, 3980):.0f}"]]
h.append(table(["方案", "坝长 km", "RCC 填方 Mm³（占 B1）", "天然库容 Mm³（占 B1 不挖）", "挖方 Mm³", "挖后库容 Mm³（占 B1）", "水位 ft（水头占 B1）", "GWh"], rows))
h.append('<p>规律：B1 的 3.9 km 里，西半环（沿 3,700 线绕鼻尖的第 1–3 段，1.9 km，RCC 15–19 Mm³）换到的天然库容只有 7 Mm³，因为那里的库底（3,700 ft 的等高线附近）本来就浅、山鼻又是凸的；东半环围的是主沟流域——同样的墙高，沟里的水深得多。填方效率：东半环 1.9（天然）/ 3.7（挖 20 后），B1 0.9 / 2.1，纯沟谷坝 3.0 / 3.3。</p>')
h.append("<h2>6. 要注意的</h2>")
h.append(f'<ul><li>墙高上限取 130 m。不设上限时优化器会把 W1 往沟下游再推（185 m 高坝、更多库容），因为沟里每延伸一段的边际效率仍 > 1；上限是设计选择，不是地形限制。</li>'
         f'<li>W2、W3 约 {f(REC["rcc_toe_out_m"]):.0f} m 的下游坝趾在宗地东界外（与 B1 的第 5–8 段同一问题），要邻地或把墙内移（内移会少几 Mm³ 库容）。</li>'
         '<li>填方按 RCC 断面；换 CFRD 是 3.5 倍（W1 跨沟口一段 1.4:1 的堆石要 26 Mm³），所以这个方案实质上是 RCC 方案，或 W1 沟口段 RCC + 其余堆石的组合。按 CFRD 成本重搜的结果见附录。</li>'
         '<li>4 m 网格估算的填方比 1 m 精算高约 30%（曲线横轴偏保守）；库容估算与精算差 5–7%。</li>'
         '<li>再挖的碗形坑规则（坑壁 0.75:1、离坝 24 m）与笔记 14 相同；坑最深处在主沟里，挖到 3,430 ft 意味着在 W1 上游 24 m 处有 60 m 深的坑壁，坝基与坑壁的稳定要一起核。</li>'
         '<li>地质假设未变：沟口的坝基是 Calaveras 千枚岩/变质火山岩，没有钻孔。</li></ul>')
cf_path = os.path.join(OUT, "wall_curves_cfrd.csv")
if os.path.exists(cf_path):
    cc = [r for r in csv.DictReader(open(cf_path)) if r["crest_ft"].startswith("4050") and f(r["fill_est_Mm3"]) > 0.05]
    h.append("<h2>附录：按 CFRD 成本重搜（坝顶 4,050）</h2>")
    first = cc[0] if cc else None
    if first: h.append(f'<p>把边界成本换成 CFRD 断面（1.4:1 两面，坡地放大更大）再扫描：λ 大于 {f(first["lambda"]):g} 时一段墙都不值得修；λ = {f(first["lambda"]):g} 时出现的仍是东半环，估算填方 {f(first["fill_est_Mm3"]):.0f} Mm³（1 m 精算 CFRD 21 Mm³、同一组墙用 RCC 7 Mm³），边际效率 {first["marg_eff"]}——每 m³ 堆石换不到 0.5 m³ 库容。也就是说，这组墙段只在 RCC（或沟口 RCC + 其余堆石）下成立；全用堆石时"墙段组合"并不比全环 B1 更划算。</p>')
    h.append(table(["λ", "累计填方（CFRD 估算）Mm³", "累计天然库容 Mm³", "墙总长 m", "天然岸长 m", "边际效率"], [[f(r["lambda"]), f"{f(r['fill_est_Mm3']):.1f}", f"{f(r['storage_est_Mm3']):.1f}", f"{f(r['dam_len_m']):,.0f}", f"{f(r['free_len_m']):,.0f}", r["marg_eff"] or "—"] for r in cc]))
    if os.path.exists(os.path.join(OUT, "wall_curves_cfrd.png")): h.append(fig(os.path.join(OUT, "wall_curves_cfrd.png"), "图 5 CFRD 成本下的曲线。"))
h.append('<p class="small">复算：<code>python3 scripts/wall_optimizer.py</code>（λ 扫描 + 精算，约 12 min）、<code>python3 scripts/wall_optimizer.py --plan=4050:2.1,4050:2.4,4000:2.1,3950:1.6,3900:1.85</code>（指定方案的平面与分段表）、<code>python3 scripts/wall_brief.py</code>。笔记：notes/16_wall_combination.md。</p>')
open(os.path.join(DEL, "D_wall_combo_brief.html"), "w").write(doc("墙段自由组合（简版）", "\n".join(h)))
print("written", os.path.getsize(os.path.join(DEL, "D_wall_combo_brief.html")) // 1024, "KB")
