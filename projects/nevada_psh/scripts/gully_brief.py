#!/usr/bin/env python3
"""跨沟谷短坝坝址搜索的简版交付（HTML）：读 outputs/gully_stage1*.csv、gully_stage2.csv 与三张图，回答"能否用 1/10 的填方换一半库容"。"""
import base64, csv, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); OUT = os.path.join(ROOT, "outputs"); DEL = os.path.join(ROOT, "deliverables")
f = lambda v: float(v) if v not in ("", None) else float("nan")
B1 = {"fill": 29.7, "need": 26.7, "net_nocut": 26.5, "net": 61.4, "nwl": 4030.0, "head_m": (4030 - 2400) * 0.3048, "gwh": 70.6, "dam_len_km": 3.9}
s2 = [r for r in csv.DictReader(open(os.path.join(OUT, "gully_stage2.csv")))]
s1 = [r for r in csv.DictReader(open(os.path.join(OUT, "gully_stage1.csv")))]
s1low = [r for r in csv.DictReader(open(os.path.join(OUT, "gully_stage1_zmin3400.csv")))] if os.path.exists(os.path.join(OUT, "gully_stage1_zmin3400.csv")) else []
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
.box{background:#f3e5f5;border-left:4px solid #7b1fa2;padding:8px 14px;margin:10px 0} .ans{background:#fdecea;border-left:4px solid #d62728;padding:10px 14px;margin:10px 0}
.f{font-family:Menlo,Consolas,monospace;background:#f1f3f5;padding:1px 5px;border-radius:3px;font-size:12.5px} .small{font-size:12px;color:#555} li{margin:4px 0}"""
doc = lambda title, body: f'<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{esc(title)}</title><style>{CSS}</style></head><body>{body}</body></html>'

ok2 = [r for r in s2 if r["spill"] != "True" and f(r["gross_Mm3"]) >= 0.3]
best = max(ok2, key=lambda r: f(r["net_rcc_Mm3"]))
best_in = max([r for r in ok2 if r["abut_in_parcel"] == "True"], key=lambda r: f(r["net_rcc_Mm3"]))
best_eff = max([r for r in ok2 if f(r["net_rcc_Mm3"]) > 1], key=lambda r: f(r["net_rcc_Mm3"]) / f(r["fill_rcc_Mm3"]))
h = []
h.append("<h1>跨沟谷短坝的坝址搜索（简版）：地块内 3,700 ft 以上有没有一条短坝能围出上库</h1>")
h.append('<div class="box"><b>问题。</b>不围整圈，找天然沟谷或窄口，用一条短坝围出上库；团队初估一条 350–400 m 长、沟底填高 130 m 的坝，填方 3–5 Mm³（全环 B1 的 1/10）。'
         '<b>核心问题：</b>有没有一个位置，能用全环方案 1/10 的填方（≈3 Mm³），换到全环一半以上的库容（≥ 31 Mm³），落差又不明显损失？</div>')
h.append(f'<div class="ans"><b>答案：没有。</b>地块内 3,700 ft 以上所有沟谷，单条短坝能围出的天然库容最大 {f(best["net_rcc_Mm3"]):.1f} Mm³（坝址 S{best["site"]}，坝顶 {f(best["crest_ft"]):.0f} ft，坝长 {f(best["dam_len_m"]):.0f} m，坝高 {f(best["H_max_m"]):.0f} m，RCC 填方 {f(best["fill_rcc_Mm3"]):.1f} / CFRD {f(best["fill_cfrd_Mm3"]):.1f} Mm³，坝肩已出宗地东界）；'
         f'坝肩全在宗地内的最好一处是 S{best_in["site"]}，坝顶 {f(best_in["crest_ft"]):.0f}，{f(best_in["net_rcc_Mm3"]):.1f} Mm³。这是 B1 库容（61.4 Mm³）的 {100*f(best["net_rcc_Mm3"])/B1["net"]:.0f}%、不挖时 B1（26.5 Mm³）的 {100*f(best["net_rcc_Mm3"])/B1["net_nocut"]:.0f}%。'
         f'填方确实只有 B1 的 {100*f(best["fill_rcc_Mm3"])/B1["fill"]:.0f}%，落差也几乎不损失（水位 {f(best["nwl_ft"]):.0f} 对 4,030 ft），但库容只有 1/10，不是 1/2。把沟口下移到 3,400 ft 也一样（最大 5–6 Mm³）。原因见第 2 节：这是地形决定的，与坝型和坝址选法无关。</div>')
h.append("<h2>1. 我做了什么</h2>")
h.append('<ul><li>在 4 m 分辨率的 LiDAR 地形上做水文分析（D8 汇流），找出地块内所有汇水面积 ≥ 1 ha 的沟谷线；沿沟谷每 30 m 放一个候选坝址，坝轴垂直于沟谷（±30° 内取坝长最短的方向）。</li>'
         '<li>每个坝址试坝顶 3,800–4,120 ft（每 50 ft 一档）：沿坝轴两侧找地面 ≥ 坝顶的坝肩 → 坝长与坝高；用"洪水填充"算坝后水位以下、与坝相连、不越过坝轴的淹没区 → 库容；若淹没区漫到窗口边缘（从别的鞍部流走）判为关不住。</li>'
         f'<li>{len(s1)} 个坝址×坝顶（沟底 ≥ 3,650 ft）通过筛选；再把沟底下限放到 3,400 ft 又试了 {len(s1low)} 个。入围的 6 个坝址在 1 m 地形上逐站放坝（CFRD 1.4/1.4 与 RCC 直立/0.8，同交付物 1），精算填方、上游楔、库容，并试了"坝后再挖碗形坑"。</li>'
         '<li>坝型没有锁定：两种断面的填方都给出；库容按 RCC 断面的净值报（CFRD 上游坡要再扣 30–40%）。</li></ul>')
h.append(fig(os.path.join(OUT, "gully_thalwegs.png"), "图 1 沟谷线（蓝）与全部候选坝址在各自最大库容坝顶下的坝轴（红），数字是毛库容 Mm³。3,700 ft 以上只有两条沟：东南角向西北流的主沟（汇水 27 ha）和它北边的小沟（7–9 ha），都是浅而陡的沟头，没有峡谷式的窄口。"))
h.append("<h2>2. 为什么围不出来：两个手算上限</h2>")
h.append('<ol><li><b>库容 ≤ 淹没面积 × 平均水深。</b>一条坝只能淹它上游的汇水区里低于水位的部分。主沟在坝址处的汇水面积 27 ha（67 acre），小沟 7–9 ha；沟谷纵坡 13–24°（tanθ 0.24–0.32），水位比沟底高 H 时水面只往上游伸 <span class="f">H/tanθ</span> ≈ 3–4 H，坝高 100 m 时水面只有 300–400 m 长、20–40 acre，平均水深 30–40 m → 3–6 Mm³。'
         ' V 形谷的量级公式 <span class="f">V ≈ H³ / (3·tanα·tanθ)</span>（α 谷坡、θ 纵坡）：H 100 m、tanα 0.5、tanθ 0.28 → 2.4 Mm³，与模型同量级。</li>'
         '<li><b>水面必须落在"比水位高的地"围出来的范围里。</b>地块内高于 3,800 ft 的地只有 98 acre，高于 3,900 ft 56 acre，高于 4,000 ft 27 acre，而且是一个凸出的山鼻，不是盆地——四周比中间低，没有天然的闭合圈。'
         ' 全环 B1 之所以有 61 Mm³，是用 3.9 km 的坝把 3,700 ft 的圈"抬"了 100 m，再在圈里挖 36 Mm³；短坝只封一个沟口，圈没有抬起来。'
         ' 要从 30 acre 的水面拿到 31 Mm³，平均要挖 250 m 深——0.75:1 的坑壁在 300 m 宽的沟里 110 m 就相遇了，做不到。</li></ol>')
h.append("<h2>3. 入围坝址的结果</h2>")
rows = []
for r in sorted(ok2, key=lambda r: (-f(r["net_rcc_Mm3"]))):
    rows.append([f"S{r['site']}", f"{f(r['z_floor_ft']):,.0f}", f"{f(r['crest_ft']):,.0f} / {f(r['nwl_ft']):,.0f}", f"{f(r['dam_len_m']):.0f}", f"{f(r['H_max_m']):.0f}",
                 "是" if r["abut_in_parcel"] == "True" else "否", f"{f(r['fill_rcc_Mm3']):.2f} / {f(r['fill_cfrd_Mm3']):.2f}", f"{f(r['gross_Mm3']):.2f}", f"{f(r['net_rcc_Mm3']):.2f} / {f(r['net_cfrd_Mm3']):.2f}",
                 f"{f(r['water_acre']):.0f}（出界 {f(r['water_out_acre']):.0f}）", f"{f(r['net_rcc_Mm3'])/f(r['fill_rcc_Mm3']):.1f}", f"{f(r['deep_net_rcc_Mm3']):.1f}（挖 {f(r['deep_cut_Mm3']):.1f}）" if r.get("deep_net_rcc_Mm3") else "—",
                 f"{f(r['net_rcc_GWh']):.1f}", f"{100*f(r['net_rcc_Mm3'])/B1['net']:.0f}% / {100*f(r['fill_rcc_Mm3'])/B1['fill']:.0f}% / {100*(f(r['nwl_ft'])-2400)/(B1['nwl']-2400):.0f}%"])
h.append(table(["坝址", "沟底 ft", "坝顶 / 水位 ft", "坝长 m", "坝高 m", "坝肩在界内", "填方 RCC / CFRD Mm³", "毛库容 Mm³", "净库容 RCC / CFRD Mm³", "水面 acre", "库容÷填方(RCC)", "坝后再挖：净库容 Mm³", "GWh", "占 B1：库容 / 填方 / 落差"], rows))
h.append('<p class="small">坝肩出界 = 坝顶高程在宗地东界内找不到那么高的地面，坝要延到界外。"坝后再挖"按笔记 14 的碗形坑（坑顶线在上游坝趾内 20 m，坑壁 0.75:1）：因为水面只有 10–40 acre、又被坝趾带占去一半，坑很小，只能多给 0.1–0.6 Mm³。GWh 按下库 2,400 ft、η 0.85。数据：outputs/gully_stage2.csv。</p>')
h.append(fig(os.path.join(OUT, "gully_sites_plan.png"), "图 2 入围坝址（各取库容最大的坝顶）：粗线坝轴、同色填充淹没区。所有淹没区都挤在东南角的两条沟头里。"))
h.append(fig(os.path.join(OUT, "gully_fill_vs_storage.png"), "图 3 每个坝址×坝顶的填方（横轴，对数）与净库容（纵轴），与 B1 对照。竖虚线是 B1 填方的 1/10，横虚线是 B1 库容的 1/2：右上角那个象限里没有任何沟谷坝。"))
h.append("<h2>4. 与团队初估的对照</h2>")
h.append('<ul><li>"350–400 m 长、沟底填高 130 m"：地块里没有这样的谷。要填高 100 m 以上（坝顶 4,000–4,050），坝长是 670–840 m，因为沟谷是开口的浅沟，两岸要走很远才到坝顶高程；坝长 350–400 m 时坝高只有 40–56 m，库容 0.4–2.6 Mm³。</li>'
         '<li>"填方 3–5 Mm³"：对，CFRD 断面在 H 90–100 m 时确实是 3–4.5 Mm³，RCC 1.0–1.6 Mm³。填方估对了，库容估错了一个数量级——因为库容随沟谷纵坡 1/tanθ 变化，这里的沟头坡度 13–24°，水面伸不远。</li>'
         f'<li>库容/填方效率：沟谷坝 {f(best_eff["net_rcc_Mm3"])/f(best_eff["fill_rcc_Mm3"]):.1f}（S{best_eff["site"]} 坝顶 {f(best_eff["crest_ft"]):.0f}），B1 不挖 0.9、挖后 2.1。效率高但绝对量小，是"小水库"而不是"半个大水库"。</li></ul>')
h.append("<h2>5. 如果目标是少填方，真正可选的是这几条</h2>")
h.append('<ul><li><b>把全环的坝顶压低，而不是换坝址。</b>组合 III 坝顶 3,900（交付物 1 的 A）：填方 20 Mm³（RCC 只有 5）、库容 25 Mm³、落差 −3%——东界和南台地已经不需要坝，坝只剩 3.1 km。</li>'
         '<li><b>台地顶的小环</b>（笔记 13 解读 C，4,000 ft 环 25 acre）：坝短、填方小，但库容也只有 4–5 Mm³，和沟谷坝一个量级。</li>'
         f'<li><b>沟谷坝当一期小库</b>：S29 坝顶 3,950–4,000（RCC 1.5–2.1 Mm³、库容 4–6 Mm³、水位 3,930–3,980），坝肩在 4,000 时出东界，需要邻地；或 S{best_in["site"]} 坝顶 {f(best_in["crest_ft"]):.0f}（全在界内，{f(best_in["net_rcc_Mm3"]):.1f} Mm³）。</li>'
         '<li>要拿到 30 Mm³ 以上，只有把圈抬起来（环形坝）或者把圈扩到宗地外的台地（笔记 13 解读 D，要用地）。</li></ul>')
h.append('<p class="small">复算：<code>python3 scripts/gully_dam_sites.py</code>（默认沟底 ≥ 3,650 ft；<code>--zmin=3400</code> 放宽到 3,400），<code>python3 scripts/gully_brief.py</code> 生成本文件。笔记：notes/15_gully_dam_sites.md。</p>')
open(os.path.join(DEL, "D_gully_sites_brief.html"), "w").write(doc("跨沟谷短坝坝址搜索（简版）", "\n".join(h)))
print("written", os.path.getsize(os.path.join(DEL, "D_gully_sites_brief.html")) // 1024, "KB; best", best["site"], best["crest_ft"], best["net_rcc_Mm3"], "| best_in", best_in["site"], best_in["crest_ft"], best_in["net_rcc_Mm3"], "| eff", best_eff["site"], best_eff["crest_ft"])
