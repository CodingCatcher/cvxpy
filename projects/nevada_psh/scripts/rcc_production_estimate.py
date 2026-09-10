#!/usr/bin/env python3
"""RCC 坝段的方量、日浇筑量与拌合站规模估算（资料 5：USBR RCC 手册 2017 + USACE EM 1110-2-2006）。

思路：层面热缝要求每天至少 2 层（单班）或 3 层（双班/连续）。一层的面积 ≈ 段长 × 该高程的断面宽，
最底层最大。于是"保住热缝"所需的日浇筑量 = 层数/天 × 层厚 × 最大层面积，
拌合站持续产能 = 日浇筑量 / (班时 × 效率)。再按 USACE 的生产率分级（小 50–230、中 230–460、大 460–1,000 yd³/h）
和 USBR 附录 B 的方量-单价带给出定位。这只是把两份资料的数字连起来，不做选型。

断面：上游面垂直，坝顶宽 b，下游坡 m:1（水平:竖直），高 H 的段落基底宽 B = b + m·H。
段落用 "L:H" 或 "L:H1:H2"（沿段长线性变化）表示，单位米。

用法：
  python3 rcc_production_estimate.py                         # 默认情况表
  python3 rcc_production_estimate.py --seg 300:20:30 --seg 150:25 --m 0.8 --b 6.1 --lifts 2
"""
import argparse

YD3_PER_M3 = 1.30795
LIFT = 0.30            # m，压实层厚 12 in
CLASSES = [(230, "小型 (50–230 yd³/h)"), (460, "中型 (230–460 yd³/h)"), (1000, "大型 (460–1,000 yd³/h)")]
COST_BANDS = [  # (最大方量 yd³, 附录 B 单价带，不含胶材，未折现)
    (20000, "$95–170/yd³（0.6–1.9 万 yd³ 级，2000s）"),
    (70000, "$30–46/yd³（1.5–6.3 万 yd³ 级，1990s–2000）"),
    (1000000, "USBR 案例无此量级（6.3 万–100 万 yd³ 之间），需另找数据"),
    (float("inf"), "$10–24/yd³（>100 万 yd³ 级，1980s）"),
]


def parse_seg(s):
    p = [float(x) for x in s.split(":")]
    if len(p) == 2:
        return p[0], p[1], p[1]
    if len(p) == 3:
        return p[0], p[1], p[2]
    raise argparse.ArgumentTypeError("段落格式 L:H 或 L:H1:H2")


def segment(L, H1, H2, b, m):
    """返回 (体积 m³, 最大层面积 m², 平均基底宽, 最大高度, 层数)。
    高度沿段长线性：面积/m = b·H + m·H²/2；沿长度积分。"""
    n = 200
    vol = 0.0
    lift_area = 0.0
    for i in range(n):
        H = H1 + (H2 - H1) * (i + 0.5) / n
        dx = L / n
        vol += (b * H + 0.5 * m * H * H) * dx
        lift_area += (b + m * H) * dx        # 最底层（各处基底宽之和）
    Hmax = max(H1, H2)
    return vol, lift_area, lift_area / L, Hmax, int(round(Hmax / LIFT))


def classify(rate):
    for cap, name in CLASSES:
        if rate <= cap:
            return name
    return "超出 USACE 分级（>1,000 yd³/h）"


def cost_band(vol_yd3):
    for cap, name in COST_BANDS:
        if vol_yd3 <= cap:
            return name
    return COST_BANDS[-1][1]


def report(segs, b, m, lifts_per_day, shift_h, eff):
    print(f"# 断面：上游垂直，坝顶宽 {b} m，下游坡 {m}:1；层厚 {LIFT} m；热缝要求 {lifts_per_day} 层/天；"
          f"班时 {shift_h} h，效率 {eff}")
    hdr = (f"{'段落':>14} {'B均(m)':>7} {'方量(m³)':>10} {'方量(yd³)':>10} {'底层面积(m²)':>12} "
           f"{'日浇(yd³/d)':>11} {'持续率(yd³/h)':>13} {'层数':>4} {'工期(d)':>7}  分级 / 单价带")
    print(hdr)
    print("-" * len(hdr))
    tot_vol = 0.0
    for (L, H1, H2) in segs:
        vol, la, Bmean, Hmax, nl = segment(L, H1, H2, b, m)
        vol_yd = vol * YD3_PER_M3
        day_m3 = lifts_per_day * LIFT * la
        day_yd = day_m3 * YD3_PER_M3
        rate = day_yd / (shift_h * eff)
        days = nl / lifts_per_day           # 若按热缝节奏，层数/每天层数
        tot_vol += vol
        name = f"{L:.0f}m×{H1:.0f}" + (f"–{H2:.0f}" if H2 != H1 else "") + "m"
        print(f"{name:>14} {Bmean:7.1f} {vol:10.0f} {vol_yd:10.0f} {la:12.0f} {day_yd:11.0f} {rate:13.0f} {nl:4d} {days:7.0f}  "
              f"{classify(rate)} / {cost_band(vol_yd)}")
    print(f"\n合计方量 {tot_vol:.0f} m³ = {tot_vol*YD3_PER_M3:.0f} yd³ → 整体单价带：{cost_band(tot_vol*YD3_PER_M3)}")
    print("# 注：'日浇'是为保住热缝在最底层必须达到的量，越往上层面积越小、要求越低；")
    print("#     '工期'仅按层数/日层数计，不含试验段、面板、缝、廊道、天气；")
    print("#     若做不到该持续率，层面退化为冷缝/施工缝 → 全靠垫层砂浆，粘聚力与扬压力假设随之改变（笔记 §8.2、§8.4）。")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seg", type=parse_seg, action="append", help="段落 L:H 或 L:H1:H2（米），可多次")
    ap.add_argument("--m", type=float, default=0.8, help="下游坡（水平:竖直）")
    ap.add_argument("--b", type=float, default=6.1, help="坝顶宽 m（USBR 施工最小 20 ft = 6.1 m）")
    ap.add_argument("--lifts", type=int, default=2, help="每天层数（2 单班 / 3 双班）")
    ap.add_argument("--shift", type=float, default=16.0, help="每天浇筑小时数")
    ap.add_argument("--eff", type=float, default=0.8, help="拌合站效率（规范脚注：按 80% 计）")
    a = ap.parse_args()
    segs = a.seg or [(100, 20, 20), (100, 30, 30), (300, 20, 20), (300, 30, 30), (300, 20, 30), (600, 25, 25)]
    report(segs, a.b, a.m, a.lifts, a.shift, a.eff)


if __name__ == "__main__":
    main()
