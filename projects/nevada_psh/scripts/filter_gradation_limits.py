#!/usr/bin/env python3
"""USBR DS-13 第 5 章 (2011) 反滤级配控制点计算（第 2–8 步）。

输入底土级配（筛孔 mm → 过筛百分数），输出：
  再分级后的曲线、底土类别、D15F 上限（点 A）、D15F 下限（点 B）、
  D5F/D100F 限（点 I/J）、D90F 上限（点 K）、带宽限（35 百分点），
  以及附录 C 的"水平滑杆"锚点 C/D/E/F。

只做第 2–8 步；第 1 步（选底土：平均 / 细侧包络 / 类别选择）和第 9 步
（在控制点内画带）需要工程判断，脚本不替代。结果按笔记
notes/08_usbr_ds13_ch10_construction_ch5_filters.md 第 16 节。

用法：
  python filter_gradation_limits.py --demo          # 附录 C 两个界面的算例复算
  python filter_gradation_limits.py --grad "75:100,37.5:85.7,19:74.6,9.5:65.9,4.75:57.9,2.36:54.6,1.18:49,0.6:42.6,0.3:32.2,0.15:19.8,0.075:13" [--dispersive]
"""
import argparse
import math

NO4 = 4.75      # mm
NO200 = 0.075   # mm


def interp_d(grad, pct):
    """从级配（按粒径升序的 (d, p) 列表）按对数粒径线性内插 D_pct。"""
    pts = sorted(grad)
    if pct <= pts[0][1]:
        return pts[0][0]
    if pct >= pts[-1][1]:
        return pts[-1][0]
    for (d0, p0), (d1, p1) in zip(pts, pts[1:]):
        if p0 <= pct <= p1:
            if p1 == p0:
                return d0
            f = (pct - p0) / (p1 - p0)
            return 10 ** (math.log10(d0) + f * (math.log10(d1) - math.log10(d0)))
    return pts[-1][0]


def pct_passing(grad, d):
    """按对数粒径内插粒径 d 的过筛百分数。"""
    pts = sorted(grad)
    if d <= pts[0][0]:
        return pts[0][1]
    if d >= pts[-1][0]:
        return pts[-1][1]
    for (d0, p0), (d1, p1) in zip(pts, pts[1:]):
        if d0 <= d <= d1:
            f = (math.log10(d) - math.log10(d0)) / (math.log10(d1) - math.log10(d0))
            return p0 + f * (p1 - p0)
    return pts[-1][1]


def coefficients(grad):
    d10, d30, d60 = (interp_d(grad, p) for p in (10, 30, 60))
    cu = d60 / d10 if d10 > 0 else float("inf")
    cz = d30 ** 2 / (d60 * d10) if d10 > 0 and d60 > 0 else float("nan")
    return d10, d30, d60, cu, cz


def needs_regrading(grad):
    """图 5.4.2-3：含砾石且不同时满足三条件则再分级。"""
    p4 = pct_passing(grad, NO4)
    if p4 >= 100 - 1e-9:
        return False, "无砾石（2a）"
    fines = pct_passing(grad, NO200)
    _, _, _, cu, cz = coefficients(grad)
    cond_fines = fines < 15
    # 图 5.4.2-3 原文："not broadly graded (i.e., Cu not > 6 and Cz not between 1 and 3)"。
    # 按字面：Cu ≤ 6 且 Cz 不在 1–3 才算"非宽级配"。附录 C 把 Cu=398–811、Cz=0.64 的土也称宽级配，
    # 即 Cu 是决定性的；这里按字面实现（偏保守：更常再分级 → 反滤更细）。
    cond_broad = (cu <= 6) and not (1 <= cz <= 3)
    # 间断级配的判定：任一段斜率比 4x 线更平（D_i < 4·D_{i-1} 上过筛率跳变 >35 点视为间断）——用手册"4x 线"近似
    cond_gap = not is_gap_graded(grad)
    if cond_fines and cond_broad and cond_gap:
        return False, "含砾石但细料<15%、非间断、非宽级配（2b 全满足）"
    why = []
    if not cond_fines:
        why.append(f"细料 {fines:.1f}% ≥15%")
    if not cond_broad:
        why.append(f"宽级配 Cu={cu:.1f} Cz={cz:.2f}")
    if not cond_gap:
        why.append("疑似间断级配")
    return True, "；".join(why)


def is_gap_graded(grad):
    """粗略判据：相邻两筛粒径比 ≤4 而过筛率差 ≥35 点 → 视为间断（手册无定量定义，此为保守近似）。"""
    pts = sorted(grad)
    for (d0, p0), (d1, p1) in zip(pts, pts[1:]):
        if d1 / d0 <= 4.0 and (p1 - p0) >= 35:
            return True
    return False


def regrade(grad):
    """第 3 步：以过 4 号筛百分数归一。返回只含 ≤4.75 mm 的曲线。"""
    p4 = pct_passing(grad, NO4)
    k = 100.0 / p4
    out = [(d, min(100.0, p * k)) for d, p in sorted(grad) if d <= NO4 + 1e-9]
    if not any(abs(d - NO4) < 1e-9 for d, _ in out):
        out.append((NO4, 100.0))
    return sorted(out)


def category(fines):
    if fines > 85:
        return 1
    if fines >= 40:
        return 2
    if fines >= 15:
        return 3
    return 4


def d15f_max(cat, d85b, fines, dispersive=False):
    """表 5.4.4-1，点 A。"""
    if cat == 1:
        return max((6.5 if dispersive else 9.0) * d85b, 0.2)
    if cat == 2:
        return 0.5 if dispersive else 0.7
    if cat == 3:
        if dispersive:
            return 0.5
        base = max(4.0 * d85b, 0.7)
        return (40.0 - fines) / (40.0 - 15.0) * (base - 0.7) + 0.7
    return 4.0 * d85b


def d15f_min(d15b):
    """第 6 步，点 B。"""
    return max(5.0 * d15b, 0.1)


def d90f_max(d10f):
    """表 5.4.6-2，点 K。"""
    for lim, val in ((0.5, 20), (1.0, 25), (2.0, 30), (5.0, 40), (10.0, 50)):
        if d10f < lim:
            return val
    return 60


def design(grad, dispersive=False, label="", force_no_regrade=False):
    res = {"label": label}
    d10, d30, d60, cu, cz = coefficients(grad)
    res["input"] = dict(fines=pct_passing(grad, NO200), p4=pct_passing(grad, NO4), cu=cu, cz=cz)
    rg, why = needs_regrading(grad)
    if force_no_regrade:
        rg, why = False, "用户指定不再分级"
    res["regrade"] = (rg, why)
    g = regrade(grad) if rg else sorted(grad)
    fines = pct_passing(g, NO200)
    cat = category(fines)
    d85b, d15b = interp_d(g, 85), interp_d(g, 15)
    A = d15f_max(cat, d85b, fines, dispersive)
    B = d15f_min(d15b)
    C, D = 0.7 * A, max(0.7 * B, 0.075)
    E, F = 6.0 * C, 2.0 * D
    res.update(dict(regraded=g, fines=fines, cat=cat, d85b=d85b, d15b=d15b,
                    A=A, B=B, I=0.075, J=51.0, K=d90f_max(D), C=C, D=D, E=E, F=F,
                    band_pts=35, dispersive=dispersive))
    return res


def report(r):
    print(f"=== {r['label']} ===")
    i = r["input"]
    print(f"输入：过200号 {i['fines']:.1f}%  过4号 {i['p4']:.1f}%  Cu={i['cu']:.1f}  Cz={i['cz']:.2f}")
    rg, why = r["regrade"]
    print(f"再分级：{'是' if rg else '否'}（{why}）")
    if rg:
        print("  再分级后：" + ", ".join(f"{d:g}:{p:.1f}" for d, p in r["regraded"]))
    print(f"底土：过200号 {r['fines']:.1f}% → {r['cat']} 类；D85B={r['d85b']:.3f} mm  D15B={r['d15b']:.4f} mm"
          + ("（分散性）" if r["dispersive"] else ""))
    print(f"点A D15F max = {r['A']:.3f} mm    点B D15F min = {r['B']:.3f} mm")
    print(f"点I D5F min = {r['I']} mm   点J D100F max = {r['J']} mm   点K D90F max = {r['K']} mm（按 D10F min={r['D']:.3f}）")
    print(f"带宽 ≤{r['band_pts']} 百分点；水平滑杆锚点 C={r['C']:.3f} D={r['D']:.3f} E={r['E']:.2f} F={r['F']:.3f} mm")
    if r["A"] < r["B"]:
        print("!! A < B：颗粒保持上限低于渗透性下限，无解——需回到第 1 步（分级/分区/两级）")
    print()


def parse_grad(s):
    out = []
    for tok in s.split(","):
        d, p = tok.split(":")
        out.append((float(d), float(p)))
    return sorted(out)


def demo():
    # 附录 C 界面 1：既有坝心墙细侧包络（再分级后 3 类，A=35%，D85B=1.71，D15B=0.005）
    # 直接给出再分级后的量以复算控制点（原文未给完整曲线）。
    print("附录 C 界面 1（3 类，A=35%，D85B=1.71 mm，D15B=0.005 mm）")
    A = d15f_max(3, 1.71, 35.0)
    B = d15f_min(0.005)
    print("  注：手册印 A=1.98，但按其自身公式 0.2×(4×1.71−0.7)+0.7 = 1.93；后续 C/E 与 1.98 一致，属手册算术误差。")
    print(f"  A = {A:.2f} mm（手册 1.98）  B = {B:.2f} mm（手册 0.10）  "
          f"C={0.7*A:.2f}（1.39） D={max(0.7*B,0.075):.3f}（0.075） E={6*0.7*A:.2f}（8.34） F={2*0.075:.2f}（0.15）  "
          f"K={d90f_max(0.075)}（20）")
    print("附录 C 界面 2（C33 砂作底土：4 类，D85B=1.18，D15B=0.18）")
    A2 = d15f_max(4, 1.18, 2.0)
    B2 = d15f_min(0.18)
    print(f"  A = {A2:.2f} mm（4.72）  B = {B2:.2f} mm（0.90）  C={0.7*A2:.2f}（3.30） D={0.7*B2:.2f}（0.63） "
          f"E={6*0.7*A2:.1f}（19.8） F={2*0.7*B2:.2f}（1.26）  K={d90f_max(0.63)}（25）")
    print(f"  放宽 9·D85B = {9*1.18:.2f}（10.62）；No.467 D10F=5.5 → K={d90f_max(5.5)}（50）")
    print()
    # 图 5.4.2-4 再分级算例
    g = parse_grad("75:100,37.5:85.7,19:74.6,9.5:65.9,4.75:57.9,2.36:54.6,1.18:49,0.6:42.6,0.3:32.2,0.15:19.8,0.075:13,0.03:9.9,0.01:5.4,0.003:2.9,0.001:1.6")
    r = design(g, label="图 5.4.2-4 算例（再分级后过200号应为 22.5%）")
    report(r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grad", help="底土级配 'd_mm:pct,...'")
    ap.add_argument("--dispersive", action="store_true")
    ap.add_argument("--no-regrade", action="store_true", help="强制不再分级（均粒砾石第二级设计）")
    ap.add_argument("--label", default="底土")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    if a.demo or not a.grad:
        demo()
        return
    report(design(parse_grad(a.grad), a.dispersive, a.label, a.no_regrade))


if __name__ == "__main__":
    main()
