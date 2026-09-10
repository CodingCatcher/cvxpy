#!/usr/bin/env python3
"""重力坝断面稳定筛查：把 EM 1110-2-2200 与 FERC 第 3 章的准则写成可复算的公式。

用途：对给定的三角形/梯形重力坝断面（上游面垂直、坝顶宽 b、下游坡 m:1），
在"正常运行（usual）"工况下计算：
  - 自重、静水压力、扬压力（可选排水线位置与排水有效率 E）
  - 基底裂缝迭代（裂缝内全水头；裂缝越过排水线则 E=0）
  - 合力位置、基底压应力、剪摩安全系数 FS = (N·tanφ + c·L_uncracked) / ΣH
  - 对照：EM 2200 Table 4-1（usual：中 1/3、FS 2.0）
          FERC Table 2（高/显著危害有粘聚力 3.0；低危害 2.0）
          FERC Table 2A（无粘聚力 1.5）
  - 反算满足各 FS 目标且静力不开裂的最缓下游坡 m_min
可选 --kh 做地震系数法初筛（仅 USACE 允许作初筛；FERC 不接受拟静力，故只作量级参考）。

不做的事：不选坝型、不选断面。只是把两份资料的数字变成可以逐段套用的工具。

用法：
  python3 gravity_section_check.py                       # 打印 H=20/25/30 m 的情况表
  python3 gravity_section_check.py --H 25 --m 0.8 --phi 40 --E 0.5 --detail
  python3 gravity_section_check.py --freeboard 0          # 水面到坝顶
"""
import argparse
import math

GAMMA_W = 9.81   # kN/m3
GAMMA_C = 23.5   # kN/m3 (~150 pcf)
M_SCAN_MIN = 0.50   # 反算最缓坡时的扫描下限；EM 2200 的实用范围是 0.7–0.8，比它陡只说明"稳定不是控制因素"
TAN_BETA = math.tan(math.radians(24.0))   # 本地块山坡 24°，用于水平基底的坝踵超挖估算


def section(H, m, b=0.0, t=0.0):
    """返回断面多边形顶点（逆时针不要求）与基底宽 B。
    坝踵 (0,0)；上游面垂直到 (0,H)；坝顶到 (b,H)；下游面先垂直 t 再按 m:1 到坝趾 (B,0)。"""
    B = b + m * (H - t)
    pts = [(0.0, 0.0), (0.0, H), (b, H), (b, H - t), (B, 0.0)]
    return pts, B


def area_centroid(pts):
    a = cx = cy = 0.0
    n = len(pts)
    for i in range(n):
        x0, y0 = pts[i]
        x1, y1 = pts[(i + 1) % n]
        cr = x0 * y1 - x1 * y0
        a += cr
        cx += (x0 + x1) * cr
        cy += (y0 + y1) * cr
    a *= 0.5
    return abs(a), cx / (6 * a), cy / (6 * a)


def uplift(hw, B, Lc, xd, E):
    """扬压力（kN/m）及其对坝踵的力臂（m）。
    水头折线：坝踵→裂缝尖端 Lc 全水头 hw；若有排水且 Lc < xd：在 xd 处水头 (1-E)·hw；坝趾 0（无尾水）。
    裂缝越过排水线 → E 视为 0（EM 2200 Fig 3-5 / FERC Fig 3）。"""
    pts = [(0.0, hw), (Lc, hw)]
    if xd is not None and Lc < xd and E > 0:
        pts.append((xd, (1 - E) * hw))
    pts.append((B, 0.0))
    U = M = 0.0
    for (x0, h0), (x1, h1) in zip(pts[:-1], pts[1:]):
        if x1 <= x0:
            continue
        w = x1 - x0
        # 梯形分成矩形 + 三角形
        rect = h0 * w if h0 <= h1 else h1 * w
        tri_h = abs(h1 - h0)
        tri = 0.5 * tri_h * w
        x_rect = x0 + w / 2
        x_tri = x0 + (2 * w / 3 if h1 > h0 else w / 3)
        U += GAMMA_W * (rect + tri)
        M += GAMMA_W * (rect * x_rect + tri * x_tri)
    return U, (M / U if U > 0 else 0.0)


def analyse(H, m, hw, phi_deg, c=0.0, E=0.0, xd_frac=None, b=0.0, t=0.0, kh=0.0):
    """返回字典：合力位置、裂缝长度、基底应力、FS 等。"""
    pts, B = section(H, m, b, t)
    A, xg, yg = area_centroid(pts)
    W = GAMMA_C * A
    P = 0.5 * GAMMA_W * hw ** 2          # 静水压力，作用点 hw/3
    xd = None if (xd_frac is None or E <= 0) else xd_frac * B

    # 地震系数法（仅初筛）：坝体惯性 kh·W 在形心高度；Westergaard 水动力 7/12·kh·γw·hw²，作用点 0.4hw
    Fi = kh * W
    Fh = (7.0 / 12.0) * kh * GAMMA_W * hw ** 2

    tan_phi = math.tan(math.radians(phi_deg))

    def solve(Lc):
        U, xu = uplift(hw, B, Lc, xd, E)
        N = W - U
        # 对坝踵取矩求合力落点：自重使合力靠坝踵，水压/惯性（向下游）与扬压力使合力向坝趾移动
        if N <= 0:
            return None
        xr = (W * xg - U * xu + P * hw / 3 + Fi * yg + Fh * 0.4 * hw) / N   # 合力到坝踵距离
        Bp = B - Lc                   # 未裂宽
        e = xr - (Lc + Bp / 2)        # 相对未裂段中心的偏心，向坝趾为正
        s_tip = (N / Bp) * (1 - 6 * e / Bp)   # 裂缝尖端（或坝踵）应力，压为正
        s_toe = (N / Bp) * (1 + 6 * e / Bp)
        return dict(U=U, N=N, xr=xr, Bp=Bp, s_tip=s_tip, s_toe=s_toe)

    r0 = solve(0.0)
    if r0 is None:
        return dict(B=B, W=W, P=P, unstable=True)
    Lc = 0.0
    r = r0
    if r0["s_tip"] < 0:
        # 二分找 s_tip = 0 的裂缝长度
        lo, hi = 0.0, B * 0.999
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            rm = solve(mid)
            if rm is None or rm["s_tip"] < 0:
                lo = mid
            else:
                hi = mid
        Lc = hi
        r = solve(Lc)
        if r is None:
            return dict(B=B, W=W, P=P, unstable=True)
    H_sum = P + Fi + Fh
    FS = (r["N"] * tan_phi + c * r["Bp"]) / H_sum
    return dict(
        B=B, W=W, P=P, U=r["U"], N=r["N"], xr=r["xr"], Lc=Lc, Bp=r["Bp"],
        s_heel=r["s_tip"], s_toe=r["s_toe"], FS=FS,
        in_middle_third=(B / 3 <= r["xr"] <= 2 * B / 3), unstable=(r["xr"] <= 0 or r["xr"] >= B),
    )


def min_slope(H, hw, phi, c, E, xd_frac, target_FS, b=0.0, t=0.0, require_no_crack=True, require_third=True):
    """扫描下游坡 m，找满足 FS ≥ target 且（可选）不开裂、合力在中 1/3 的最小 m。"""
    m = M_SCAN_MIN
    while m <= 2.50:
        r = analyse(H, m, hw, phi, c, E, xd_frac, b, t)
        ok = (not r.get("unstable")) and r["FS"] >= target_FS
        if ok and require_no_crack:
            ok = r["Lc"] <= 1e-6
        if ok and require_third:
            ok = r["in_middle_third"]
        if ok:
            return m
        m += 0.01
    return float("nan")


def print_detail(H, m, hw, phi, c, E, xd_frac, b, t, kh):
    r = analyse(H, m, hw, phi, c, E, xd_frac, b, t, kh)
    print(f"断面：H={H} m, 水深 hw={hw:.2f} m, 下游坡 {m}:1, 坝顶宽 {b} m, 基底宽 B={r['B']:.2f} m")
    print(f"参数：φ={phi}°, c={c} kPa, 排水有效率 E={E}, 排水线位置 {xd_frac if E>0 else '无'}·B, kh={kh}")
    if r.get("unstable"):
        print("  → 合力出基底 / 净竖向力 ≤ 0：不稳定")
        return
    print(f"  自重 W={r['W']:.0f} kN/m  水压 P={r['P']:.0f} kN/m  扬压 U={r['U']:.0f} kN/m  净竖向 N={r['N']:.0f} kN/m")
    print(f"  合力距坝踵 {r['xr']:.2f} m（中1/3：{r['B']/3:.2f}–{2*r['B']/3:.2f}）  裂缝长 Lc={r['Lc']:.2f} m")
    print(f"  基底应力：坝踵/裂缝尖端 {r['s_heel']:.0f} kPa，坝趾 {r['s_toe']:.0f} kPa")
    print(f"  剪摩 FS={r['FS']:.2f}")
    checks = [
        ("EM 2200 usual：中1/3 且 FS≥2.0", r["in_middle_third"] and r["FS"] >= 2.0),
        ("FERC Table 2 高/显著危害（有 c）：FS≥3.0 且不开裂", r["FS"] >= 3.0 and r["Lc"] <= 1e-6),
        ("FERC Table 2 低危害：FS≥2.0 且不开裂", r["FS"] >= 2.0 and r["Lc"] <= 1e-6),
        ("FERC Table 2A 无 c：FS≥1.5 且不开裂", r["FS"] >= 1.5 and r["Lc"] <= 1e-6),
    ]
    for name, ok in checks:
        print(f"  [{'OK' if ok else '--'}] {name}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--H", type=float, help="坝高 m（单个断面详算）")
    ap.add_argument("--m", type=float, default=0.8, help="下游坡 m:1（水平:竖直）")
    ap.add_argument("--phi", type=float, default=40.0, help="接触面摩擦角（°）")
    ap.add_argument("--c", type=float, default=0.0, help="粘聚力 kPa（默认 0，即 FERC 2A 假设）")
    ap.add_argument("--E", type=float, default=0.0, help="排水有效率 0–0.67")
    ap.add_argument("--xd", type=float, default=0.10, help="排水线距坝踵 / 基底宽（默认 0.10）")
    ap.add_argument("--b", type=float, default=4.0, help="坝顶宽 m")
    ap.add_argument("--t", type=float, default=0.0, help="下游面顶部垂直段高 m")
    ap.add_argument("--freeboard", type=float, default=6.1, help="坝顶到正常水面 m（本地块 4,120−4,100 ft = 6.1 m）")
    ap.add_argument("--kh", type=float, default=0.0, help="地震系数（仅初筛）")
    ap.add_argument("--detail", action="store_true")
    a = ap.parse_args()

    if a.H:
        hw = a.H - a.freeboard
        print_detail(a.H, a.m, hw, a.phi, a.c, a.E, a.xd, a.b, a.t, a.kh)
        return

    print(f"# 最缓下游坡 m_min（水平:竖直），使正常运行工况满足目标 FS、静力不开裂、合力在中 1/3")
    print(f"# 上游面垂直，坝顶宽 {a.b} m，超高 {a.freeboard} m，c=0，排水线在 {a.xd}·B 处")
    print()
    hdr = f"{'H(m)':>5} {'hw(m)':>6} {'φ':>4} {'E':>5} | {'FS1.5(FERC 2A)':>15} {'FS2.0(EM2200/FERC低)':>21} {'FS3.0(FERC高)':>14}"
    print(hdr)
    print("-" * len(hdr))
    for H in (20.0, 25.0, 30.0):
        hw = H - a.freeboard
        for phi in (35.0, 40.0, 45.0):
            for E in (0.0, 0.5):
                ms = [min_slope(H, hw, phi, 0.0, E, a.xd, fs, a.b, a.t) for fs in (1.5, 2.0, 3.0)]
                cells = [("≤%.2f" % M_SCAN_MIN) if abs(x - M_SCAN_MIN) < 1e-9 else ("%.2f" % x) for x in ms]
                print(f"{H:5.0f} {hw:6.1f} {phi:4.0f} {E:5.2f} | {cells[0]:>15} {cells[1]:>21} {cells[2]:>14}")
        print()
    print("# 注：'≤0.50' 表示稳定不是控制因素，实际取 EM 2200 的 0.7–0.8 即可；")
    print("#     FS 3.0 一列在 c=0 时通常需要极缓坡——FERC Table 2 的 3.0 是配合试验取得的粘聚力使用的，")
    print("#     无粘聚力应看 Table 2A 的 1.5。两列并排是为了显示'试验投入 vs 断面宽度'的取舍。")
    print("#     下游坡 > 1.0 意味着重力坝已不经济，此时的比较对象是堆石坝，见资料 1–3。")
    print()
    print("# 几何：基底宽 B、混凝土断面积、水平基底时坝踵处的超挖深度 tan24°·B（见 hillside_toe_geometry.py）")
    hdr2 = f"{'H(m)':>5} {'m':>5} {'B(m)':>6} {'面积(m²/m)':>10} {'坝踵超挖(m)':>11}"
    print(hdr2)
    print("-" * len(hdr2))
    for H in (20.0, 25.0, 30.0):
        for m in (0.7, 0.8, 1.0):
            pts, B = section(H, m, a.b, a.t)
            A, _, _ = area_centroid(pts)
            print(f"{H:5.0f} {m:5.2f} {B:6.2f} {A:10.1f} {TAN_BETA*B:11.2f}")


if __name__ == "__main__":
    main()
