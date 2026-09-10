#!/usr/bin/env python3
"""智能压实（IC）数据判据（资料 11：FHWA HIF-24-097、HIF-17-036/037/046、HIF-12-002 附录 D）。

把 FHWA 文件里可编程的部分写成函数：
- ICMV 公式：CMV = 300·A2Ω/AΩ；MDP = Pg − W·v·(sinα + A'/g) − (m·v + b)；MDP80 = 150 − 1.37·MDP。
- 压实曲线取目标 ICMV：遍间增幅 < 5% 的第一遍（通用规范附录 D）。
- 生产判据：评价区 ≥ 25,000 ft²；≥ 90% 面积达最优遍数；≥ 70% 面积达目标 ICMV（FHWA 通用规范）；
  EU CCC：≥ 80% 面积达目标；弱区 = 10 百分位（均值 − 1.28σ）；CoV ≤ 20%（NCHRP 24-45 ≤ 25%）。
- 相关性门槛：R ≥ 0.7（R² ≥ 0.5）——2024 TechBrief 同时说明影响深度不同时相关差不等于 ICMV 错。
- 色图阈值：< 20% 目标 = 红，20%–目标 = 黄，目标–150% = 绿，> 150% = 蓝。
- 经验半变异函数与指数模型（Veta 用指数模型；sill 小、range 长 = 更均匀）。
- GPS 校验：压路机与流动站坐标差 ≤ 40 mm（FHWA 通用规范；各州 6–12 in）。

只做判据与统计，不选设备。用法：
  python3 ic_data_analysis.py --demo
"""
from __future__ import annotations

import argparse
import math
import random

FT2_PER_M2 = 10.7639
MIN_EVAL_AREA_FT2 = 25000.0
CMV_C = 300.0


# ---------------------------------------------------------------------------
# 1. ICMV 公式
# ---------------------------------------------------------------------------
def cmv(a_2omega: float, a_omega: float, c: float = CMV_C) -> float:
    """CMV = C × A2Ω / AΩ（二次谐波 / 基波加速度幅值）。"""
    return c * a_2omega / a_omega


def bouncing_value(a_half_omega: float, a_omega: float, c: float = CMV_C) -> float:
    """RMV/BV = C × A0.5Ω / AΩ；> 0 表示双跳（Dynapac 在 BV≈14 时自动降幅）。"""
    return c * a_half_omega / a_omega


def mdp_kj_s(pg_kw: float, weight_kn: float, v_m_s: float, slope_deg: float,
             accel_m_s2: float, m_kj_m: float, b_kj_s: float) -> float:
    """MDP = Pg − W·v·(sinα + A'/g) − (m·v + b)；标定硬面 MDP = 0，正值更松、负值更密。"""
    g = 9.81
    return pg_kw - weight_kn * v_m_s * (math.sin(math.radians(slope_deg)) + accel_m_s2 / g) - (m_kj_m * v_m_s + b_kj_s)


def mdp80(mdp: float) -> float:
    """MDP80 = 150 − 1.37·MDP（标定面 150；MDP = 108.47 kJ/s → 1）。"""
    return 150.0 - 1.37 * mdp


def mdp40(mdp: float) -> float:
    return 150.0 - 2.75 * mdp


# ---------------------------------------------------------------------------
# 2. 压实曲线 → 最优遍数与目标 ICMV
# ---------------------------------------------------------------------------
def target_from_curve(icmv_by_pass: dict[int, float], rel_change: float = 0.05) -> tuple[int, float]:
    """遍间增幅第一次 < rel_change 的那一遍为最优遍数，其 ICMV 为目标（通用规范附录 D）。
    若 ICMV 先升后降（HMA 常见），取峰值遍。"""
    ks = sorted(icmv_by_pass)
    for a, b in zip(ks, ks[1:]):
        va, vb = icmv_by_pass[a], icmv_by_pass[b]
        if vb <= va:                     # 到顶或开始下降
            return a, va
        if (vb - va) / va < rel_change:
            return b, vb
    return ks[-1], icmv_by_pass[ks[-1]]


# ---------------------------------------------------------------------------
# 3. 生产区判据
# ---------------------------------------------------------------------------
def mean_std(xs: list[float]) -> tuple[float, float]:
    n = len(xs)
    m = sum(xs) / n
    s = math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1)) if n > 1 else 0.0
    return m, s


def coverage_criteria(pass_counts: list[int], icmvs: list[float], target_passes: int, target_icmv: float,
                      area_m2: float, pass_frac_req: float = 0.90, icmv_frac_req: float = 0.70) -> dict:
    """FHWA 通用规范：评价区 ≥ 25,000 ft²；≥ 90% 面积达最优遍数；≥ 70% 面积达目标 ICMV。"""
    n = len(pass_counts)
    pass_ok = sum(1 for p in pass_counts if p >= target_passes) / n
    icmv_ok = sum(1 for v in icmvs if v >= target_icmv) / n
    area_ok = area_m2 * FT2_PER_M2 >= MIN_EVAL_AREA_FT2
    return {"area_ok": area_ok, "pass_frac": pass_ok, "pass_ok": pass_ok >= pass_frac_req,
            "icmv_frac": icmv_ok, "icmv_ok": icmv_ok >= icmv_frac_req,
            "accept_qc": area_ok and pass_ok >= pass_frac_req and icmv_ok >= icmv_frac_req}


def eu_ccc_criteria(icmvs: list[float], target_icmv: float, area_frac_req: float = 0.80,
                    cov_max: float = 0.20) -> dict:
    """CEN/TS 17006 要点：≥ 80% 面积达目标；弱区阈值 = 10 百分位 ≈ 均值 − 1.28σ；CoV ≤ 20%。"""
    m, s = mean_std(icmvs)
    frac = sum(1 for v in icmvs if v >= target_icmv) / len(icmvs)
    weak_thr = m - 1.28 * s
    weak_frac = sum(1 for v in icmvs if v < weak_thr) / len(icmvs)
    cov = s / m if m else float("inf")
    return {"mean": m, "std": s, "cov": cov, "frac_at_target": frac, "weak_threshold": weak_thr,
            "weak_frac": weak_frac, "area_ok": frac >= area_frac_req, "cov_ok": cov <= cov_max}


def issmge_criteria(icmvs: list[float], spec_min: float) -> dict:
    """ISSMGE 2005 / 奥地利 RVS：低于规定最小值的轨迹 ≤ 10%；实测最小 ≥ 80% 规定最小；
    最大 ≤ 150% 规定最小；一遍内 σ/均值 ≤ 20%。"""
    m, s = mean_std(icmvs)
    below = sum(1 for v in icmvs if v < spec_min) / len(icmvs)
    return {"below_min_frac": below, "below_ok": below <= 0.10,
            "min_ok": min(icmvs) >= 0.8 * spec_min, "max_ok": max(icmvs) <= 1.5 * spec_min,
            "std_ok": (s / m) <= 0.20}


def color_class(v: float, target: float) -> str:
    """HIF-17-036 定制四色：< 20% 目标 红（疑似失效）；20%–目标 黄（软）；目标–150% 绿；> 150% 蓝（硬）。"""
    if v < 0.2 * target:
        return "red"
    if v < target:
        return "yellow"
    if v <= 1.5 * target:
        return "green"
    return "blue"


# ---------------------------------------------------------------------------
# 4. 相关性
# ---------------------------------------------------------------------------
def pearson_r(xs: list[float], ys: list[float]) -> float:
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return sxy / math.sqrt(sxx * syy) if sxx and syy else 0.0


def correlation_ok(r: float, r_min: float = 0.7) -> bool:
    return r >= r_min


# ---------------------------------------------------------------------------
# 5. 半变异函数（一维沿碾压方向，指数模型）
# ---------------------------------------------------------------------------
def experimental_semivariogram(values: list[float], spacing_m: float, max_lag_m: float,
                               lag_step_m: float) -> list[tuple[float, float]]:
    """γ(h) = 1/(2N(h)) Σ [Z(u+h) − Z(u)]²，等间距一维序列。"""
    out = []
    lag = lag_step_m
    while lag <= max_lag_m:
        k = int(round(lag / spacing_m))
        pairs = [(values[i + k] - values[i]) ** 2 for i in range(len(values) - k)]
        if pairs:
            out.append((lag, sum(pairs) / (2 * len(pairs))))
        lag += lag_step_m
    return out


def fit_exponential(semivar: list[tuple[float, float]]) -> dict:
    """γ(h) = c·[1 − exp(−3h/a)]，用网格搜索拟合 sill c 与 range a（nugget 取 0）。"""
    hs = [h for h, _ in semivar]
    gs = [g for _, g in semivar]
    best = None
    c_candidates = [max(gs) * f for f in (0.6, 0.8, 1.0, 1.2, 1.5)]
    a_candidates = [max(hs) * f for f in (0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0)]
    for c in c_candidates:
        for a in a_candidates:
            err = sum((c * (1 - math.exp(-3 * h / a)) - g) ** 2 for h, g in semivar)
            if best is None or err < best[0]:
                best = (err, c, a)
    return {"sill": best[1], "range_m": best[2], "sse": best[0]}


# ---------------------------------------------------------------------------
# 6. GPS 校验
# ---------------------------------------------------------------------------
def gps_check(roller_xy: tuple[float, float], rover_xy: tuple[float, float], tol_m: float = 0.040) -> dict:
    d = math.hypot(roller_xy[0] - rover_xy[0], roller_xy[1] - rover_xy[1])
    return {"diff_m": d, "ok": d <= tol_m}


# ---------------------------------------------------------------------------
def demo() -> None:
    rng = random.Random(11)
    print("== 1. ICMV 公式 ==")
    print(f"  CMV(A2Ω=0.6 g, AΩ=6 g) = {cmv(0.6, 6.0):.1f}；BV(A0.5Ω=0.2 g) = {bouncing_value(0.2, 6.0):.1f}（>0 → 双跳）")
    m = mdp_kj_s(pg_kw=60, weight_kn=120, v_m_s=1.0, slope_deg=2, accel_m_s2=0.0, m_kj_m=10, b_kj_s=15)
    print(f"  MDP 示意 = {m:.1f} kJ/s → MDP80 = {mdp80(m):.0f}（标定面 150；MDP 108.47 → 1）；校核 mdp80(108.47) = {mdp80(108.47):.1f}")

    print("\n== 2. 压实曲线 → 最优遍数与目标 ICMV（通用规范附录 D：遍间增幅 < 5%）==")
    curve = {1: 18, 2: 24, 3: 28, 4: 30.5, 5: 31.8, 6: 32.3, 7: 32.5, 8: 32.4}
    p, t = target_from_curve(curve)
    print(f"  曲线 {curve} → 最优遍数 {p}，目标 ICMV {t}")

    print("\n== 3. 生产区判据（合成 2,000 格；目标 ICMV 32、最优 5 遍）==")
    n = 2000
    passes = [rng.choice([4, 5, 5, 5, 6, 6, 7]) for _ in range(n)]
    icmvs = [max(1.0, rng.gauss(36, 7)) for _ in range(n)]
    r = coverage_criteria(passes, icmvs, 5, 32.0, area_m2=3000)
    print(f"  FHWA：面积达标 {r['area_ok']}，≥5 遍 {r['pass_frac']*100:.1f}%（要 90），≥目标 ICMV {r['icmv_frac']*100:.1f}%（要 70）→ QC {'通过' if r['accept_qc'] else '不通过'}")
    e = eu_ccc_criteria(icmvs, 32.0)
    print(f"  EU CCC：均值 {e['mean']:.1f}，σ {e['std']:.1f}，CoV {e['cov']*100:.1f}%（≤20），达目标 {e['frac_at_target']*100:.1f}%（≥80），"
          f"弱区阈 {e['weak_threshold']:.1f}，弱区占 {e['weak_frac']*100:.1f}%")
    i = issmge_criteria(icmvs, spec_min=28.0)
    print(f"  ISSMGE（规定最小 28）：低于最小 {i['below_min_frac']*100:.1f}%（≤10）→ {i['below_ok']}；最小≥80% {i['min_ok']}；最大≤150% {i['max_ok']}；σ≤20% {i['std_ok']}")
    cls = {}
    for v in icmvs:
        cls[color_class(v, 32.0)] = cls.get(color_class(v, 32.0), 0) + 1
    print(f"  色图分档（红/黄/绿/蓝）：{ {k: round(v/n*100,1) for k, v in cls.items()} } %")

    print("\n== 4. 相关性（合成：均匀地基 vs 非均匀地基）==")
    x = [rng.uniform(20, 60) for _ in range(30)]
    y_uniform = [1.5 * v + rng.gauss(0, 6) for v in x]
    y_hetero = [1.5 * v + rng.gauss(0, 35) for v in x]
    for name, y in (("均匀地基", y_uniform), ("非均匀地基", y_hetero)):
        rr = pearson_r(x, y)
        print(f"  {name}：R = {rr:.2f} → {'≥0.7 可用' if correlation_ok(rr) else '<0.7（2024 TechBrief：相关差不等于 ICMV 错，是影响深度不同）'}")

    print("\n== 5. 半变异函数（沿碾压方向 200 m、每 0.5 m 一值）==")
    for name, amp, wavelength in (("均匀（只有 σ=2 的噪声）", 0.0, 1.0), ("斑块状（叠加 ±8、波长 80 m 的软硬带）", 8.0, 80.0)):
        vals = [30.0 + amp * math.sin(2 * math.pi * (k * 0.5) / wavelength) + rng.gauss(0, 2.0) for k in range(400)]
        sv = experimental_semivariogram(vals, 0.5, 60, 2.5)
        fit = fit_exponential(sv)
        print(f"  {name}：sill {fit['sill']:.1f}，range {fit['range_m']:.0f} m（sill 越大越不均匀；range 是软硬带的空间尺度）")

    print("\n== 6. GPS 校验（FHWA 40 mm）==")
    for d in ((0.0, 0.03), (0.05, 0.02)):
        g = gps_check((100.0, 200.0), (100.0 + d[0], 200.0 + d[1]))
        print(f"  差 {g['diff_m']*1000:.0f} mm → {'合格' if g['ok'] else '不合格，重校'}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--demo", action="store_true")
    ap.parse_args()
    demo()


if __name__ == "__main__":
    main()
