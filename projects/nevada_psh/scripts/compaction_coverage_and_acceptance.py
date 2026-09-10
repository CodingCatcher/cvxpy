#!/usr/bin/env python3
"""压实回数网格管理、摊铺厚换算、竣工面（出来形）面管理判定与 DEM 差分方量不确定度（资料 10）。

依据：
- 国土交通省《TS・GNSS を用いた盛土の締固め管理要領》（2020）：施工范围切成正方形管理块
  （推土机 0.25 m、压路机 0.50 m），碾轮/履带每经过一块计 1 回；只取 FIX 解；GNSS 精度
  水平 ±20 mm、垂直 ±30 mm；本施工まき出し厚 = 所定仕上り厚 × (试验まき出し厚/试验仕上り厚)。
- 《3 次元計測技術を用いた出来形管理要領（案）》（2025）与《施工管理基準及び規格値（案）》：
  面管理 = 全面每点与设计面的标高较差，平均值与个别值各有规格值，99.7% 的点满足个别值
  即视为"全部满足"（弃却 ≤0.3%）；计测密度 ≥1 点/m²；个别值规格内含 ±50 mm 计测精度。
- 无人机/机载 DEM 差分方量：高程误差 σz 在面积 A 上的方量不确定度 ≈ A·σz（系统性）或
  A·σz/√N（随机、N 个独立像元）——两者都给，供交付物 4/6 标注精度。

只做规则与几何，不选系统。用法：
  python3 compaction_coverage_and_acceptance.py --demo
"""
from __future__ import annotations

import argparse
import math
import random

BLOCK_SIZE = {"bulldozer": 0.25, "tire_roller": 0.50, "vibratory_roller": 0.50}
GNSS_TOL = {"xy_mm": 20.0, "z_mm": 30.0}
INDIVIDUAL_INCLUDES_MEAS_MM = 50.0

# 面管理规格值（mm）：(平均值, 个别值)。负号表示只允许偏低（河川盛土）；道路为 ±。
SURFACE_TOL = {
    "river_fill_crest":        (-50, -150),
    "river_fill_slope_steep":  (-50, -170),   # 坡比陡于 1:4
    "river_fill_slope_flat":   (-60, -170),   # 坡比缓于 1:4（含小段）
    "road_fill_crest":         (50, 150),     # ±
    "road_fill_slope":         (80, 190),     # ±
}


# ---------------------------------------------------------------------------
# 1. 压实回数网格：直线平行碾压的覆盖率与"踏み残し"
# ---------------------------------------------------------------------------
def pass_count_grid(width_m: float, length_m: float, block: float,
                    lane_centers: list[float], drum_width: float, passes_per_lane: int,
                    xy_noise_m: float = 0.0, seed: int = 0) -> list[list[int]]:
    """沿 y 方向直线碾压：每条车道中心 x=c，滚轮宽 drum_width，来回 passes_per_lane 次。
    每次经过时给车道中心加一次横向 GNSS 噪声（正态，σ=xy_noise_m），模拟定位误差造成的错计。
    返回 nx × ny 的回数网格。"""
    rng = random.Random(seed)
    nx, ny = int(math.ceil(width_m / block)), int(math.ceil(length_m / block))
    grid = [[0] * ny for _ in range(nx)]
    for c in lane_centers:
        for _ in range(passes_per_lane):
            cc = c + (rng.gauss(0.0, xy_noise_m) if xy_noise_m > 0 else 0.0)
            x0, x1 = cc - drum_width / 2, cc + drum_width / 2
            i0 = max(0, int(math.floor(x0 / block)))
            i1 = min(nx - 1, int(math.floor((x1 - 1e-9) / block)))
            for i in range(i0, i1 + 1):
                # 块中心在滚轮范围内才计（要领：碾轮通过该块即计 1 回；这里取块中心判据）
                xc = (i + 0.5) * block
                if x0 <= xc <= x1:
                    for j in range(ny):
                        grid[i][j] += 1
    return grid


def coverage_stats(grid: list[list[int]], required: int, over_limit: int | None = None) -> dict:
    n = sum(len(r) for r in grid)
    ok = sum(1 for r in grid for v in r if v >= required)
    under = n - ok
    over = sum(1 for r in grid for v in r if over_limit is not None and v >= over_limit)
    return {"blocks": n, "ok_frac": ok / n, "under_frac": under / n, "over_frac": over / n if over_limit else 0.0}


def lane_centers(width_m: float, drum_width: float, overlap_m: float) -> list[float]:
    """按滚轮宽与搭接宽排车道，从一侧排到另一侧（最后一道贴边）。"""
    step = drum_width - overlap_m
    cs, c = [], drum_width / 2
    while c + drum_width / 2 < width_m - 1e-9:
        cs.append(c)
        c += step
    cs.append(width_m - drum_width / 2)
    return cs


# ---------------------------------------------------------------------------
# 2. 摊铺厚换算与试验施工回数（沉降收敛）
# ---------------------------------------------------------------------------
def spread_thickness(finished_thickness_m: float, trial_spread_m: float, trial_finished_m: float) -> float:
    """本施工まき出し厚 = 所定仕上り厚 × (试验まき出し厚 / 试验仕上り厚)。"""
    return finished_thickness_m * trial_spread_m / trial_finished_m


def passes_from_settlement(settlement_by_pass: dict[int, float], rel_tol: float = 0.05) -> int:
    """岩块料（不能用密度）：取表面沉降随回数收敛的拐点——相邻回数增量小于总沉降 rel_tol 的第一处。"""
    ks = sorted(settlement_by_pass)
    total = settlement_by_pass[ks[-1]]
    for a, b in zip(ks, ks[1:]):
        if settlement_by_pass[b] - settlement_by_pass[a] <= rel_tol * total:
            return b
    return ks[-1]


# ---------------------------------------------------------------------------
# 3. 出来形面管理判定
# ---------------------------------------------------------------------------
def surface_acceptance(dz_mm: list[float], mean_tol_mm: float, indiv_tol_mm: float,
                       one_sided: bool | None = None, reject_frac_max: float = 0.003) -> dict:
    """dz = 实测高程 − 设计高程（mm）。one_sided=True 表示只允许负偏（河川盛土规格值为负号）。
    返回平均、最大、最小、超出个别值的点数与比例、判定。"""
    if one_sided is None:
        one_sided = mean_tol_mm < 0
    n = len(dz_mm)
    mean = sum(dz_mm) / n
    mx, mn = max(dz_mm), min(dz_mm)
    if one_sided:
        # 河川盛土规格值只写负号：只限制偏低（平均 ≥ −50、个别 ≥ −150），偏高不限
        mean_ok = mean >= mean_tol_mm
        bad = [d for d in dz_mm if d < indiv_tol_mm]
    else:
        mean_ok = abs(mean) <= abs(mean_tol_mm)
        bad = [d for d in dz_mm if abs(d) > abs(indiv_tol_mm)]
    reject_frac = len(bad) / n
    return {"n": n, "mean": mean, "max": mx, "min": mn, "reject_n": len(bad), "reject_frac": reject_frac,
            "mean_ok": mean_ok, "indiv_ok": reject_frac <= reject_frac_max,
            "accept": mean_ok and reject_frac <= reject_frac_max}


def required_points(plan_area_m2: float, density_per_m2: float = 1.0) -> int:
    return int(math.ceil(plan_area_m2 * density_per_m2))


# ---------------------------------------------------------------------------
# 4. DEM 差分方量的不确定度
# ---------------------------------------------------------------------------
def volume_uncertainty(area_m2: float, sigma_z_m: float, cell_m: float | None = None) -> dict:
    """系统性：A·σz；随机（像元独立）：A·σz/√N，N = A/cell²。真实误差介于两者之间（空间相关）。"""
    systematic = area_m2 * sigma_z_m
    out = {"systematic_m3": systematic}
    if cell_m:
        n = area_m2 / cell_m ** 2
        out["random_m3"] = systematic / math.sqrt(n)
    return out


# ---------------------------------------------------------------------------
def demo() -> None:
    print("== 1. 压实回数网格：20 m × 60 m 填筑面，2.1 m 宽振动碾，规定 8 回（每道来回 4 趟 × 2 = 8）==")
    W, L, drum = 20.0, 60.0, 2.1
    for overlap in (0.0, 0.2, 0.4):
        cs = lane_centers(W, drum, overlap)
        for noise in (0.0, 0.02, 0.10):
            g = pass_count_grid(W, L, BLOCK_SIZE["vibratory_roller"], cs, drum, 8, xy_noise_m=noise, seed=1)
            s = coverage_stats(g, required=8, over_limit=12)
            print(f"  搭接 {overlap:3.1f} m，{len(cs):2d} 道，横向噪声 σ={noise*1000:3.0f} mm → "
                  f"达标块 {s['ok_frac']*100:5.1f}%，欠压 {s['under_frac']*100:4.1f}%，≥过压上限 {s['over_frac']*100:4.1f}%")
    print("  注：要领要求 GNSS 水平 ±20 mm；σ=100 mm 是 FLOAT 解混入的量级，说明为什么只取 FIX 解。")

    print("\n== 2. まき出し厚：试验 35 cm 摊铺 → 压后 30 cm；要 30 cm 仕上り厚 → 本施工摊铺", end=" ")
    print(f"{spread_thickness(0.30, 0.35, 0.30)*100:.1f} cm；岩块料要 60 cm 层 → 摊铺 {spread_thickness(0.60, 0.35, 0.30)*100:.1f} cm ==")
    settle = {0: 0, 2: 40, 4: 62, 6: 72, 8: 76, 10: 78, 12: 79}
    print(f"   沉降收敛（mm，随回数）{settle} → 拐点回数 {passes_from_settlement(settle)}（增量 ≤5% 总沉降）")

    print("\n== 3. 出来形面管理判定（合成数据）==")
    rng = random.Random(7)
    area = 5000.0
    n = required_points(area)
    for label, bias, sd, key in (("道路路体 天端，偏差 −10 mm、σ 40 mm", -10, 40, "road_fill_crest"),
                                 ("道路路体 法面，偏差 +30 mm、σ 70 mm", 30, 70, "road_fill_slope"),
                                 ("河川盛土 天端，偏差 −20 mm、σ 40 mm", -20, 40, "river_fill_crest"),
                                 ("河川盛土 天端，偏差 −60 mm、σ 40 mm（偏低）", -60, 40, "river_fill_crest")):
        dz = [rng.gauss(bias, sd) for _ in range(n)]
        mt, it = SURFACE_TOL[key]
        r = surface_acceptance(dz, mt, it)
        print(f"  {label}：n={r['n']}，平均 {r['mean']:6.1f}，最大 {r['max']:6.1f}，最小 {r['min']:7.1f}，"
              f"弃却 {r['reject_n']} ({r['reject_frac']*100:.2f}%) → 平均{'合' if r['mean_ok'] else '不合'}/"
              f"个别{'合' if r['indiv_ok'] else '不合'} → {'合格' if r['accept'] else '不合格'}")
    print("  注：河川盛土规格值只有负号（只许偏低不许偏高），道路为 ±；个别值规格内含 ±50 mm 计测精度。")

    print("\n== 4. DEM 差分方量不确定度 ==")
    for A, sz, cell in ((10_000, 0.05, 0.05), (100_000, 0.05, 0.05), (100_000, 0.03, 0.05)):
        u = volume_uncertainty(A, sz, cell)
        print(f"  面积 {A:>7,} m²，σz {sz*100:.0f} cm，像元 {cell*100:.0f} cm → 系统性 ±{u['systematic_m3']:,.0f} m³，"
              f"纯随机 ±{u['random_m3']:.1f} m³（真实在两者之间；料堆经验 ±2–3%）")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    demo()


if __name__ == "__main__":
    main()
