#!/usr/bin/env python3
"""土石方施工产量与阻力计算（资料 9：Caterpillar Performance Handbook 第 49 版，2019）。

把手册"矿山与土石方基础"(第 28 节)、"表格"(第 30 节)、"履带推土机/松土器"(第 19 节)、
"压实机"(第 23 节)、"液压挖掘机"(第 7 节)、"铰接卡车/矿用卡车"(第 1、10 节) 里
可编程的公式写成函数，并用手册自带算例复算（--demo）。所有单位为公制；
手册的英制常数在注释里给出，方便对照。

只做"公式 + 手册数字"，不做选型。曲线（推土产量、牵引力-速度）只在图上有，
这里用从图上读的少数点做插值/拟合，误差 ±10–15%，只供量级判断。

用法：
  python3 earthmoving_production.py --demo     # 复算手册算例（D8T 推土、D10T2 松土、825G 压实、320 挖掘、631G 铲运链）
  python3 earthmoving_production.py --site     # 用本地块的坡度/岩性/海拔算几组"情况"数字
"""
from __future__ import annotations

import argparse
import math

# ---------------------------------------------------------------------------
# 常数
# ---------------------------------------------------------------------------
YD3_PER_M3 = 1.30795
LB_PER_KG = 2.20462
FT_PER_M = 3.28084
HP_PER_KW = 1.34102
KG_PER_T_PER_PCT = 10.0        # 手册：每 1% 坡度 = 10 kg/t（英制 20 lb/ton）
RR_BASE_PCT = 2.0              # 硬路面基础滚动阻力 2% GMW
RR_PER_CM_PCT = 0.6            # 轮胎每下陷 1 cm 加 0.6%（英制每 1 in 加 1.5%）
GRADE_HP_CONST_METRIC = 273.75 # HP = GMW(kg) × TR × km/h / 273.75（英制 375）
COMPACT_CONST_ENGLISH = 16.3   # CCY/h = W(ft) × S(mph) × L(in) × 16.3 / P

# ---------------------------------------------------------------------------
# 第 30 节 材料表（松方 kg/m³, 原岩方 kg/m³, 荷载系数 LF = 松/原）
# ---------------------------------------------------------------------------
MATERIALS = {
    "basalt":                 (1960, 2970, 0.67),
    "caliche":                (1250, 2260, 0.55),
    "clay_natural_bed":       (1660, 2020, 0.82),
    "clay_dry":               (1480, 1840, 0.81),
    "clay_wet":               (1660, 2080, 0.80),
    "clay_gravel_dry":        (1420, 1660, 0.85),
    "clay_gravel_wet":        (1540, 1840, 0.85),
    "decomposed_rock_75_25":  (1960, 2790, 0.70),
    "decomposed_rock_50_50":  (1720, 2280, 0.75),
    "decomposed_rock_25_75":  (1570, 1960, 0.80),
    "earth_dry_packed":       (1510, 1900, 0.80),
    "earth_wet_excavated":    (1600, 2020, 0.79),
    "earth_loam":             (1250, 1540, 0.81),
    "granite_broken":         (1660, 2730, 0.61),
    "gravel_pitrun":          (1930, 2170, 0.89),
    "gravel_dry":             (1510, 1690, 0.89),
    "gravel_dry_6_50mm":      (1690, 1900, 0.89),
    "gravel_wet_6_50mm":      (2020, 2260, 0.89),
    "limestone_broken":       (1540, 2610, 0.59),
    "sand_dry_loose":         (1420, 1600, 0.89),
    "sand_damp":              (1690, 1900, 0.89),
    "sand_wet":               (1840, 2080, 0.89),
    "sand_clay_loose":        (1600, 2020, 0.79),
    "sand_gravel_dry":        (1720, 1930, 0.89),
    "sand_gravel_wet":        (2020, 2230, 0.91),
    "sandstone":              (1510, 2520, 0.60),
    "shale":                  (1250, 1660, 0.75),
    "slag_broken":            (1750, 2940, 0.60),
    "stone_crushed":          (1600, 2670, 0.60),
    "top_soil":               (950, 1370, 0.70),
    "traprock_broken":        (1750, 2610, 0.67),
}

# 第 30 节 典型滚动阻力（% GMW）
ROLLING_RESISTANCE_TABLE = {
    "hard_smooth_stabilized_no_penetration": 1.5,
    "firm_smooth_maintained_flexes_slightly": 3.0,
    "snow_packed": 2.5,
    "snow_loose": 4.5,
    "dirt_rutted_25mm": 4.0,
    "dirt_rutted_50mm": 5.0,
    "dirt_rutted_100mm": 8.0,
    "loose_sand_or_gravel": 10.0,
    "rutted_soft_200mm": 14.0,
    "mud_300mm": 20.0,
}

# 第 30 节 附着系数（轮胎, 履带）
TRACTION = {
    "concrete": (0.90, 0.45),
    "clay_loam_dry": (0.55, 0.90),
    "clay_loam_wet": (0.45, 0.70),
    "rutted_clay_loam": (0.40, 0.70),
    "dry_sand": (0.20, 0.30),
    "wet_sand": (0.40, 0.50),
    "quarry_pit": (0.65, 0.55),
    "gravel_road_loose": (0.36, 0.50),
    "packed_snow": (0.20, 0.27),
    "ice": (0.12, 0.12),
    "firm_earth": (0.55, 0.90),
    "loose_earth": (0.45, 0.60),
    "coal_stockpiled": (0.45, 0.60),
}

# ---------------------------------------------------------------------------
# 方量换算（第 28 节）
# ---------------------------------------------------------------------------
def load_factor_from_swell(swell_pct: float) -> float:
    """LF = 100 / (100 + 膨胀%)。膨胀 25% → LF 0.80。"""
    return 100.0 / (100.0 + swell_pct)


def swell_from_load_factor(lf: float) -> float:
    return 100.0 / lf - 100.0


def bank_to_loose(bcm: float, lf: float) -> float:
    return bcm / lf


def loose_to_bank(lcm: float, lf: float) -> float:
    return lcm * lf


def bank_to_compacted(bcm: float, shrinkage_factor: float) -> float:
    """SF = 压实方 / 原方（手册 631G 算例取 0.85；堆石可能 >1，见资料 7/8）。"""
    return bcm * shrinkage_factor


# ---------------------------------------------------------------------------
# 阻力、有效坡度、牵引（第 28 节）
# ---------------------------------------------------------------------------
def rolling_resistance_pct(penetration_cm: float = 0.0, base_pct: float = RR_BASE_PCT) -> float:
    """RR% = 2% + 0.6%/cm 下陷（手册经验式；轮胎；履带机另有表）。"""
    return base_pct + RR_PER_CM_PCT * penetration_cm


def rr_factor_kg_per_t(penetration_cm: float = 0.0) -> float:
    """RR 系数 kg/t = 20 + 6/cm（= 10 kg/t × RR%）。"""
    return KG_PER_T_PER_PCT * rolling_resistance_pct(penetration_cm)


def effective_grade_pct(rr_pct: float, grade_pct: float) -> float:
    """总阻力（有效坡度）% = 滚动阻力% + 坡度%（下坡取负；下坡时又叫 grade assistance）。"""
    return rr_pct + grade_pct


def total_resistance_kg(gmw_t: float, rr_pct: float, grade_pct: float) -> float:
    """所需牵引力 kg = GMW(t) × 10 kg/t × 有效坡度%。"""
    return gmw_t * KG_PER_T_PER_PCT * effective_grade_pct(rr_pct, grade_pct)


def usable_rimpull_kg(weight_on_drivers_kg: float, traction_coeff: float) -> float:
    """可用牵引力受附着限制：驱动轮/履带上的重量 × 附着系数。"""
    return weight_on_drivers_kg * traction_coeff


def grade_power_kw(gmw_kg: float, total_resistance_frac: float, speed_kmh: float) -> float:
    """坡道功率（第 10 节）：HP = GMW(kg) × TR × km/h / 273.75；返回 kW。
    手册算例：317 520 kg × 0.10 × 13.2 km/h → 1530 HP。"""
    hp = gmw_kg * total_resistance_frac * speed_kmh / GRADE_HP_CONST_METRIC
    return hp / HP_PER_KW


def slope_deg_to_pct(deg: float) -> float:
    return math.tan(math.radians(deg)) * 100.0


def slope_pct_to_deg(pct: float) -> float:
    return math.degrees(math.atan(pct / 100.0))


# ---------------------------------------------------------------------------
# 牵引力-速度曲线（第 1 节铰接卡车，从图上读点：速度 km/h → 牵引力 t）
# 读图误差 ±10–15%。用法：所需牵引力 = GMW(t) × 总阻力%/100，再反查速度。
# ---------------------------------------------------------------------------
RIMPULL_CURVES = {
    # 740C EJ, 29.5R25, 海平面；空车 36.0 t，满载 74.0 t，最高约 55 km/h
    "740C_EJ": {"empty_t": 36.0, "loaded_t": 74.0,
                "pts": [(0, 45.0), (3, 30.0), (5, 22.0), (7, 20.0), (10, 14.5), (14, 10.5),
                        (20, 7.2), (30, 5.0), (40, 3.9), (50, 3.2), (55, 2.5)]},
    # 730C2, 23.5R25；空车 23.7 t，满载 51.7 t，最高约 55 km/h
    "730C2": {"empty_t": 23.7, "loaded_t": 51.7,
              "pts": [(0, 35.0), (3, 22.0), (5, 17.0), (8, 11.5), (10, 10.0), (15, 7.5),
                      (20, 6.0), (25, 5.0), (30, 4.2), (40, 3.2), (50, 2.5), (55, 2.0)]},
}

# 740C EJ 制动/缓速曲线读数（第 1 节 p.1-26）：有利有效坡度% → 最高安全下坡速度 km/h
RETARDER_740C_EJ = {
    "loaded_74t": {5: 55, 10: 40, 15: 28, 20: 22, 25: 17, 30: 13, 35: 10},
    "empty_36t":  {10: 55, 15: 48, 20: 38, 25: 30, 30: 24, 35: 19},
}


def retarder_speed_kmh(table: dict, favorable_grade_pct: float, top_speed: float = 55.0) -> float:
    """缓速曲线：有利有效坡度% → 最高安全下坡速度（表内线性插值；低于最小键取最高车速）。"""
    keys = sorted(table)
    if favorable_grade_pct <= keys[0]:
        return min(top_speed, table[keys[0]]) if favorable_grade_pct >= keys[0] else top_speed
    if favorable_grade_pct >= keys[-1]:
        return float(table[keys[-1]])
    for k1, k2 in zip(keys, keys[1:]):
        if k1 <= favorable_grade_pct <= k2:
            v1, v2 = table[k1], table[k2]
            return v1 + (v2 - v1) * (favorable_grade_pct - k1) / (k2 - k1)
    return top_speed


def rimpull_speed_kmh(truck: str, gmw_t: float, total_resistance_pct: float) -> float:
    """按牵引力曲线查某总阻力下能维持的速度（线性插值；超出曲线端点截断）。"""
    c = RIMPULL_CURVES[truck]
    need_t = gmw_t * total_resistance_pct / 100.0
    pts = c["pts"]
    if need_t >= pts[0][1]:
        return 0.0
    if need_t <= pts[-1][1]:
        return float(pts[-1][0])
    for (s1, r1), (s2, r2) in zip(pts, pts[1:]):
        if r2 <= need_t <= r1:
            return s1 + (s2 - s1) * (r1 - need_t) / (r1 - r2)
    return float(pts[-1][0])


# ---------------------------------------------------------------------------
# 推土机产量（第 19 节）：产量 = 图上最大值 × 各修正系数
# ---------------------------------------------------------------------------
DOZER_CORRECTIONS = {
    "operator_excellent": 1.00, "operator_average": 0.75, "operator_poor": 0.60,
    "loose_stockpile": 1.20,
    "hard_to_cut_with_tilt": 0.80, "hard_to_cut_without_tilt": 0.70,
    "hard_to_drift": 0.80,
    "rock_ripped_or_blasted_low": 0.60, "rock_ripped_or_blasted_high": 0.80,
    "slot_dozing": 1.20,
    "side_by_side_low": 1.15, "side_by_side_high": 1.25,
    "visibility_poor": 0.80,
    "efficiency_50min": 0.83, "efficiency_40min": 0.67,
}
BASE_LOOSE_DENSITY = 1370.0   # kg/Lm³（2300 lb/LCY），图的基准密度

# 推土产量曲线读点：(15 m 处 LCY/h, 91 m 处 LCY/h)。U 刀与 SU 刀。图读数 ±10%。
DOZER_CURVE_POINTS_LCY = {
    "D11CD_U": (6200, 1500), "D11_U": (5600, 1300), "D10T2_U": (3000, 800),
    "D9T_U": (2000, 600), "D8T_U": (1400, 400), "D7E_U": (1200, 350),
    "D11_SU": (4400, 1100), "D10T2_SU": (2500, 700), "D9T_SU": (1700, 500),
    "D8T_SU": (1000, 350), "D7E_SU": (900, 300), "D8R_SU": (1050, 300),
    "D7R_SU": (850, 250), "D6T_SU": (720, 180), "D6N_SU": (520, 130),
}


def dozer_grade_factor(grade_pct: float) -> float:
    """% 坡度 vs 推土系数图：下坡为负。图上近似直线：0%→1.0，每 10% 变 0.2；
    手册算例 −15%（下坡）取 1.30。"""
    return 1.0 - 0.02 * grade_pct


def dozer_curve_lcm_h(model_blade: str, distance_m: float) -> float:
    """用两点 1/P = a + b·d（循环时间 = 固定时间 + 行程/速度）拟合曲线，返回 Lm³/h。"""
    p1, p2 = DOZER_CURVE_POINTS_LCY[model_blade]
    d1, d2 = 15.24, 91.44
    b = (1.0 / p2 - 1.0 / p1) / (d2 - d1)
    a = 1.0 / p1 - b * d1
    lcy_h = 1.0 / (a + b * distance_m)
    return lcy_h / YD3_PER_M3


def dozer_production(max_uncorrected: float, factors: list[float],
                     loose_density_kg_m3: float | None = None) -> float:
    """产量 = 未修正最大产量 × Π(系数) × (1370 / 实际松方密度)。单位随输入。"""
    p = max_uncorrected
    for f in factors:
        p *= f
    if loose_density_kg_m3:
        p *= BASE_LOOSE_DENSITY / loose_density_kg_m3
    return p


# ---------------------------------------------------------------------------
# 松土产量（第 19 节，测距法；结果通常比实际高 10–20%）
# ---------------------------------------------------------------------------
def ripper_production_bcm_h(pass_len_m: float, spacing_m: float, depth_m: float,
                            speed_kmh: float, turn_min: float = 0.25,
                            work_min_per_h: float = 60.0) -> dict:
    v_m_min = speed_kmh * 1000.0 / 60.0
    t_pass = pass_len_m / v_m_min + turn_min
    passes_h = work_min_per_h / t_pass
    vol_pass = pass_len_m * spacing_m * depth_m
    gross = vol_pass * passes_h
    return {"min_per_pass": t_pass, "passes_per_h": passes_h, "bcm_per_pass": vol_pass,
            "bcm_h_method": gross, "bcm_h_80pct": 0.8 * gross, "bcm_h_90pct": 0.9 * gross}


def loosening_cost_per_bcm(oo_cost_per_h: float, bcm_h: float) -> float:
    return oo_cost_per_h / bcm_h


# ---------------------------------------------------------------------------
# 挖掘机产量（第 7 节）
# ---------------------------------------------------------------------------
def excavator_production_lcm_h(cycle_s: float, bucket_m3: float, fill_factor: float = 1.0,
                               efficiency: float = 1.0) -> float:
    """产量 = 3600/循环秒 × 斗容 × 充满系数 × 效率。"""
    return 3600.0 / cycle_s * bucket_m3 * fill_factor * efficiency


def excavator_effective_cycles_h(cycle_min: float, avail: float = 0.90, util: float = 0.95,
                                 eff: float = 0.83) -> float:
    """大型挖掘机选型第 3 步：有效循环/h = 60/CT × 可用率 × 利用率 × 作业效率。"""
    return 60.0 / cycle_min * avail * util * eff


# ---------------------------------------------------------------------------
# 运输循环（第 10 节固定时间；第 1/10 节速度）
# ---------------------------------------------------------------------------
def truck_cycle_min(load_min: float, haul_segments: list[tuple[float, float]],
                    return_segments: list[tuple[float, float]],
                    exchange_min: float = 0.7, dump_min: float = 1.1) -> dict:
    """segments: [(距离 m, 速度 km/h), ...]。固定时间：装车交换 0.6–0.8，卸车 1.0–1.2 min。"""
    def travel(segs):
        return sum(d / (v * 1000.0 / 60.0) for d, v in segs)
    haul = travel(haul_segments)
    ret = travel(return_segments)
    total = load_min + exchange_min + dump_min + haul + ret
    return {"load": load_min, "haul": haul, "return": ret, "fixed": exchange_min + dump_min,
            "total": total}


def trucks_to_match_loader(truck_cycle_min_: float, load_min: float) -> float:
    """理论所需车数 = 卡车循环 / 装车时间（不含等待）。"""
    return truck_cycle_min_ / load_min


def pusher_cycle_min(load_min: float) -> float:
    """推土助铲循环 = 1.4 × 装载时间 + 0.25 min（第 28 节）。"""
    return 1.4 * load_min + 0.25


# ---------------------------------------------------------------------------
# 压实产量（第 23/28 节）
# ---------------------------------------------------------------------------
def compactor_production_ccm_h(width_m: float, speed_kmh: float, lift_mm: float, passes: int) -> float:
    """压实方 m³/h = W(m) × S(km/h) × L(mm) / P。英制：CCY/h = W(ft) × S(mph) × L(in) × 16.3 / P。"""
    return width_m * speed_kmh * lift_mm / passes


def compactor_production_ccy_h_english(width_ft: float, speed_mph: float, lift_in: float, passes: int) -> float:
    return width_ft * speed_mph * lift_in * COMPACT_CONST_ENGLISH / passes


# ---------------------------------------------------------------------------
# 拥有成本（第 25 节）
# ---------------------------------------------------------------------------
def hourly_owning_cost(price: float, salvage: float, years: float, hours_per_year: float,
                       rate_interest_insurance_tax: float) -> dict:
    """折旧 = (P − S)/总小时；利息+保险+税 = 平均投资 × 年率 / 年小时，
    平均投资 = [P(N+1) + S(N−1)] / (2N)。"""
    total_h = years * hours_per_year
    depreciation = (price - salvage) / total_h
    avg_inv = (price * (years + 1) + salvage * (years - 1)) / (2 * years)
    iit = avg_inv * rate_interest_insurance_tax / hours_per_year
    return {"depreciation": depreciation, "interest_ins_tax": iit, "owning": depreciation + iit}


# ---------------------------------------------------------------------------
# 手册算例复算
# ---------------------------------------------------------------------------
def check(label: str, got: float, want: float, tol: float = 0.01) -> None:
    ok = abs(got - want) <= tol * max(abs(want), 1e-9)
    print(f"  [{'OK' if ok else 'DIFF'}] {label}: 得 {got:,.2f}，手册 {want:,.2f}")


def demo() -> None:
    print("== 1. 第 28 节：方量换算与 631G 铲运机算例链 ==")
    check("膨胀 25% → LF", load_factor_from_swell(25), 0.80)
    check("20,000 BCY × 1/0.8 = 25,000 LCY", bank_to_loose(20000, 0.8), 25000)
    check("推土助铲循环 1.4×0.7+0.25", pusher_cycle_min(0.7), 1.23)
    check("每台助铲机可带铲运机数 6.95/1.23", 6.95 / 1.23, 5.65, 0.002)
    cycles = 60 / 6.95
    unit = 31 * 0.80 * cycles
    check("单机产量 31 LCY×0.8×8.6", unit, 213, 0.01)
    fleet = unit * 0.83 * 11
    check("车队 0.83×213×11", fleet, 1947, 0.01)
    check("压实需求 0.85×1947 CCY/h", fleet * 0.85, 1655, 0.01)
    check("825G 压实能力 7.4 ft×6 mph×7 in×16.3/3", compactor_production_ccy_h_english(7.4, 6, 7, 3), 1688, 0.001)
    metric = compactor_production_ccm_h(7.4 / FT_PER_M, 6 * 1.609344, 7 * 25.4, 3)
    check("同题公制 m³/h → 换成 CCY/h", metric * YD3_PER_M3, 1688, 0.005)

    print("== 2. 第 28 节：阻力与坡道功率 ==")
    check("RR 系数：下陷 5 cm → 20+6×5 kg/t", rr_factor_kg_per_t(5), 50)
    check("有利坡 20%、RR 50 kg/t → 总有效坡 15%", effective_grade_pct(5, -20), -15)
    check("317 520 kg × 0.10 × 13.2 km/h → 1530 HP", grade_power_kw(317520, 0.10, 13.2) * HP_PER_KW, 1530, 0.002)
    check("631G 满载 88.4 t、RR 200 lb/ton = 10%：TR 17,686 lb",
          total_resistance_kg(88.4 * 0.907185, 10, 0) * LB_PER_KG, 17686, 0.01)

    print("== 3. 第 19 节：D8T/8SU 推土算例（硬粘土、下坡 15%、45 m、槽推、一般操作手、50 min/h）==")
    f = [DOZER_CORRECTIONS["hard_to_cut_with_tilt"], dozer_grade_factor(-15),
         DOZER_CORRECTIONS["slot_dozing"], DOZER_CORRECTIONS["operator_average"],
         DOZER_CORRECTIONS["efficiency_50min"]]
    check("坡度系数 −15% → 1.30", dozer_grade_factor(-15), 1.30)
    check("600 LCY/h × 系数 × 0.87（手册把 2300/2650=0.868 取 0.87）", dozer_production(600, f + [0.87]), 405.5, 0.002)
    check("458 Lm³/h × 同一组系数（手册公制例也沿用 0.87）", dozer_production(458, f + [0.87]), 309.6, 0.002)
    print(f"       严格按 1370/1600=0.856 算则为 {dozer_production(458, f, 1600):.1f} Lm³/h（差 1.6%，是手册的取整）")
    check("曲线拟合：D8T_SU 45 m", dozer_curve_lcm_h("D8T_SU", 45.7) * YD3_PER_M3, 600, 0.15)

    print("== 4. 第 19 节：D10T2 单齿松土算例（91 m、0.9 m 间距、0.6 m 深、1.6 km/h、转向 0.25 min）==")
    r60 = ripper_production_bcm_h(91, 0.9, 0.6, 1.6, 0.25, 60)
    r45 = ripper_production_bcm_h(91, 0.9, 0.6, 1.6, 0.25, 45)
    check("每趟时间 3.66 min", r60["min_per_pass"], 3.66, 0.002)
    check("每趟方量 49.1 BCM", r60["bcm_per_pass"], 49.14, 0.001)
    check("60 min/h → 805 BCM/h", r60["bcm_h_method"], 805, 0.003)
    check("45 min/h → 604 BCM/h", r45["bcm_h_method"], 604, 0.003)
    check("实际 80% → 644", r60["bcm_h_80pct"], 644, 0.003)
    check("松动成本 $115/h ÷ 644 → $0.179/BCM", loosening_cost_per_bcm(115, 644), 0.179, 0.005)

    print("== 5. 第 7 节：320 挖掘机产量表（1.0 m³ 斗、充满 0.9、17.1 s 循环 → 190 Lm³/60 min）==")
    check("3600/17.1 × 0.9", excavator_production_lcm_h(17.1, 1.0, 0.9), 189.5, 0.01)
    check("有效循环/h：CT 0.5 min × 0.9×0.95×0.83", excavator_effective_cycles_h(0.5), 85.2, 0.01)

    print("== 6. 第 25 节：拥有成本公式（示意：P=100 万, S=20 万, 5 年, 2000 h/年, 年率 10%）==")
    o = hourly_owning_cost(1_000_000, 200_000, 5, 2000, 0.10)
    print(f"  折旧 {o['depreciation']:.1f}/h，利息保险税 {o['interest_ins_tax']:.1f}/h，合计 {o['owning']:.1f}/h")


# ---------------------------------------------------------------------------
# 本地块"情况"数字（不选型）
# ---------------------------------------------------------------------------
def site() -> None:
    print("== A. 地块坡度与 Cat 的坡度语言 ==")
    for deg in (20, 24.4, 25, 30, 34, 45):
        print(f"  {deg:5.1f}° = {slope_deg_to_pct(deg):6.1f}% 坡度")
    print("  Cat：'极端坡' 定义 >25°(47%)；静态前后倾极限 45°(100%)；推土产量图坡度轴到 ±30%。")

    print("\n== B. 运输道路：坡度 × 路面下陷 → 有效坡度 → 满载速度（740C EJ 74 t / 730C2 51.7 t，读图 ±15%）==")
    print("  路面(下陷cm)  坡度%  RR%   有效坡%  740C满载km/h  730C2满载km/h  740C空车km/h(下坡有效坡%→缓速极限)")
    for pen in (1, 3, 5):
        rr = rolling_resistance_pct(pen)
        for g in (6, 8, 10, 12, 15):
            eg = effective_grade_pct(rr, g)
            v740 = rimpull_speed_kmh("740C_EJ", 74.0, eg)
            v730 = rimpull_speed_kmh("730C2", 51.7, eg)
            fav = g - rr   # 空车下坡：有利坡 − RR
            print(f"  {pen:>6} cm     {g:>4}  {rr:4.1f}  {eg:6.1f}   {v740:8.1f}      {v730:8.1f}        (有利 {fav:4.1f}%)")
    print("  注：740C EJ 满载缓速极限（有利有效坡→km/h）：", RETARDER_740C_EJ["loaded_74t"])
    print("      重车下坡（挖方在坝线之上、往下运）时用缓速曲线；轻车上坡用牵引曲线。")

    print("\n== C. 松土：D9T/D10T2 按测距法在几种趟长下的 BCM/h（手册法偏高 10–20%，已按 80% 折）==")
    print("  趟长 m  间距 m  深 m  速度 km/h  →  BCM/h(80%)   每 BCM 松动费 @O&O $115/h")
    for L, sp, d, v in ((91, 0.9, 0.6, 1.6), (50, 0.9, 0.6, 1.6), (30, 0.9, 0.5, 1.4),
                        (91, 1.2, 0.8, 1.6), (50, 1.0, 0.7, 1.2)):
        r = ripper_production_bcm_h(L, sp, d, v, 0.25, 50)
        print(f"  {L:>6}  {sp:5.1f}  {d:4.1f}  {v:8.1f}      {r['bcm_h_80pct']:8.0f}      ${loosening_cost_per_bcm(115, r['bcm_h_80pct']):.3f}")
    print("  地震波速界限（手册图，fps）：D9T 花岗岩可松 ≤~7,000/边缘 ~8,500；片岩 ~7,500/~10,000；板岩 ~8,000/~10,000；"
          "D10T2 玄武岩 ~8,500/~9,500、花岗岩 ~8,000/~9,000。地块三种岩性均无实测波速——这是开口。")

    print("\n== D. 推土：向下坡短距离推松散岩（已松土/爆破 0.7、一般操作手 0.75、50 min 0.83、松方密度 1,960 kg/m³ 玄武岩类）==")
    print("  机型/刀    距离 m   坡度%   Lm³/h")
    for mb in ("D8T_SU", "D9T_SU", "D10T2_SU", "D9T_U"):
        for dist in (30, 60, 100):
            for g in (-15, -25):
                base = dozer_curve_lcm_h(mb, dist)
                p = dozer_production(base, [0.7, 0.75, 0.83, dozer_grade_factor(g)], 1960)
                print(f"  {mb:9s} {dist:6d}   {g:5d}   {p:7.0f}")
    print("  注：坡度系数图只画到 −30%；地块 44% 的自然坡超出图外，在图外的数字没有依据。")

    print("\n== E. 压实：W×S×L/P（Cat 夯脚 815K 4.2 m / 825K 5.3 m 双遍覆盖宽；振动光轮不在已读章节）==")
    print("  机型   宽 m  速度 km/h  层厚 mm  遍数  →  Cm³/h")
    for name, w in (("815K", 4.2), ("825K", 5.3)):
        for s, L, P in ((9.5, 100, 5), (9.5, 150, 4), (6.0, 300, 6), (6.0, 200, 8)):
            print(f"  {name}  {w:4.1f}  {s:8.1f}  {L:7d}  {P:4d}      {compactor_production_ccm_h(w, s, L, P):7.0f}")

    print("\n== F. 装运匹配（示意，不选型）：大型挖掘机 4.0 m³ 斗 × 0.85 × 30 s 循环装 740C（23 m³ 松方），路面下陷 3 cm ==")
    passes = math.ceil(23 / (4.0 * 0.85))
    load_min = passes * 0.5
    print(f"  装车 {passes} 斗 ≈ {load_min:.1f} min（手册大型挖掘机选型建议 4–6 斗，7 斗略多）")
    rr = rolling_resistance_pct(3)
    print("  (1) 料场在坝线之下：满载上坡用牵引曲线，空车下坡用缓速曲线")
    for haul_m, g in ((300, 10), (800, 10), (1500, 8), (3000, 8)):
        v_up = rimpull_speed_kmh("740C_EJ", 74.0, effective_grade_pct(rr, g))
        v_dn = retarder_speed_kmh(RETARDER_740C_EJ["empty_36t"], g - rr)
        c = truck_cycle_min(load_min, [(haul_m, v_up)], [(haul_m, v_dn)])
        n = trucks_to_match_loader(c["total"], c["load"])
        prod = 23 * 0.67 * 60 / c["total"] * 0.83
        print(f"    运距 {haul_m:5d} m 坡 {g:2d}%：满载上坡 {v_up:4.1f} km/h，空车下坡 {v_dn:4.1f} km/h，循环 {c['total']:5.1f} min，"
              f"匹配 {n:4.1f} 车/挖机，单车 {prod:5.0f} BCM/h(LF 0.67, 50 min)")
    print("  (2) 料场在坝线之上：满载下坡受缓速曲线限制，空车上坡用牵引曲线")
    for haul_m, g in ((300, 10), (800, 10), (1500, 8), (3000, 8)):
        v_dn = retarder_speed_kmh(RETARDER_740C_EJ["loaded_74t"], g - rr)
        v_up = rimpull_speed_kmh("740C_EJ", 36.0, effective_grade_pct(rr, g))
        c = truck_cycle_min(load_min, [(haul_m, v_dn)], [(haul_m, v_up)])
        n = trucks_to_match_loader(c["total"], c["load"])
        prod = 23 * 0.67 * 60 / c["total"] * 0.83
        print(f"    运距 {haul_m:5d} m 坡 {g:2d}%：满载下坡 {v_dn:4.1f} km/h，空车上坡 {v_up:4.1f} km/h，循环 {c['total']:5.1f} min，"
              f"匹配 {n:4.1f} 车/挖机，单车 {prod:5.0f} BCM/h(LF 0.67, 50 min)")
    print("  两种布置都要留着；坝线以上挖方向下运在能量上占便宜，但受制动/缓速能力约束。")

    print("\n== G. 海拔：地块 1,130–1,290 m。手册：多数机型 1,500–2,290 m 内不降功率；D8T/D9T/D10T2 T4F 到 4,600 m 100%。→ 不降额。")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--demo", action="store_true", help="复算手册算例")
    ap.add_argument("--site", action="store_true", help="本地块情况数字")
    a = ap.parse_args()
    if not (a.demo or a.site):
        a.demo = a.site = True
    if a.demo:
        demo()
    if a.site:
        print()
        site()


if __name__ == "__main__":
    main()
