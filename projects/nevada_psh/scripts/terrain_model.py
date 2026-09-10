#!/usr/bin/env python3
"""地形模型（任务 12）：从 3DEP 1 m DEM 出发，对任意坝线做纵剖面、横断面、分段坝体方量、库容曲线、
库内开挖方量，并用独立方法互相核对。

数据（data/）：dem_3dep_1m_utm10.tif（+ .json：范围、像元、nodata）、parcel_014160001000_utm10.json（宗地边界，UTM 10N）、
dem_3dep_2m_context_utm10.tif（更大范围的 2 m 数据，只作核对与看外围）。

约定：所有坐标 UTM 10N（EPSG:32610，m）；高程 m；英尺换算 FT = 0.3048；库容 m³（另给 acre-ft）。
坝线是一条折线（可闭合）；"内侧"= 库侧，由闭合多边形的内部决定（开口线用 --inside-point 指定）。

功能：
  --build-lines          从 3,700 ft 等高线与宗地边界生成候选坝线，写 data/dam_line_candidates.json
  --run NAME             对候选坝线 NAME 计算：纵剖面、横断面、分段断面与方量、库容曲线、开挖情景、独立核对
  --crest-ft / --nwl-ft  坝顶与正常水位（默认 4120 / 4100 ft）
  --section MODE         rockfill | gravity | auto（auto：按高度分段示意：H ≤ 30 m 用重力断面，其余堆石；只是演示分段逻辑）
  --floor-ft F           库内开挖情景：把库内高于 F 的地面挖到 F（默认不挖）

断面预设（可改）：堆石 顶宽 10 m、上下游 1.6:1；重力（RCC）顶宽 6 m、上游直立、下游 0.8:1。
这些只是几何演示，不是选型；各断面的稳定性由 gravity_section_check.py 等另算。
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np
import tifffile
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import substring
import shapely

FT = 0.3048
ACRE = 4046.8564224
ACRE_FT = ACRE * FT
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "outputs")

SECTIONS = {
    "rockfill": {"crest_w": 10.0, "m_us": 1.6, "m_ds": 1.6},
    "gravity":  {"crest_w": 6.0,  "m_us": 0.0, "m_ds": 0.8},
}
AUTO_GRAVITY_MAX_H = 30.0   # m：auto 模式下 H ≤ 30 m 的段用重力断面（只是演示）


# ---------------------------------------------------------------------------
# DEM
# ---------------------------------------------------------------------------
class Terrain:
    def __init__(self, tif: str, meta_json: str):
        meta = json.load(open(meta_json))
        z = tifffile.imread(tif).astype("float64")
        z[z <= meta["nodata"] + 1] = np.nan
        self.z = z
        self.px = float(meta["pixel_size_m"][0])
        self.xmin = float(meta["extent"]["xmin"])
        self.ymax = float(meta["extent"]["ymax"])
        self.ny, self.nx = z.shape
        self.xmax = self.xmin + self.nx * self.px
        self.ymin = self.ymax - self.ny * self.px
        self.xc = self.xmin + self.px * (0.5 + np.arange(self.nx))
        self.yc = self.ymax - self.px * (0.5 + np.arange(self.ny))
        self.name = os.path.basename(tif)

    def elev(self, E, N):
        """双线性插值；越界返回 nan。E、N 可为数组。"""
        E = np.asarray(E, dtype=float); N = np.asarray(N, dtype=float)
        c = (E - self.xmin) / self.px - 0.5
        r = (self.ymax - N) / self.px - 0.5
        c0 = np.floor(c).astype(int); r0 = np.floor(r).astype(int)
        fc = c - c0; fr = r - r0
        ok = (c0 >= 0) & (c0 < self.nx - 1) & (r0 >= 0) & (r0 < self.ny - 1)
        c0c = np.clip(c0, 0, self.nx - 2); r0c = np.clip(r0, 0, self.ny - 2)
        z = self.z
        v = (z[r0c, c0c] * (1 - fr) * (1 - fc) + z[r0c + 1, c0c] * fr * (1 - fc)
             + z[r0c, c0c + 1] * (1 - fr) * fc + z[r0c + 1, c0c + 1] * fr * fc)
        return np.where(ok, v, np.nan)

    def slope_deg(self):
        dy, dx = np.gradient(self.z, self.px)
        return np.degrees(np.arctan(np.hypot(dx, dy)))

    def contour_lines(self, level_m: float):
        import contourpy
        gen = contourpy.contour_generator(x=self.xc, y=self.yc, z=self.z, name="serial")
        return [LineString(l) for l in gen.lines(level_m) if len(l) > 2]

    def mask_inside(self, poly: Polygon):
        """多边形内像元的 (rows, cols) 布尔掩膜（只在多边形外接框内计算）。"""
        minx, miny, maxx, maxy = poly.bounds
        c0 = max(int((minx - self.xmin) / self.px) - 1, 0); c1 = min(int((maxx - self.xmin) / self.px) + 2, self.nx)
        r0 = max(int((self.ymax - maxy) / self.px) - 1, 0); r1 = min(int((self.ymax - miny) / self.px) + 2, self.ny)
        X, Y = np.meshgrid(self.xc[c0:c1], self.yc[r0:r1])
        m = shapely.contains_xy(poly, X.ravel(), Y.ravel()).reshape(X.shape)
        return (slice(r0, r1), slice(c0, c1)), m


# ---------------------------------------------------------------------------
# 坝线
# ---------------------------------------------------------------------------
def load_parcel() -> Polygon:
    return Polygon(json.load(open(os.path.join(DATA, "parcel_014160001000_utm10.json")))["parcel_utm10"])


def densify(line: LineString, step: float) -> np.ndarray:
    n = max(int(math.ceil(line.length / step)), 1)
    s = np.linspace(0, line.length, n + 1)
    pts = np.array([line.interpolate(v).coords[0] for v in s])
    return pts


def close_along_boundary(line: LineString, poly: Polygon) -> LineString:
    """把端点落在多边形边界上的开口线，沿边界较短的一侧接成闭合环。"""
    ring = LineString(poly.exterior.coords)
    a, b = Point(line.coords[0]), Point(line.coords[-1])
    sa, sb = ring.project(a), ring.project(b)
    if sa > sb:
        sa, sb = sb, sa
        line = LineString(line.coords[::-1])   # 现在 line 从 sa 端走到 sb 端
    seg_in = substring(ring, sa, sb)                       # 边界上 sa→sb
    seg_out = LineString(list(substring(ring, sb, ring.length).coords) + list(substring(ring, 0, sa).coords))
    bpath = seg_in if seg_in.length < seg_out.length else seg_out
    # line: sa端 → sb端；bpath 应从 sb 端回到 sa 端
    bp = list(bpath.coords)
    if Point(bp[0]).distance(Point(line.coords[-1])) > Point(bp[-1]).distance(Point(line.coords[-1])):
        bp = bp[::-1]
    coords = list(line.coords) + bp[1:]
    if Point(coords[0]).distance(Point(coords[-1])) > 1e-6:
        coords.append(coords[0])
    return LineString(coords)


def build_candidate_lines(t: Terrain, level_ft: float = 3700.0, simplify_m: float = 2.0) -> dict:
    """候选 A：宗地内的 3,700 ft 等高线 + 沿宗地边界（东界、南界）闭合；候选 B：不受宗地约束的 3,700 ft 等高线（开口，只作对比）。"""
    P = load_parcel()
    lines = t.contour_lines(level_ft * FT)
    # 与宗地相交的等高线段中，取落在宗地内最长的一段
    best = None
    for l in lines:
        inter = l.intersection(P)
        geoms = [inter] if inter.geom_type == "LineString" else list(getattr(inter, "geoms", []))
        for g in geoms:
            if g.geom_type == "LineString" and (best is None or g.length > best.length):
                best = g
    if best is None:
        raise SystemExit("no contour segment inside parcel")
    # 简化到 5 m，避免锯齿
    best = best.simplify(simplify_m)
    ringA = close_along_boundary(best, P)
    # 候选 B：包含 best 的原始等高线（整条）
    full = max(lines, key=lambda l: l.length if l.intersects(best.buffer(1.0)) else -1)
    path = os.path.join(DATA, "dam_line_candidates.json")
    out = json.load(open(path)) if os.path.exists(path) else {"crs": "EPSG:32610"}
    out["note"] = "A_parcel_ring_<ft>: contour at <ft> inside the parcel, closed along the parcel boundary (shorter side). B_contour<ft>_open: same contour without parcel constraint (open; ends at DEM edge)."
    lv = int(round(level_ft))
    suffix = f"_s{int(simplify_m)}" if simplify_m > 2.0 else ""
    out[f"A_parcel_ring_{lv}{suffix}"] = {"closed": True, "level_ft": lv, "simplify_m": simplify_m,
                                         "coords": [list(map(float, c)) for c in ringA.coords]}
    out[f"B_contour{lv}_open"] = {"closed": False, "level_ft": lv, "coords": [list(map(float, c)) for c in full.simplify(2.0).coords],
                                  "inside_point": [686478.0, 4355908.0]}
    out.pop("A_parcel_ring", None); out.pop("B_contour3700_open", None)
    json.dump(out, open(path, "w"), indent=1)
    return out


# ---------------------------------------------------------------------------
# 剖面与断面
# ---------------------------------------------------------------------------
def longitudinal_profile(t: Terrain, line: LineString, step: float = 10.0):
    pts = densify(line, step)
    s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(pts[:, 0]), np.diff(pts[:, 1])))])
    zg = t.elev(pts[:, 0], pts[:, 1])
    return s, pts, zg


def inward_normals(pts: np.ndarray, inside_poly: Polygon | None, inside_point=None) -> np.ndarray:
    """每个站点的单位法向，指向库侧。"""
    d = np.gradient(pts, axis=0)
    d /= np.maximum(np.hypot(d[:, 0], d[:, 1])[:, None], 1e-9)
    n = np.stack([-d[:, 1], d[:, 0]], axis=1)      # 左法向
    for i in range(len(pts)):
        test = pts[i] + 5.0 * n[i]
        if inside_poly is not None:
            if not inside_poly.contains(Point(test)):
                n[i] = -n[i]
        else:
            # 开口线：指向 inside_point 的一侧
            v = np.array(inside_point) - pts[i]
            if np.dot(v, n[i]) < 0:
                n[i] = -n[i]
    return n


def cross_section(t: Terrain, p: np.ndarray, n: np.ndarray, half_width: float = 900.0, step: float = 1.0):
    """沿法向取地面线：offset < 0 = 库外（下游），> 0 = 库内（上游）。"""
    off = np.arange(-half_width, half_width + step, step)
    E = p[0] + off * n[0]; N = p[1] + off * n[1]
    return off, t.elev(E, N)


def dam_section(off, zg, crest: float, spec: dict, nwl: float | None = None) -> dict:
    """在一条横断面上放坝：坝轴在 offset 0，坝顶高程 crest；上游（库侧，offset>0）坡 m_us，下游 m_ds。
    返回面积、两侧坝趾、最大结构高度、上游侧水位以下的坝体面积（算库容时扣除）。"""
    b = spec["crest_w"]; mus = spec["m_us"]; mds = spec["m_ds"]
    zc = np.interp(0.0, off, zg)
    if not np.isfinite(zc) or zc >= crest:
        return {"area": 0.0, "H_axis": max(crest - zc, 0.0) if np.isfinite(zc) else float("nan"),
                "toe_us": None, "toe_ds": None, "H_max": 0.0, "area_us_below_nwl": 0.0, "note": "ground at/above crest"}
    # 坝面高程函数
    def zface(o):
        a = np.abs(o)
        if a <= b / 2:
            return crest
        m = mus if o > 0 else mds
        if m == 0:
            return -np.inf   # 直立面：越过坝顶边缘就没有坝体
        return crest - (a - b / 2) / m
    zf = np.array([zface(o) for o in off])
    inside = zf > zg
    # 从轴向两侧找第一次坝面低于地面的位置（坝趾）
    i0 = int(np.argmin(np.abs(off)))
    def toe(direction):
        i = i0
        while 0 <= i < len(off) and inside[i]:
            i += direction
        if i < 0 or i >= len(off):
            return None
        return i
    ius, ids = toe(+1), toe(-1)
    if ius is None or ids is None:
        return {"area": float("nan"), "H_axis": crest - zc, "toe_us": None if ius is None else float(off[ius]),
                "toe_ds": None if ids is None else float(off[ids]), "H_max": float("nan"), "area_us_below_nwl": float("nan"),
                "note": "toe not found within half-width (ground slope >= face slope?)"}
    sl = slice(ids, ius + 1)
    thick = np.maximum(zf[sl] - zg[sl], 0.0)
    d = float(off[1] - off[0])
    area = float(np.sum(thick) * d)
    hmax = float(np.nanmax(thick))
    a_us = 0.0
    if nwl is not None:
        slu = slice(i0, ius + 1)
        a_us = float(np.sum(np.maximum(np.minimum(zf[slu], nwl) - zg[slu], 0.0)) * d)
    return {"area": area, "H_axis": float(crest - zc), "toe_us": float(off[ius]), "toe_ds": float(off[ids]),
            "H_max": hmax, "area_us_below_nwl": a_us, "z_toe_us": float(zg[ius]), "z_toe_ds": float(zg[ids]),
            "base_width": float(off[ius] - off[ids]), "note": ""}


def flattest_closing_slope(off, zg, crest: float, crest_w: float, side: int, ms=(1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2)):
    """在给定断面上，从缓到陡试坡比 m（水平:竖直），返回第一个能在半宽内与地面相交的 m；None = 直立面也不闭合。
    side=-1 下游（库外），+1 上游（库内）。"""
    i0 = int(np.argmin(np.abs(off)))
    for m in ms:
        i = i0
        while 0 <= i < len(off):
            a = abs(off[i])
            zf = crest if a <= crest_w / 2 else crest - (a - crest_w / 2) / m
            if zf <= zg[i]:
                return m
            i += side
    return None


def choose_spec(mode: str, H_axis: float) -> tuple[str, dict]:
    if mode == "auto":
        return ("gravity", SECTIONS["gravity"]) if H_axis <= AUTO_GRAVITY_MAX_H else ("rockfill", SECTIONS["rockfill"])
    return mode, SECTIONS[mode]


# ---------------------------------------------------------------------------
# 库容与开挖
# ---------------------------------------------------------------------------
def storage_curve(t: Terrain, poly: Polygon, levels_m: np.ndarray, z_override=None):
    win, m = t.mask_inside(poly)
    z = (t.z if z_override is None else z_override)[win]
    zin = z[m]
    zin = zin[np.isfinite(zin)]
    cell = t.px ** 2
    areas = np.array([(zin < L).sum() * cell for L in levels_m])
    vols = np.array([np.sum(np.maximum(L - zin, 0.0)) * cell for L in levels_m])
    return areas, vols, zin


def excavate(t: Terrain, poly: Polygon, floor_m: float):
    """把多边形内高于 floor 的地面削到 floor；返回新 DEM 与挖方（原岩方）。"""
    win, m = t.mask_inside(poly)
    z2 = t.z.copy()
    sub = z2[win]
    cut = np.where(m & (sub > floor_m), sub - floor_m, 0.0)
    sub[m & (sub > floor_m)] = floor_m
    z2[win] = sub
    return z2, float(np.nansum(cut) * t.px ** 2)


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def run(name: str, crest_ft: float, nwl_ft: float, mode: str, floor_ft: float | None, step: float = 10.0, tag: str = ""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = Terrain(os.path.join(DATA, "dem_3dep_1m_utm10.tif"), os.path.join(DATA, "dem_3dep_1m_utm10.json"))
    cand = json.load(open(os.path.join(DATA, "dam_line_candidates.json")))[name]
    line = LineString(cand["coords"])
    closed = cand["closed"]
    poly = Polygon(cand["coords"]) if closed else None
    inside_pt = cand.get("inside_point")
    crest = crest_ft * FT; nwl = nwl_ft * FT
    os.makedirs(OUT, exist_ok=True)
    oname = name + (("_" + tag) if tag else "")

    # 开挖情景（先改 DEM，再算一切）
    cut_bcm = 0.0
    if floor_ft is not None and poly is not None:
        t.z, cut_bcm = excavate(t, poly, floor_ft * FT)

    # 1. 纵剖面
    s, pts, zg = longitudinal_profile(t, line, step)
    H_axis = crest - zg
    nrm = inward_normals(pts, poly, inside_pt)

    # 2. 逐站断面
    rows = []
    for i in range(len(pts)):
        off, zsec = cross_section(t, pts[i], nrm[i])
        kind, spec = choose_spec(mode, float(H_axis[i]))
        d = dam_section(off, zsec, crest, spec, nwl)
        d.update({"station": float(s[i]), "E": float(pts[i, 0]), "N": float(pts[i, 1]), "z_ground": float(zg[i]), "kind": kind,
                  "feasible": bool(np.isfinite(d["area"]))})
        if not d["feasible"]:
            d["m_ds_closing"] = flattest_closing_slope(off, zsec, crest, spec["crest_w"], -1)
            d["m_us_closing"] = flattest_closing_slope(off, zsec, crest, spec["crest_w"], +1)
        rows.append(d)

    # 3. 分段方量（平均断面法）
    def end_area_volume(rows, key="area"):
        v = 0.0
        for a, b in zip(rows, rows[1:]):
            A1, A2 = a[key], b[key]
            if not (np.isfinite(A1) and np.isfinite(A2)):
                continue
            v += 0.5 * (A1 + A2) * (b["station"] - a["station"])
        return v
    V_dam = end_area_volume(rows)
    V_us_below = end_area_volume(rows, "area_us_below_nwl")
    by_kind = {}
    for a, b in zip(rows, rows[1:]):
        if np.isfinite(a["area"]) and np.isfinite(b["area"]):
            by_kind[a["kind"]] = by_kind.get(a["kind"], 0.0) + 0.5 * (a["area"] + b["area"]) * (b["station"] - a["station"])
    missing = sum(1 for r in rows if not np.isfinite(r["area"]))

    # 4. 库容（闭合线才算）
    storage = None
    if poly is not None:
        levels = np.arange(math.floor(np.nanmin(t.elev(pts[:, 0], pts[:, 1]))) - 5, crest + 1, 1.0)
        areas, vols, zin = storage_curve(t, poly, levels)
        V_gross = float(np.interp(nwl, levels, vols)); A_nwl = float(np.interp(nwl, levels, areas))
        # 核对 1：面积-高程曲线梯形积分
        V_trap = float(np.trapezoid(areas[levels <= nwl], levels[levels <= nwl]))
        # 核对 2：2 m 外围 DEM
        try:
            t2 = Terrain(os.path.join(DATA, "dem_3dep_2m_context_utm10.tif"), os.path.join(DATA, "dem_3dep_2m_context_utm10.json"))
            if floor_ft is not None:
                t2.z, _ = excavate(t2, poly, floor_ft * FT)
            _, v2, _ = storage_curve(t2, poly, np.array([nwl]))
            V_2m = float(v2[0])
        except Exception as e:  # noqa
            V_2m = float("nan")
        storage = {"levels_m": levels.tolist(), "areas_m2": areas.tolist(), "vols_m3": vols.tolist(),
                   "V_gross_nwl_m3": V_gross, "V_gross_trapezoid_m3": V_trap, "V_gross_2m_dem_m3": V_2m,
                   "V_dam_inside_below_nwl_m3": V_us_below, "V_net_nwl_m3": V_gross - V_us_below,
                   "A_nwl_m2": A_nwl, "ring_area_m2": float(poly.area),
                   "ground_inside_min_m": float(zin.min()), "ground_inside_max_m": float(zin.max()),
                   "area_inside_above_nwl_m2": float((zin >= nwl).sum() * t.px ** 2)}

    # 5. 独立核对：坝体方量的网格法（最近站的断面几何铺到每个像元）
    from scipy.spatial import cKDTree
    dense = densify(line, 2.0)
    dn = inward_normals(dense, poly, inside_pt)
    tree = cKDTree(dense)
    buf = line.buffer(950.0)
    win, m = t.mask_inside(buf)
    X, Y = np.meshgrid(t.xc[win[1]], t.yc[win[0]])
    zz = t.z[win]
    sel = m & np.isfinite(zz)
    dd, idx = tree.query(np.column_stack([X[sel], Y[sel]]))
    vec = np.column_stack([X[sel], Y[sel]]) - dense[idx]
    offv = np.einsum("ij,ij->i", vec, dn[idx])            # 带号偏移：>0 库内
    # 每个像元用最近站的断面类型
    st_dense = np.concatenate([[0], np.cumsum(np.hypot(np.diff(dense[:, 0]), np.diff(dense[:, 1])))])
    H_dense = crest - t.elev(dense[:, 0], dense[:, 1])
    kinds = np.array([choose_spec(mode, float(h))[0] for h in H_dense])
    zf = np.full(offv.shape, -np.inf)
    for kname, spec in SECTIONS.items():
        k = kinds[idx] == kname
        a = np.abs(offv[k]); b = spec["crest_w"]
        mm = np.where(offv[k] > 0, spec["m_us"], spec["m_ds"])
        with np.errstate(divide="ignore", invalid="ignore"):
            face = np.where(a <= b / 2, crest, np.where(mm > 0, crest - (a - b / 2) / np.where(mm > 0, mm, 1), -np.inf))
        zf[k] = face
    thick = np.maximum(zf - zz[sel], 0.0)
    # 只计入"从轴连续"的坝体：用断面法的坝趾范围限制（近似：偏移在各站坝趾之间）
    toe_us = np.array([r["toe_us"] if r["toe_us"] is not None else np.nan for r in rows])
    toe_ds = np.array([r["toe_ds"] if r["toe_ds"] is not None else np.nan for r in rows])
    st_rows = np.array([r["station"] for r in rows])
    feas = np.array([r["feasible"] for r in rows], dtype=float)
    feas_dense = np.interp(st_dense[idx], st_rows, feas) >= 0.999     # 只在可闭合站之间计
    tu = np.interp(st_dense[idx], st_rows, np.nan_to_num(toe_us, nan=0.0))
    td = np.interp(st_dense[idx], st_rows, np.nan_to_num(toe_ds, nan=0.0))
    within = (offv <= tu) & (offv >= td) & feas_dense
    V_grid = float(np.sum(thick[within]) * t.px ** 2)
    infeasible_len = float(sum(1 for r in rows if not r["feasible"]) * step)
    closing = [r.get("m_ds_closing") for r in rows if not r["feasible"]]
    closing_hist = {str(k): closing.count(k) * step for k in sorted(set(closing), key=lambda v: (v is None, v))}

    # 6. 分段汇总（按 25 m 站）
    Hs = np.array([r["H_axis"] for r in rows])
    seg = []
    for lo, hi in ((0, 10), (10, 20), (20, 30), (30, 50), (50, 80), (80, 100), (100, 130), (130, 999)):
        k = (Hs > lo) & (Hs <= hi)
        if k.any():
            seg.append({"H_range_m": [lo, hi], "length_m": float(k.sum() * step),
                        "volume_m3": float(sum(0.5 * (rows[i]["area"] + rows[i + 1]["area"]) * step
                                              for i in np.where(k)[0] if i + 1 < len(rows)
                                              and np.isfinite(rows[i]["area"]) and np.isfinite(rows[i + 1]["area"])))})

    summary = {
        "line": name, "closed": closed, "length_m": float(s[-1]), "crest_ft": crest_ft, "nwl_ft": nwl_ft, "section_mode": mode,
        "floor_ft": floor_ft, "cut_inside_bcm": cut_bcm,
        "ground_along_axis_ft": {"min": float(np.nanmin(zg) / FT), "max": float(np.nanmax(zg) / FT)},
        "H_axis_m": {"max": float(np.nanmax(Hs)), "mean": float(np.nanmean(Hs)), "length_H_gt_0": float((Hs > 0).sum() * step)},
        "dam_volume_m3": {"end_area_10m": V_dam, "grid_method": V_grid, "by_kind": by_kind, "stations_without_toe": missing,
                          "note": "end-area ignores plan curvature (underestimates on convex outer side); grid method assigns each cell to the nearest axis point. Both exclude stations where the face never meets the ground."},
        "infeasible": {"length_m": infeasible_len, "flattest_closing_m_ds_length_m": closing_hist,
                       "meaning": "stations where the chosen downstream face slope never meets the ground within 900 m; histogram gives the flattest slope (H:V) that would close, None = even a vertical face does not close"},
        "segments_by_height": seg,
        "storage": storage,
    }
    # 5 m 与 25 m 站距的平均断面法核对
    for st2 in (5.0, 25.0):
        s2, p2, z2 = longitudinal_profile(t, line, st2)
        n2 = inward_normals(p2, poly, inside_pt)
        areas2 = []
        for i in range(len(p2)):
            off, zsec = cross_section(t, p2[i], n2[i])
            kind, spec = choose_spec(mode, float(crest - z2[i]))
            areas2.append(dam_section(off, zsec, crest, spec, nwl)["area"])
        rr = [{"area": a, "station": float(v)} for a, v in zip(areas2, s2)]
        summary["dam_volume_m3"][f"end_area_{int(st2)}m"] = end_area_volume(rr)

    json.dump(summary, open(os.path.join(OUT, f"terrain_{oname}_summary.json"), "w"), indent=1)
    # 站表 CSV
    with open(os.path.join(OUT, f"terrain_{oname}_stations.csv"), "w") as f:
        f.write("station_m,E,N,z_ground_ft,H_axis_m,kind,area_m2,H_max_m,toe_us_m,toe_ds_m,base_width_m,area_us_below_nwl_m2,note\n")
        for r in rows:
            f.write(f"{r['station']:.0f},{r['E']:.1f},{r['N']:.1f},{r['z_ground']/FT:.0f},{r['H_axis']:.1f},{r['kind']},"
                    f"{r['area'] if np.isfinite(r['area']) else ''},{r.get('H_max','')},{r['toe_us'] if r['toe_us'] is not None else ''},"
                    f"{r['toe_ds'] if r['toe_ds'] is not None else ''},{r.get('base_width','')},{r['area_us_below_nwl']},{r['note']}\n")

    # 图 1：平面 + 图 2：纵剖面 + 图 3：代表断面 + 图 4：库容曲线
    fig, ax = plt.subplots(figsize=(10, 8))
    sl = t.slope_deg()
    win, mk = t.mask_inside(line.buffer(900))
    ax.imshow(sl[win], cmap="Greys", extent=[t.xc[win[1]][0] - 0.5, t.xc[win[1]][-1] + 0.5, t.yc[win[0]][-1] - 0.5, t.yc[win[0]][0] + 0.5], vmin=0, vmax=50)
    zz = t.z
    ax.contour(t.xc, t.yc, zz / FT, levels=np.arange(2600, 4400, 100), colors="y", linewidths=0.3)
    ax.contour(t.xc, t.yc, zz / FT, levels=[nwl_ft], colors="red", linewidths=0.8)
    P = load_parcel(); px, py = P.exterior.xy; ax.plot(px, py, "w-", lw=1.2, label="parcel")
    hc = ax.scatter(pts[:, 0], pts[:, 1], c=Hs, cmap="viridis", s=6, vmin=0, vmax=max(1, np.nanmax(Hs)))
    plt.colorbar(hc, ax=ax, label="H at axis (m)")
    # 坝趾
    for r in rows[::3]:
        i = int(r["station"] / step)
        for key, col in (("toe_us", "cyan"), ("toe_ds", "magenta")):
            if r[key] is not None:
                q = pts[i] + r[key] * nrm[i]; ax.plot(q[0], q[1], ".", color=col, ms=2)
    ax.set_xlim(t.xc[win[1]][0], t.xc[win[1]][-1]); ax.set_ylim(t.yc[win[0]][-1], t.yc[win[0]][0])
    ax.set_aspect("equal"); ax.set_title(f"{name}: axis colored by H; cyan=upstream toe, magenta=downstream toe; red={nwl_ft:.0f} ft")
    ax.legend(loc="lower left"); fig.tight_layout(); fig.savefig(os.path.join(OUT, f"terrain_{oname}_plan.png"), dpi=120); plt.close(fig)

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.plot(s, zg / FT, "k-", lw=1, label="ground along axis")
    ax.axhline(crest_ft, color="r", lw=1, label=f"crest {crest_ft:.0f} ft"); ax.axhline(nwl_ft, color="b", lw=0.8, ls="--", label=f"NWL {nwl_ft:.0f} ft")
    zt_us = np.array([r.get("z_toe_us", np.nan) for r in rows]); zt_ds = np.array([r.get("z_toe_ds", np.nan) for r in rows])
    ax.plot(s, zt_us / FT, "c.", ms=2, label="upstream toe elev"); ax.plot(s, zt_ds / FT, "m.", ms=2, label="downstream toe elev")
    kinds_arr = np.array([r["kind"] for r in rows])
    for kname, col in (("gravity", "orange"), ("rockfill", "green")):
        k = kinds_arr == kname
        if k.any():
            ax.fill_between(s, zg / FT, crest_ft, where=k, color=col, alpha=0.25, label=f"{kname} segments")
    ax.set_xlabel("station along axis (m)"); ax.set_ylabel("elevation (ft)"); ax.legend(loc="lower right", fontsize=8)
    ax.set_title(f"{name}: longitudinal profile (L={s[-1]:.0f} m, Hmax={np.nanmax(Hs):.0f} m)")
    ax.grid(alpha=0.3); fig.tight_layout(); fig.savefig(os.path.join(OUT, f"terrain_{oname}_profile.png"), dpi=120); plt.close(fig)

    # 代表断面：H 最大、H 中位、H≈20 m、H≈50 m
    picks = []
    order = np.argsort(Hs)
    for target in (np.nanmax(Hs), np.nanmedian(Hs), 20.0, 50.0):
        i = int(np.nanargmin(np.abs(Hs - target))); picks.append(i)
    picks = sorted(set(picks))
    fig, axs = plt.subplots(len(picks), 1, figsize=(12, 3.2 * len(picks)))
    axs = np.atleast_1d(axs)
    for ax, i in zip(axs, picks):
        off, zsec = cross_section(t, pts[i], nrm[i])
        kind, spec = choose_spec(mode, float(Hs[i]))
        d = dam_section(off, zsec, crest, spec, nwl)
        ax.plot(off, zsec / FT, "k-", lw=1)
        if d["toe_us"] is not None and d["toe_ds"] is not None:
            b = spec["crest_w"]
            xs = [d["toe_ds"], -b / 2, b / 2, d["toe_us"]]
            zs = [d["z_toe_ds"] / FT, crest_ft, crest_ft, d["z_toe_us"] / FT]
            ax.fill(xs + [d["toe_us"], d["toe_ds"]], zs + [d["z_toe_us"] / FT, d["z_toe_ds"] / FT], color="orange" if kind == "gravity" else "green", alpha=0.35)
            ax.plot(xs, zs, "r-", lw=1)
        ax.axhline(nwl_ft, color="b", ls="--", lw=0.7)
        ax.set_title(f"station {s[i]:.0f} m  H_axis={Hs[i]:.0f} m  {kind}  area={d['area']:.0f} m²  base={d.get('base_width', float('nan')):.0f} m  {d['note']}", fontsize=9)
        ax.set_xlabel("offset from axis (m): + = reservoir side"); ax.set_ylabel("ft"); ax.grid(alpha=0.3)
        ax.set_xlim(-400, 300)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f"terrain_{oname}_sections.png"), dpi=120); plt.close(fig)

    if storage is not None:
        fig, ax1 = plt.subplots(figsize=(8, 5))
        lv = np.array(storage["levels_m"]) / FT
        ax1.plot(np.array(storage["vols_m3"]) / 1e6, lv, "b-", label="storage (10^6 m³)")
        ax1.set_xlabel("gross storage inside ring (10^6 m³)"); ax1.set_ylabel("water level (ft)")
        ax2 = ax1.twiny(); ax2.plot(np.array(storage["areas_m2"]) / ACRE, lv, "g--", label="area (acre)"); ax2.set_xlabel("water area (acre)")
        ax1.axhline(nwl_ft, color="r", lw=0.8); ax1.grid(alpha=0.3); ax1.set_title(f"{name}: elevation-area-storage (gross, before dam-body deduction)")
        fig.tight_layout(); fig.savefig(os.path.join(OUT, f"terrain_{oname}_storage.png"), dpi=120); plt.close(fig)

    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-lines", action="store_true")
    ap.add_argument("--level-ft", type=float, default=3700.0, help="contour level for --build-lines")
    ap.add_argument("--simplify-m", type=float, default=2.0, help="Douglas-Peucker tolerance (m) to straighten the contour into a dam line before closing")
    ap.add_argument("--run", type=str)
    ap.add_argument("--crest-ft", type=float, default=4120.0)
    ap.add_argument("--nwl-ft", type=float, default=4100.0)
    ap.add_argument("--section", choices=["rockfill", "gravity", "auto"], default="auto")
    ap.add_argument("--floor-ft", type=float, default=None)
    ap.add_argument("--step", type=float, default=10.0)
    ap.add_argument("--tag", type=str, default="", help="suffix for output file names")
    a = ap.parse_args()
    if a.build_lines:
        t = Terrain(os.path.join(DATA, "dem_3dep_1m_utm10.tif"), os.path.join(DATA, "dem_3dep_1m_utm10.json"))
        c = build_candidate_lines(t, a.level_ft, a.simplify_m)
        for k, v in c.items():
            if isinstance(v, dict) and "coords" in v:
                print(k, "closed" if v["closed"] else "open", f"{LineString(v['coords']).length:.0f} m")
    if a.run:
        sm = run(a.run, a.crest_ft, a.nwl_ft, a.section, a.floor_ft, a.step, a.tag)
        print(json.dumps({k: v for k, v in sm.items() if k != "storage"}, indent=1, ensure_ascii=False))
        if sm["storage"]:
            st = sm["storage"]
            print("storage:", {k: v for k, v in st.items() if not isinstance(v, list)})


if __name__ == "__main__":
    main()
