#!/usr/bin/env python3
"""Where does a dam's downstream face meet a uniformly sloping hillside?

First-order geometry used while reading USBR *Design of Small Dams* (ch. 7/8).
It is NOT a design; it only shows how strongly ground slope controls the
downstream toe position, the true structural height at the toe, and the
cross-section area per metre of dam for a ring/hillside reservoir.

Coordinate system (2-D section, x positive downstream, z positive up):
  * origin = downstream edge of the crest, crest elevation z = 0
  * ground line: z_g(x) = -H - tan(beta) * x   (falls downstream at beta)
  * H = nominal dam height = crest minus ground elevation under the
    downstream crest edge
  * downstream face:  z = -x / m_ds            (m_ds = horizontal : 1 vertical)
  * upstream face:    z = (x + b) / m_us       from the upstream crest edge
    (m_us = 0 means a vertical upstream face)
  * b = crest width

Run:  python3 hillside_toe_geometry.py [ground_slope_deg]
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass


@dataclass
class Section:
    name: str
    m_us: float  # upstream face slope, H:V (0 = vertical)
    m_ds: float  # downstream face slope, H:V
    crest_b: float  # crest width, m


def polygon_area(pts):
    a = 0.0
    n = len(pts)
    for i in range(n):
        x1, z1 = pts[i]
        x2, z2 = pts[(i + 1) % n]
        a += x1 * z2 - x2 * z1
    return abs(a) / 2.0


def solve(sec: Section, H: float, beta_deg: float):
    t = math.tan(math.radians(beta_deg))
    denom_ds = 1.0 / sec.m_ds - t
    if denom_ds <= 0:
        return None  # face flatter than ground: never daylights
    x_ds = H / denom_ds
    z_ds = -x_ds / sec.m_ds
    depth_ds = -z_ds  # true height crest -> downstream toe
    b = sec.crest_b
    if sec.m_us == 0:
        x_us = -b
        z_us = -H + t * b
    else:
        x_us = -(H + b / sec.m_us) / (1.0 / sec.m_us + t)
        z_us = (x_us + b) / sec.m_us
    depth_us = -z_us
    poly = [(-b, 0.0), (0.0, 0.0), (x_ds, z_ds), (x_us, z_us)]
    area = polygon_area(poly)
    base_len = x_ds - x_us
    return dict(x_ds=x_ds, depth_ds=depth_ds, x_us=x_us, depth_us=depth_us,
                area=area, base=base_len)


SECTIONS = [
    # Slopes quoted in Design of Small Dams ch. 7 (decked rockfill) and ch. 8
    Section("CFRD 1.4 / 1.5 (seismic-conservative rockfill)", 1.4, 1.5, 10.0),
    Section("CFRD 1.3 / 1.4 (steep rockfill, good rock)", 1.3, 1.4, 10.0),
    Section("Earth-core rockfill 2.0 / 2.0", 2.0, 2.0, 10.0),
    Section("RCC gravity vertical / 0.8", 0.0, 0.8, 6.0),
    Section("RCC gravity vertical / 0.7 (ch.8 trial section)", 0.0, 0.7, 6.0),
]

HEIGHTS = [10, 20, 30, 40, 60, 80, 100, 128]


def main():
    beta = float(sys.argv[1]) if len(sys.argv) > 1 else 24.0
    print(f"Ground slope beta = {beta:.1f} deg (tan = {math.tan(math.radians(beta)):.3f})\n")
    print("A. Toe geometry multipliers (independent of H, crest width ignored):")
    print(f"{'section':52s} {'x_toe/H':>8s} {'toe depth/H':>12s}")
    t = math.tan(math.radians(beta))
    for s in SECTIONS:
        d = 1.0 / s.m_ds - t
        if d <= 0:
            print(f"{s.name:52s} {'never':>8s} {'daylights':>12s}")
        else:
            print(f"{s.name:52s} {1/d:8.2f} {(1/s.m_ds)/d:12.2f}")
    print()
    print("B. Cross-section area per metre of dam, m^2/m  (hillside vs flat ground)")
    hdr = f"{'H (m)':>6s}"
    for s in SECTIONS:
        hdr += f" | {s.name[:22]:>22s}"
    print(hdr)
    for H in HEIGHTS:
        row = f"{H:6d}"
        for s in SECTIONS:
            r = solve(s, H, beta)
            f = solve(s, H, 0.0)
            if r is None:
                row += f" | {'no toe':>22s}"
            else:
                row += f" | {r['area']:9.0f} ({f['area']:7.0f})"
        print(row)
    print("\n   value = area on the hillside; (value) = same section on flat ground")
    print("\nC. Detail for H = 128 m (crest 4,120 ft over ground at 3,700 ft) and H = 30 m:")
    for H in (128, 30):
        print(f"  H = {H} m")
        for s in SECTIONS:
            r = solve(s, H, beta)
            if r is None:
                print(f"    {s.name:52s} downstream face never meets the slope")
                continue
            print(f"    {s.name:52s} DS toe {r['x_ds']:6.0f} m out, {r['depth_ds']:5.0f} m below crest;"
                  f" US heel {-r['x_us']:5.0f} m back, {r['depth_us']:4.0f} m below crest;"
                  f" base {r['base']:5.0f} m; area {r['area']:7.0f} m2/m")


if __name__ == "__main__":
    main()
