#!/usr/bin/env python3
"""把 outputs/balanced_<line>.csv 画成图：坝顶高程 vs 净库容电量（平衡点 / 库容最大点 / 不挖）。"""
import csv, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(os.path.dirname(HERE), "outputs")

def plot_balanced(line):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]; matplotlib.rcParams["axes.unicode_minus"] = False
    rows = list(csv.DictReader(open(os.path.join(OUT, f"balanced_{line}.csv"))))
    f = lambda v: float(v) if v not in ("", None) else np.nan
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    NAME = {"cfrd": "CFRD 全环", "rcc": "RCC 全环", "combo": "组合III（1、4 段 RCC，其余 CFRD）"}
    for ty, col in (("cfrd", "green"), ("rcc", "orange"), ("combo", "purple")):
        rr = [r for r in rows if r["type"] == ty]
        if not rr: continue
        xs = [f(r["crest_ft"]) for r in rr]
        bal = [r["balanced"] == "True" for r in rr]
        gwh_nocut = [f(r["net_nocut_Mm3"]) * 1e6 * 1000 * 9.81 * f(r["head_m"]) * 0.85 / 3.6e12 for r in rr]
        ax.plot(xs, gwh_nocut, "x--", color=col, alpha=0.5, label=f"{NAME[ty]} 不挖")
        ax.plot(xs, [f(r["net_GWh"]) if b else np.nan for r, b in zip(rr, bal)], "o-", color=col, label=f"{NAME[ty]} 最少开挖的平衡点（料 = 需求）")
        ax.plot(xs, [f(r["maxnet_GWh"]) for r in rr], ":", color=col, alpha=0.8)
        ax.scatter([x for x, b in zip(xs, bal) if b], [f(r["maxnet_GWh"]) for r, b in zip(rr, bal) if b], marker="s", color=col, label=f"{NAME[ty]} 挖到库容最大（有余料）")
        ax.scatter([x for x, b in zip(xs, bal) if not b], [f(r["maxnet_GWh"]) for r, b in zip(rr, bal) if not b], marker="s", facecolors="none", edgecolors=col, label=f"{NAME[ty]} 库容最大但料不够（空心）")
        for r, b in zip(rr, bal):
            if b:
                ax.annotate(f"底{f(r['floor_ft']):.0f}/挖{f(r['cut_Mm3']):.0f}", (f(r["crest_ft"]), f(r["net_GWh"])), fontsize=7, textcoords="offset points", xytext=(4, -10))
            ax.annotate(f"底{f(r['maxnet_floor_ft']):.0f}/挖{f(r['maxnet_cut_Mm3']):.0f}/x{f(r['maxnet_surplus_ratio']):.1f}", (f(r["crest_ft"]), f(r["maxnet_GWh"])), fontsize=6, color=col, textcoords="offset points", xytext=(4, 4))
    ax.set_xlabel("坝顶高程 ft（水位 = 坝顶 − 20 ft）"); ax.set_ylabel("净库容对应电量 GWh（下库 2,400 ft，η 0.85）")
    ax.grid(alpha=0.3); ax.legend(fontsize=7, loc="upper left")
    ax.set_title(f"{line}：所有坝料来自库内开挖时的平衡设计（堆石 1 m³ 原岩 → 1.2 m³；RCC 骨料 0.9）\n标注：底=库底 ft / 挖=挖方 Mm³ / x=挖方÷需求")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f"balanced_{line}.png"), dpi=130); plt.close(fig)

if __name__ == "__main__":
    plot_balanced(sys.argv[1] if len(sys.argv) > 1 else "A_parcel_ring_3700")
