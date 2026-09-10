#!/usr/bin/env python3
"""把 outputs/terrain_*_summary.json 汇总成一张情况表（Markdown + CSV）。只汇总，不评价。"""
import glob, json, os
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(os.path.dirname(HERE), "outputs")
rows = []
for f in sorted(glob.glob(os.path.join(OUT, "terrain_*_summary.json"))):
    s = json.load(open(f)); dv = s["dam_volume_m3"]; st = s["storage"] or {}; inf = s["infeasible"]
    name = os.path.basename(f).replace("terrain_", "").replace("_summary.json", "")
    rows.append({
        "case": name, "line": s["line"], "crest_ft": s["crest_ft"], "nwl_ft": s["nwl_ft"], "section": s["section_mode"],
        "floor_ft": s["floor_ft"] if s["floor_ft"] is not None else "", "L_m": round(s["length_m"]), "Hmax_m": round(s["H_axis_m"]["max"], 1),
        "L_H>100_m": sum(x["length_m"] for x in s["segments_by_height"] if x["H_range_m"][0] >= 100),
        "L_H<=30_m": sum(x["length_m"] for x in s["segments_by_height"] if x["H_range_m"][1] <= 30),
        "dam_endarea_Mm3": round(dv["end_area_10m"] / 1e6, 1), "dam_grid_Mm3": round(dv["grid_method"] / 1e6, 1),
        "rockfill_Mm3": round(dv["by_kind"].get("rockfill", 0) / 1e6, 1), "gravity_Mm3": round(dv["by_kind"].get("gravity", 0) / 1e6, 1),
        "infeasible_m": inf["length_m"],
        "cut_inside_Mm3": round(s["cut_inside_bcm"] / 1e6, 2),
        "gross_Mm3": round(st.get("V_gross_nwl_m3", float("nan")) / 1e6, 2), "gross_2mDEM_Mm3": round(st.get("V_gross_2m_dem_m3", float("nan")) / 1e6, 2),
        "us_wedge_Mm3": round(st.get("V_dam_inside_below_nwl_m3", float("nan")) / 1e6, 2), "net_Mm3": round(st.get("V_net_nwl_m3", float("nan")) / 1e6, 2),
        "A_nwl_acre": round(st.get("A_nwl_m2", float("nan")) / 4046.86, 1),
        "dam_per_net_water": round(dv["grid_method"] / st["V_net_nwl_m3"], 1) if st and st.get("V_net_nwl_m3", 0) > 0 else "",
        "net_GWh_at_518m": round(st.get("V_net_nwl_m3", 0) * 1000 * 9.81 * 518 * 0.85 / 3.6e12, 1) if st else "",
    })
cols = list(rows[0].keys())
with open(os.path.join(OUT, "terrain_cases_table.csv"), "w") as f:
    f.write(",".join(cols) + "\n")
    for r in rows: f.write(",".join(str(r[c]) for c in cols) + "\n")
with open(os.path.join(OUT, "terrain_cases_table.md"), "w") as f:
    f.write("| " + " | ".join(cols) + " |\n|" + "---|" * len(cols) + "\n")
    for r in rows: f.write("| " + " | ".join(str(r[c]) for c in cols) + " |\n")
print(open(os.path.join(OUT, "terrain_cases_table.md")).read())
