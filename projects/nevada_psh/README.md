# Nevada County PSH — 上库坝体设计学习与计算

地块：Nevada County, CA，APN 014-160-001-000，约 626 acre 山地。
方案：抽水蓄能（PSH）上库位于东南角高地，坝线大致沿 3,700 ft 等高线，
坝顶 4,120 ft，正常水面 4,100 ft；山坡约 24°，基岩浅。

目标：按资料清单逐份学习（每份一个笔记），然后基于真实地形逐段计算，
比较混凝土 / 堆石 / 分段组合方案，产出六项交付物。

## 目录

| 路径 | 内容 |
|---|---|
| `notes/` | 每份资料一个读书笔记：学到什么、对本地块意味着什么、可直接引用的数值准则 |
| `scripts/` | 可复算的脚本（几何、方量、库容、稳定、挖填平衡）。已有：`hillside_toe_geometry.py`（斜坡上各坝型坝趾位置与断面积）、`gravity_section_check.py`（重力坝断面按 EM 2200 / FERC 准则筛查）、`rcc_production_estimate.py`（RCC 段方量、热缝所需日浇量与拌合站规模）、`filter_gradation_limits.py`（DS-13 第 5 章反滤级配控制点：再分级、底土类别、D15F 上下限、D90F、带宽；含附录 C 算例复算）、`earthmoving_production.py`（Cat 手册 49：方量三态换算、滚动/坡度阻力与有效坡度、坡道功率、推土产量修正系数与曲线、松土测距法、挖掘机与压实机产量、卡车循环与牵引/缓速曲线读数；`--demo` 复算手册算例，`--site` 给本地块情况数字）、`compaction_coverage_and_acceptance.py`（国交省压实回数网格法的覆盖率/踏み残し与 GNSS 噪声模拟、まき出し厚换算、岩块料沉降收敛定回数、3D 出来形面管理"平均值 + 个别值 + 99.7%"判定、DEM 差分方量不确定度）、`terrain_model.py`（地形模型：任意坝线的纵剖面、横断面放坝找坝趾、分段方量（平均断面法 + 网格法）、库容曲线与净库容、库内开挖情景、可闭合性诊断、候选坝线生成）、`terrain_cases_table.py`（汇总各情况）、`ic_data_analysis.py`（FHWA 智能压实判据：CMV/MDP80 公式、压实曲线取目标 ICMV、90%/70% 覆盖判据、EU CCC 80%/10 百分位/CoV、ISSMGE 判据、四色阈值、相关系数门槛、半变异函数指数拟合、GPS 40 mm 校验） |
| `data/` | 原始地形与站点数据：`parcel_014160001000.geojson`（宗地边界，WGS84，Nevada County 宗地服务）、`parcel_014160001000_utm10.json`（同，UTM 10N）、`dem_3dep_1m_utm10.tif` + `.json`（USGS 3DEP 1 m DEM，EPSG:32610，范围与像元见 JSON）、`site_context.json`（地质单元、断层距离、地震危险性、高程/坡度统计及其来源）、`dem_3dep_2m_context_utm10.tif` + `.json`（外围 5 × 5.5 km 2 m DEM，只作看外围与核对）、`dam_line_candidates.json`（候选坝线：3,700/3,800/3,900 ft 等高线在宗地内沿边界闭合的环，及拉直版） |
| `outputs/` | 地形模型输出：`overview_hillshade_contours.png`、`context_2m_hillshade.png`（宗地及外围地形）、`terrain_<case>_{summary.json,stations.csv,plan.png,profile.png,sections.png,storage.png}`（每个坝线/断面/水位/开挖情况一套）、`terrain_cases_table.{md,csv}`（情况汇总表） |

## 资料进度

| # | 资料 | 笔记 | 状态 |
|---|---|---|---|
| 1 | USBR Design of Small Dams (1987) | `notes/01_design_of_small_dams.md` | 已读，待讨论 |
| 2 | USACE EM 1110-2-2300 (2004) | `notes/02_em_1110-2-2300.md` | 已读，待讨论 |
| 3 | USBR DS-13 (ch.1/2/3/4/9/13) | `notes/03_usbr_ds13_embankment.md` | 已读，待讨论 |
| 4 | USACE EM 1110-2-2200 (1995) + FERC 导则第 3 章重力坝 (2016) | `notes/04_em_1110-2-2200_gravity.md` | 已读，待讨论 |
| 5 | USBR RCC Manual 2017 + USACE EM 1110-2-2006 (2000) | `notes/05_usbr_rcc_manual_2017.md` | 已读，待讨论 |
| 6 | USBR Engineering Geology Field Manual (卷 I/II) + 地块地质/地震/地形初查 | `notes/06_usbr_geology_field_manual.md` | 已读，待讨论 |
| 7 | USACE EM 1110-2-2301 试验采石场与试验填筑 (1994) | `notes/07_em_1110-2-2301_test_quarries_fills.md` | 已读，待讨论 |
| 8 | USBR DS-13 ch.10 坝体施工 (2012) + ch.5 保护性反滤 (2011) | `notes/08_usbr_ds13_ch10_construction_ch5_filters.md` | 已读，待讨论 |
| 9 | Caterpillar Performance Handbook 49 (2019) | `notes/09_caterpillar_performance_handbook_49.md` | 已读，待讨论 |
| 10 | Komatsu Smart Construction（官方手册/产品页 + 日本国交省压实回数与 3D 出来形管理要领 + 独立实测） | `notes/10_komatsu_smart_construction.md` | 已读，待讨论 |
| 11 | FHWA HIF-24-097 智能压实（2024）+ HIF-17-046/037/036、HIF-13-052、HIF-12-002 通用规范 + 堆石 CCC 独立文献 | `notes/11_fhwa_intelligent_compaction.md` | 已读，待讨论 |

## 地形模型（任务 12）

- 笔记：`notes/12_terrain_model_reasoning.md`（模型给出的事实、坝趾几何定理、11 种情况的数字、以及"地形 → 坝线 → 分段高度 → 断面 → 方量与库容 → 料源 → 运输 → 设备 → 压实/验收 → 交付物"的推理链）。
- 复算：见该笔记第 6 节的命令。

## 交付物（读完后）

1. 坝型比较表
2. 整条坝线纵剖面与代表断面
3. 石料分类表
4. 挖填平衡表（原岩方 / 松方 / 压实方）
5. 施工配套表
6. 可复算的方量与库容结果（地形、参数、脚本、独立核对）
