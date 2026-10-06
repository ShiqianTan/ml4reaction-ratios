# P2 报告：温度依赖反应率热力学

> 研究问题：**反应率 `r` 随温度如何变化？`ln r = −ΔΔH‡/(RT) + ΔΔS‡/R`
> 能否从多温度序列中稳定拟合，且 ΔΔH‡/ΔΔS‡ 是否能从单体结构预测？**

---

## 1. 数据与方法

- **数据**：RRMaster（12,992 行，5,742 去重），筛选 `ln r` 与 `temperature_K`
  均已知的测量。
- **序列定义**：
  - *loose*：按 (canonical_pair_id, solvent_family) 分组；
  - *strict*：再按 source_dataset 分组（同一对+溶剂+来源）。
- **拟合**：`ln r` 对 `1/T` 线性回归（scipy.stats.linregress），
  `ΔΔH‡ = −R·slope`，`ΔΔS‡ = R·intercept`。
- **门槛**：≥3 个不同温度、跨度 ≥20 K。

## 2. 结果

| 指标 | 值 |
|---|---:|
| loose 序列（对+溶剂） | 178（67 对） |
| **strict 序列（对+溶剂+来源）** | **172（66 对）** |
| 显著温度响应（p<0.05） | 92（51.7%） |
| 拟合 R² 中位数 | 0.401 |
| 温度跨度中位数 | 50 K |
| ΔΔH‡ 中位数 | 0.0 kJ/mol（IQR −7.4 … +8.5） |
| ΔΔS‡ 中位数 | −2.4 J/mol/K |
| 补偿斜率 | 326 K（r²=0.986，未中心化） |

**P2 Go/No-Go：STRONG GO**（172 strict ≥ 150）。

## 3. 关键结论

> **最强发现：温度响应不可约化（irreducible）** —— 不存在物理的 ΔΔH‡–ΔΔS‡ 补偿
> （实为拟合伪影），且 ΔΔH‡/ΔΔS‡ 无法从单体结构预测，只能靠实验测量。

1. **约一半的序列有显著温度响应**，典型 |ΔΔH‡| ≈ 5–8 kJ/mol，说明温度对
   反应率的影响普遍且不可忽略。
2. **ΔΔH‡–ΔΔS‡ 补偿是拟合伪影**（`p2_compensation.py`）：
   - 未中心化的 ΔΔH‡ vs ΔΔS‡ 回归斜率 326 K、r²=0.986，看似强补偿；
   - 但**中心化回归**（x = 1/T − mean(1/T)，可去除 OLS 斜率-截距的统计耦合）
     使相关性从 r²=0.986 **塌到 0.033**（Pearson −0.18）；
   - 补偿温度 326 K 恰好跟踪数据集平均温度 320 K —— 这正是
     Krug–Hunter–Grieger 描述的伪影特征。
   - **结论：不存在物理的 ΔΔH‡–ΔΔS‡ 补偿**；原始焓变是真实的，但未被熵抵消。
3. **负面结论：ΔΔH‡/ΔΔS‡ 无法从单体结构预测**（`p2_thermo_model.py`）：
   leave-pair-out 的 Ridge 和 RF 均 R² ≤ 0；ΔΔH‡ 的符号预测准确率 0.50（=随机）。
   即使只用 89 个显著序列，R² 仍为 −0.5。
   → **温度响应本质不可约化，必须靠实验测量**——这正是 P3 主动实验设计的动机。

## 4. 局限与下一步

- 多温度序列主要来自 Coxwell / PolyCARP，样本量与温度跨度仍偏小（中位 50 K），
  ΔΔH‡ 的个体估计误差较大。
- 「不可预测」结论基于结构特征（ECFP + 描述符）；未尝试把**溶剂/机理/压力**等
  条件作为特征，也未用物理描述符（CopDDB 势垒随温度变化）做先验。
- **下一步**：① 用 errors-in-variables 正式建模 ΔΔH‡ 不确定度；② 把温度响应
  作为 P3 主动设计中「必须测量的方向」的论据；③ 与 P2+P3 合并成闭环论文。

## 5. 复现

```bash
python3 src/p2_temperature.py     # 多温度序列拟合
python3 src/p2_compensation.py    # 补偿伪影检验
python3 src/p2_thermo_model.py    # 结构→ΔΔH‡/ΔΔS‡ 可预测性
```

输出：`projects/P2_temperature_rr/temperature_series_fits.csv`、
`p2_summary.json`、`p2_compensation_check.json`、`p2_thermo_model.json`。
