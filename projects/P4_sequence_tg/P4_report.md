# P4 报告：反应率(RR) → 序列 → 玻璃化转变温度(Tg) 多尺度链路

> 研究问题：**单体反应率 `(r1, r2)` 能否通过 Mayo-Lewis 序列分布传导到共聚物 Tg，
> 并给出可预测的定量修正？**
>
> 本文档汇总 P4 的全部主要实验与结论。数据与代码见本仓库
> `src/p4_*.py`、`projects/P4_sequence_tg/`。

---

## 1. 目标与总体框架

P4 把「反应率」这一动力学量接到「材料性质」上，构成三段链路：

```
(r1, r2)  ──Mayo-Lewis──▶  序列 (二元组/嵌段度/交替度)  ──二元组可加性──▶  Tg
```

零阶基线是**组成规则**（Fox / Gordon-Taylor，只依赖共聚物组成，与序列无关）。
P4 的核心问题是：**在组成之上，序列信息（由 RR 决定）还能贡献多少？**
答案分两半：

- 对**普通无规自由基共聚**（r1·r2 接近 1），序列 ≈ 理想无规，**RR 几乎不贡献**；
- 对**强交替/强序列**体系，RR 通过二元组分布**显著**修正 Tg，且方向是**体系特异**的。

---

## 2. 数据来源

| 数据 | 内容 | 用途 |
|---|---|---|
| RRMaster（本仓库） | 12,992 原始 → 5,742 去重，**2,633 单体对**，1,218 单体 | 提供 `(r1, r2)` |
| PolyInfo 均聚物 Tg | 11,380 条均聚物 Tg + SMILES | 均聚物 Tg 模型、Tg_AA/Tg_BB |
| `DeltaLearningCopolymerTg/data` | tgvalue 351 行 + styrenic 48 行 | 近无规自由基共聚 Tg（真实） |
| `SequenceTg/data`（sibling repo） | Johnston 1973/76 31 样本、Barton 1970、NMR 二元组 | 强序列自由基体系（真实） |

---

## 3. 实验与结论

### 实验 1 — Mayo-Lewis 序列引擎 + 蒙特卡洛校验（`p4_sequence.py`）

实现末端模型，由 `(r1, r2, f1)` 计算瞬时组成 `F1` 与序列描述符
（`P11/P12/P21/P22`、二元组 `dyad_11/12/21/22`、嵌段度、交替度、链长 `L1/L2`、
序列熵）。

**蒙特卡洛校验**（r1=1, r2=2, f1=0.5）：

| 量 | 解析解 | MC |
|---|---:|---:|
| F1 | 0.4000 | 0.3997 |
| P11 | 0.5000 | 0.5015 |
| P22 | 0.6667 | 0.6681 |

✅ 引擎正确。已对 **2,633 个真实 RR 对**生成 50,046 行组成-序列描述符。

### 实验 2 — RR → 序列：Mayo-Lewis 对照 NMR 二元组（`p4_sequence_tg_validation.py`）

6 个有 NMR 二元组实测的体系（AN-EA、BVE-CEVE、NVP-CEVE、S-Acrylamide、MMA-S、VDF-HFP）。

**结论：RR → 序列这一半链路在真实数据上成立，二元组 MAE = 0.055**（杂二元组 0.043）。

| 体系 | NMR 杂二元组 | 预测 | 偏差 |
|---|---:|---:|---:|
| AN-EA | 0.48 | 0.485 | +0.005 |
| S-Acrylamide | 0.40 | 0.402 | +0.002 |
| NVP-CEVE | 0.18 | 0.196 | +0.016 |
| VDF-HFP | 0.11 | 0.092 | −0.018 |
| MMA-S（真无规） | 0.50 | 0.671 | +0.171 |

唯一系统性偏差是 MMA-St「无规」——该样品由 RAFT 制得、真正无规，而
`r_MMA·r_S=0.24` 暗示轻微交替；这说明**当样品实际序列偏离 RR 预测时，误差方向可解释**，
进一步佐证（而非反驳）链路有效性。

### 实验 3 — RR → 序列 → Tg：Barton 二元组可加性 vs Fox（Johnston 5 体系，31 样本）

Barton 模型 `Tg = n_AA·Tg_AA + n_BB·Tg_BB + (n_AB+n_BA)·Tg_AB`，二元组由
Mayo-Lewis（有 RR 时）或理想无规（无 RR 时）给出。

**描述性结果：二元组模型在这 31 个样本上的 Fox 误差约减半；这不是独立体系预测验证。**

该对比使用了体系级 `Tg_AB` 参数。两条标为 alternating 的样品直接预测为输入 `Tg_AB`，观测 Tg 也等于该参数，因而两条样品误差为零；其余 random 样品仍使用同一体系提供的 `Tg_AB`。因此 15.3→7.8 K 是参数已知时的样本级比较，不能作为新体系的端到端预测性能。按类型拆分及补充分析见 `P4_supplementary_experiments.md`。

| | MAE (K) | R² |
|---|---:|---:|
| Fox（仅组成） | 15.3 | 0.754 |
| Barton 二元组 | **7.8** | **0.881** |

分体系：

| 体系 | Fox MAE | 二元组 MAE |
|---|---:|---:|
| S-AN | 6.6 K | **0.9 K** |
| AN-MMA | 13.6 K | **1.6 K** |
| aMS-AN | 25.0 K | **4.5 K** |
| BMA-VCl | 11.1 K | 7.7 K |
| MMA-VCl | 20.4 K | 22.2 K（变差，见实验 4） |

### 实验 4 — Barton 符号规则（交替二元组的方向）

交替二元组 Tg（Tg_AB）相对均聚物均值 `(Tg_AA+Tg_BB)/2` 的偏移决定序列效应的**符号**：

| 体系 | Tg_AB − 均值 (K) | 方向 |
|---|---:|---|
| MMA-VCl | −41.4 | 交替**降** Tg |
| BMA-VCl | −35.9 | 交替**降** Tg |
| AN-MMA | −26.0 | 交替**降** Tg |
| aMS-AN | −19.0 | 交替**降** Tg |
| S-AN | **+9.0** | 交替**升** Tg |

**结论：交替对 Tg 的作用方向是体系特异的（4/5 降、1/5 升）。**
这解释了 MMA-VCl 在实验 3 中变差——其 Tg_AB 远低于均聚均值，单一 Tg_AB 的
二元组模型会过度修正低 MMA 组成的样品。它同时也解释了为什么「交替→升 Tg」的
通用启发式会失效。

### 实验 5 — 近无规真实共聚物：序列几乎不贡献（`p4_real_sequence_tg.py`）

`DeltaLearningCopolymerTg` 的 210 个测量（40 对，r1·r2 中位数 0.53，绝大多数
0.2–0.9）上做 leave-one-pair-out 的 Δ-learning：

| 模型 | MAE (K) | R² |
|---|---:|---:|
| Fox（摩尔组成）基线 | 8.60 | 0.925 |
| Δ（仅组成） | 10.0 | 0.902 |
| Δ（+ 反应率） | 9.90 | 0.896 |
| Δ（+ 序列描述符） | 10.36 | 0.900 |
| Δ（+ rr + 结构 ECFP） | **8.14** | **0.923** |

**结论（零阶）：对普通无规自由基共聚物，Fox 已解释 ~92% 方差，RR 序列信息
无额外预测力**（`r1·r2 ↔ 残差` 相关仅 0.06；唯一增益来自单体 ECFP 结构，反映
单体特异性化学如苯乙烯-丙烯酰胺 +25 K 的氢键效应）。这正是实验 3 的另一端——
在近无规区，二元组 ≈ 理想无规，Fox ≈ Barton。

### 实验 6 — 均聚物 Tg 模型 + Δ-learning 方法论自检（`p4_tg_model.py`, `p4_delta_learning.py`）

- **均聚物 Tg 从结构可预测且可外推**：RF 随机划分 MAE 21.2 °C / R² 0.887；
  scaffold 化学划分 MAE 34.4 °C / R² 0.752。
- **Δ-learning 能恢复已知的序列信号**（受控模拟，γ·(嵌段度−1)·ΔT 生成目标）：
  Fox 基线 4.19 °C → Δ-RF **3.02 °C**（残差 R² 0.51），且学到的 top 特征正是
  alternation/blockiness/P21/P22/dyad_22。

---

## 4. 总纲（P4 定论）

1. **RR→序列有初步真实数据校验**：选取的 6 个 NMR 体系中，dyad MAE 为 0.055；
   这支持序列引擎在这些案例上的可行性，但样本数和来源追溯仍需补强。
2. **Johnston 子集尚未独立验证 RR→序列→Tg**：当前 Barton 与 Fox 的误差比较使用了
   已知体系 `Tg_AB`，且没有按化学体系留出。该结果只能作为经典模型的样本级描述性比较。
3. **近无规数据未显示稳健的序列预测增益**：210 个测量、40 个单体对上，序列模型
   没有稳定优于 Fox；但数据主要覆盖近无规体系，组成也不是已知进料，不能据此断言
   序列效应普遍不存在。
4. **候选论文主张**：评估 RR 导出的序列统计何时能在组成基线之外改善 Tg 预测，并报告
   这种增益的适用范围与不确定性。该主张需要独立体系留出、来源级验证和更可靠的组成输入。

### 最新补充：RR 序列的正向案例与失效边界

- 新近文献筛查发现 Cipullo 等 2026 年 E/V 工作：论文报告 Sample 1/3 的 V 含量分别约 17.9%/18.1%；按 Table 2 已舍入 dyad 反算则为 18.05%/18.25%，应保留为不同来源的数值。其 `[EV]` 为 35.1%/26.5%，`[VV]` 为 0.5%/5.0%；Tg 约 −23/−14 °C，Mn 为 268/401 kDa。论文报告相近组成、可比分子量而序列不同的样品 Tg 相差约 10 °C。分子量并非完全相同，因此这是重要关联证据而非完全隔离的因果效应。Table 4 还报告了逐样品 RR：C1 Sample 1 `rE=19,rV=0.003`；C2 Samples 2/3 分别为 `14,0.14` 和 `16,0.11`。SI 条件支持用 Sample 2→3 做序列检验；现已完成以 Sample 3 报告组成作为输入的预测。论文自身 RR 是从样品 NMR 序列推得。[开放获取原文](https://doi.org/10.1016/j.polymer.2026.130510)
- 继续追查 Brar 与 Sunita 1992 年 AN–EA 论文：摘要明确同时有 EVM RR、NMR 三元/二元序列和 DSC Tg，且称 RR 预测的 triad 与实测值吻合；但出版商当前可见页面只有摘要，未找到可复核的逐样品表格，所以列为高优先级待取全文案例，尚未加入定量比较。[期刊记录](https://doi.org/10.1016/0014-3057(92)90087-I)
- 扩展文献审计又转录了 S/MA 的逐样本匹配数据，与 SAN 案例构成一个有信息量的正负对照。SAN 微乳液中 RR 推 heterodyad 对 NMR 的 MAE 为 0.082；五个 DSC Tg 样本内留一时，Fox MAE 3.22 °C，RR-dyad 修正 1.07 °C。S/MA 本体聚合样本中，按报告 feed 输入的 terminal RR 对 NMR heterodyad 的 MAE 为 0.266，Tg 留一 MAE 3.70 K，高于 Fox 的 2.78 K。S/MA 只有 3 个实测 DSC Tg。
- 若用 NMR 推得的目标样品聚合物组成反解进料，S/MA 的 RR dyad Tg 留一 MAE 降至 1.36 K；但此分析使用目标组成、每折仅两点拟合交叉项，不能当作独立预测增益。它提示组成/聚合位点条件与 RR 必须共同恢复真实序列，才能解释差异。
- 因此更新后的主线是：**什么时候 terminal RR 和其条件输入足够准确，以致序列统计能改善 Fox；当序列重建失配时，这种修正何时会变差？**已有一个正案例和一个失效案例，但都来自单一单体对/来源，样本内交互项拟合也很小，尚没有跨体系泛化证据。
- 条件审计表明此前报告的 0.32 heterodyad 宽度混合了 RR 来源和 Tg 样品组成。固定 `f=0.5` 的跨文献来源宽度中位数为 0.155；只改变报告组成时宽度中位数为 0.278。温度、溶剂和拟合方法元数据不完整；同 pair、同引用来源内的条件变化非常稀疏，无法可靠分摊来源间离散。
- 严格留一 Tg 来源并剔除跨源重复单体对后，tgvalue 测试上 Fox MAE 7.66 K，RF 组成+序列模型 12.61 K；styrenic 测试只有 2 个独立单体对，序列相对 Fox 仅改善 0.04 K。现有数据不支持稳定的来源外序列增益。

因此，P4 仍能发展成一篇聚焦明确的论文，但当前证据更适合定位为**范围与失效条件的系统评估**，而不是宣称已证明普适预测提升。核心需要补足更多独立来源的逐样本匹配表，并采用同时留出来源和单体对的验证。文献案例提取和来源/条件分析见 `sequence_tg_literature_case_audit.md`、`P4_supplementary_experiments.md`；S/MA 数据与复算脚本为 `styrene_methyl_acrylate_matched_case.csv`、`src/p4_literature_SMA_case.py`。

### E/V SI 核对与序列留出更新（2026-10-06）

已读取本地 E/V 补充信息并录入 Table S1。实验统一在 80 °C、乙烯 80 psi 和
MAO/BHT=1:0.4 下进行；SI 分别列出 3 个样品的催化剂量、VCH 体积、时间和产率。
用 Sample 2 的 C2 RR (`rE=14,rV=0.14`) 并固定 Sample 3 正文报告的聚合物组成
`xV=18.1 mol%`，预测 Sample 3 NMR 序列，dyad MAE/RMSE 为 0.159/0.176 个百分点，
triad MAE/RMSE 为 0.725/0.790 个百分点。这是条件于已知最终组成的跨样品序列预测，
不是只凭初始进料的全盲预测。SI 未报告总液体/溶剂体积或乙烯 uptake 时间历程，
故无法从初始 VCH 和反应器设定重建半间歇聚合中乙烯的有效浓度轨迹。RR 也由同一研究
的 NMR 序列拟合得到，验证仍须标注这一依赖。条件数据、逐项误差和代码见
`ev_2026_s1_conditions.csv`、`ev_2026_c2_sample2_to_3_sequence_holdout.csv`、
`p4_ev_c2_sequence_holdout.json` 和 `src/p4_ev_c2_sequence_holdout.py`。

### 两阶段证据的新整理

当前证据需分开回答“序列能否重建”和“序列项是否改善 Tg”。选取的六个 NMR 检查中，
RR→heterodyad MAE 为 0.043；但六例不构成独立 Tg 验证集。SAN 的序列误差较小且小样本
留一 Tg 好于 Fox；S/MA 按 feed 推序列时误差较大且 Tg 低于 Fox。E/V 条件组成下预测
序列很准，但没有建立外部 Tg 增益。直接 NMR 序列与 Tg 的 PLGA 2024 结果在组成内有
统计关联，不过不测试乙烯基 RR 链路。

新补入的文献候选包括 S/EA 1993、AN/EA 1992 和 MMA/EA 2001。尤其 MMA/EA 2001 的
摘要称，转化率升高时 NMR dyad 比 RR-derived dyad 更精确地改善 Tg 计算，并报告加入
dyad tacticity 可再略微改善。这是“聚合历史/序列表征质量会控制热学预测价值”的直接
线索，但需取得逐样本表格后才能作定量证据。所有案例状态与评价边界集中在
`p4_two_stage_evidence_matrix.csv` 和 `sequence_tg_literature_case_audit.md`。

### 新增分析：序列效应的二维杠杆图

把 RR/报告序列相对同组成 Bernoulli 链的异二元组变化 `ΔH`，与体系异二元组 Tg
相对均聚二元组均值的偏移 `ΔT_AB` 相乘，得到 Barton 近似下的序列位移：
`ΔTg_seq = ΔH × ΔT_AB`。在 Johnston 的 31 行、5 个体系诊断数据中，样本级计算
序列位移范围为 −9.5 至 +4.5 K。aMS–AN 的中位绝对位移约 6.0 K；MMA–VCl 尽管
`|ΔT_AB|=41.4 K` 最大，但 RR 序列只轻度偏离 Bernoulli，位移中位数约 2.4 K。
S–AN 的交互项为正，其余四个已有交互参数的体系为负。

这支持将主线表述为**二维适用域**：RR 必须产生足够大的序列偏离，且该体系的 dyad–Tg
交互必须足够强，序列修正才可能明显；交互项的正负决定方向。当前计算使用来源级
`Tg_AB` 和 RR，且 Fox 残差远大于计算出的序列位移，故它是机理诊断而不是独立的 Tg
预测证据。逐样本数据与限制见 `p4_sequence_regime_map.csv/json` 和补充实验 S10。

### Novelty 压力测试后的主线收窄

文献核查确认，Barton 1970 已使用反应率推导 dyad 分布并报告 11 个体系的 Tg 一致性；
SAN 2000、S/EA 1993、AN/EA 1992 也明确结合 RR/NMR 序列与 Tg。MMA/EA 2001 还报告
高转化时 NMR dyad 比经典 RR dyad 更有利于 Tg 计算。因此，“序列影响 Tg”或“RR 序列可用于
Tg 修正”不宜作为新颖性主张。

可辩护的论文主线现在是**两阶段预测链路的误差传播与适用边界**：先问 RR 加上实际可得
的进料/组成信息能否预测 NMR dyad；再问这一序列误差和热学 dyad 杠杆共同决定的修正，
能否在 leave-source/leave-pair 外推中稳定超过 Fox。当前矩阵显示 SAN 为小样本正例、
S/MA 为报告 feed 下反例、E/V 为固定最终组成下的条件序列留出；尚无无目标信息泄漏的
跨单体端到端成功案例。这个更窄的定位及最低数据门槛见补充实验 S12。

---

## 5. 局限与下一步

- **强序列样本仍然少**：交替 vs 无规的直接对照主要依赖 Johnston 5 体系与少量
  Matched Groups（如 S-MMA +21 K、aMS-AN −3 K）。
- **Johnston Tg 指标存在体系参数复用**：`Tg_AB` 既是 Barton 模型输入，也来自同一体系资料；
  目前的总体误差不能解释为留出体系预测性能。详见 `P4_supplementary_experiments.md`。
- **组成 vs 投料**：`DeltaLearningCopolymerTg` 的组成列是共聚物组成（非投料），
  需显式区分 feed vs copolymer 组成，否则序列描述符的 `f1` 有系统误差。
- **MMA-VCl 的边界**：单一 Tg_AB 的二元组可加性在「Tg_AB 严重偏离均值 + 组成偏
  一端的样品」上失效，可引入三元组（Barton/Schneider）或连续链模型修正。
- **下一步**：① 收集更多强交替惰性共聚单体体系（St/fumaronitrile、St/maleimide）；
  ② 把 Matched Groups（S-MMA 嵌段/无规/交替三态）纳入 Barton/三元组模型；
  ③ 用 RRMaster 的大规模 RR 预测 ΔTg，做「哪些新单体能通过序列设计最大化 ΔTg」的
  正向筛选。

补充计算分析（PLGA 组成内序列关联、进料组成敏感性、RR 跨来源离散度传播）及复现脚本见
`projects/P4_sequence_tg/P4_supplementary_experiments.md`。

---

## 6. 复现

```bash
python3 src/p4_sequence.py                 # 引擎 + MC 校验 + 全对序列描述符
python3 src/p4_sequence_tg_validation.py   # 实验 2/3/4（NMR 二元组 + Barton + 符号规则）
python3 src/p4_real_sequence_tg.py         # 实验 5（近无规真实共聚物 Δ-learning）
python3 src/p4_tg_model.py                 # 实验 6 均聚物 Tg 模型
python3 src/p4_delta_learning.py           # 实验 6 Δ-learning 自检
```

输出文件位于 `projects/P4_sequence_tg/`：
`rr_sequence_dyad_validation.csv`、`rr_sequence_tg_validation.csv`、
`barton_sign_rule.csv`、`p4_sequence_tg_validation.json`、
`real_copolymer_tg_dataset.parquet`、`p4_real_sequence_tg.json`、
`p4_tg_model_results.json`、`p4_delta_learning.json`。
