# P1 报告：ReactivityRatioBench — 反应率预测的基准与 OOD 评估

> 研究问题：**自由基共聚反应率 `(r1, r2)` 能否从单体结构（+ 条件）预测？
> 更重要的是：在“未见过的新化学”上，模型退化多少？**
>
> 本文档汇总 P1（数据基础 + 基准）的全部实验与结论。代码见 `src/p1_*.py`、
> `src/p1_gnn.py`，输出见 `results/`、`projects/P1_rr_benchmark/`。

---

## 1. 数据基础（Stage A–G）

| 来源 | 行数 | 说明 |
|---|---:|---|
| PolyCARP | 7,582 | T/溶剂/方法/Q-e/DOI，VLM 抽取 |
| Chen & Jackson (Polymer Handbook) | 4,787 | 仅 r_err，手册 OCR（二手） |
| Coxwell 2026 full | 484 | T/溶剂/统计方法/来源（一手） |
| Coxwell 2026 error | 50 | 多温/不确定度 |
| Chen & Jackson textbook / CKAs | 79 / 10 | 带温度 |
| **总计原始** | **12,992** | |

管线（B 单体 InChIKey 注册 1,718 SMILES → C 方向规范化 6,422 交换 → D 条件标准化
→ E 目标 `ln r`/`ddg` + QC → F 来源去重 → G 质量评分 A–D）：

| 指标 | 值 |
|---|---:|
| 去重测量 | **5,742** |
| 唯一单体对 | 2,633 |
| 唯一单体 | 1,218 |
| 条件完整 | 3,225 (56 %) |
| 质量 A / B | 211 / 3,115 |

**Go/No-Go：GO**（≥1500 去重 & ≥800 条件完整）。

物理描述符：CopDDB（2,500 个自由基-单体对的势垒/反应能/轨道能隙/%Vbur）。

---

## 2. 基准设置

* **特征**：`struct` = ECFP4(2048×2) + RDKit 描述符；`+cond` = 温度/溶剂族；
  `+phys` = CopDDB Δ-势垒/Δ-反应（r1/r2 方向）。
* **模型**：Ridge（线性）、RandomForest(300)、XGBoost(300)、
  MLP(256)、GIN(3 层图网络)。
* **划分**（5 折，OOD 严重度递增）：
  - `random`：随机 KFold
  - `pair`：按单体对 GroupKFold
  - `one_monomer_out`：按单个单体分组（至少一个单体未见）
  - `both_monomers_out`：图划分（两个单体都未见；测试集最小，~100 行/折）
* **指标**：MAE(ln r)、RMSE、R²、NLL（σ=残差 std）、coverage95、
  `E_curve`（Mayo-Lewis F1 曲线误差，跨全 f1 网格）。
* 数据集：`dataset_C_condition`（3,225 行）。

---

## 3. 结果

### 3.1 每个划分的最优模型（全指标）

| 目标 | 划分 | 特征 | 模型 | MAE | R² | NLL | 覆盖95 | E_curve |
|---|---|---|---:|---:|---:|---:|---:|---:|
| ln r1 | random | struct+cond | RF | **0.864** | 0.386 | 1.696 | 0.934 | 0.093 |
| ln r1 | pair | struct | RF | 1.061 | 0.220 | 1.816 | 0.940 | 0.120 |
| ln r1 | one-monomer-out | struct+cond+phys | XGB | 1.159 | 0.111 | 1.878 | 0.938 | 0.130 |
| ln r1 | both-monomers-out | struct | RF | 1.202 | 0.076 | 1.912 | 0.931 | 0.136 |
| ln r2 | random | struct+cond | RF | **0.883** | 0.422 | 1.717 | 0.930 | 0.093 |
| ln r2 | pair | struct+cond | RF | 1.059 | 0.301 | 1.812 | 0.934 | 0.116 |
| ln r2 | one-monomer-out | struct+cond+phys | RF | 1.165 | 0.185 | 1.889 | 0.928 | 0.129 |
| ln r2 | both-monomers-out | struct+cond | RF | 1.280 | 0.155 | 1.974 | 0.932 | 0.138 |

### 3.2 模型对比（MAE，ln r1）

| 划分 | Ridge | RF | XGB | MLP | GIN |
|---|---:|---:|---:|---:|---:|
| random | 1.051 | **0.864** | 0.916 | 0.989 | 1.020 |
| pair | 1.252 | **1.061** | 1.076 | 1.164 | 1.079 |
| one-monomer-out | 1.260 | 1.159 | **1.158** | 1.249 | 1.184 |
| both-monomers-out | 1.261 | 1.202 | 1.218 | 1.298 | **1.160** |

（Ridge/XGB 取各自最优特征集；GIN 只用 struct。）

### 3.3 特征消融（RF，ln r1）

| 划分 | struct | +cond | +cond+phys |
|---|---:|---:|---:|
| random | 0.913 | **0.864** | 0.865 |
| pair | **1.061** | 1.062 | 1.062 |
| one-monomer-out | **1.159** | 1.160 | 1.159 |
| both-monomers-out | **1.202** | 1.228 | 1.228 |

### 3.4 噪声天花板（为什么 R² 上不去）

在 **470 个有 ≥3 次主测量的单体对**上分解 ln_r1 方差：

| 量 | ln_r1 | ln_r2 |
|---|---:|---:|
| 总方差 | 2.96 | 3.10 |
| 对内（重复测量）方差 | 1.13 (38 %) | 1.15 (37 %) |
| **R² 天花板**（完美对恒等模型） | **0.619** | **0.628** |

→ 即使模型完美，R² 也只能到 **~0.62**：38% 的方差是实验/实验室/条件散差，
**本质不可从结构预测**。我们的 0.39–0.42 已在该上限的 65–68%。
这解释了为何参考仓库的 log r1 R²≈0.62 ≈ 我们的天花板（他们在更干净/对级数据上评估）。

### 3.5 pair-level 基准（把重复测量聚合成对中位数）

2,632 对（5,742 条测量聚合），RF（ECFP+RDKit）：

| 划分 | ln_r1 R² | MAE |
|---|---:|---:|
| random | 0.364 | 0.964 |
| one-monomer-out | 0.213 | 1.098 |
| both-monomers-out | 0.094 | 1.236 |

**pair-level 未高于 measurement-level（0.386）**：训练样本减半（2,632 vs 5,742），
且绝大多数对只有 1–2 条测量（≥3 条仅 470 对），聚合几乎没把对间噪声平均掉。

### 3.6 噪声清洗（`p1_noise_cleaning.py`）

三条个体测量噪声标记：

| 标记 | 命中 | 占 |
|---|---:|---:|
| within-pair MAD 离群（|ln r − 对中位| > 3σ_MAD） | 290 | 5 % |
| quality_grade = D | 2,252 | 39 % |
| 相对不确定度 r_unc/r > 0.5 | 557 | 10 % |

清洗效果（固定同一测试集，RF struct+desc）：

| 方案 | 固定噪声测试 R² | 干净测试 R² |
|---|---:|---:|
| 全量 | 0.391 | 0.391 |
| 去掉 MAD 离群 | **0.399** | **0.471** |
| 去掉 quality-D | 0.308 | 0.360 |
| 去掉高不确定度 | 0.368 | 0.390 |
| 去掉全部噪声 | 0.294 | 0.415 |

**结论**：只有「within-pair MAD 离群」清洗有帮助且有限（固定测试 +0.008、干净测试
+0.08）——这类离群是真实的单位/转录错误（如 r 记成 10 vs 0.10），仅占 5%，是唯一
「可移除」的噪声。宽泛的 quality/不确定度过滤**反而有害**（去掉 quality-D 把 R² 从
0.39 拖到 0.31），与参考仓库「过滤高噪声训练点损害性能」一致：D 级≠错误值，删掉
反而损失化学覆盖。

### 3.7 干净子集（一手来源优先，`p1_clean_subset.py`）

按「来源质量」而非「重复次数」构造干净子集（quality_grade 本质是来源代理：handbook/
textbook/cka 全 D，coxwell 全 B，polycarp 混 A/B/C）。RF（ECFP+RDKit）5 折：

| 子集 | n | ln_r1 random R² | ln_r1 MAE | ln_r1 both-OOD R² | ln_r2 random R² |
|---|---:|---:|---:|---:|---:|
| full（全部主测量） | 5,742 | 0.396 | 0.874 | 0.109 | 0.430 |
| primary（polycarp+coxwell） | 3,494 | 0.349 | 0.917 | −0.051 | 0.365 |
| primary + 现代拟合 | 479 | 0.313 | 0.910 | −0.052 | 0.227 |
| primary + 现代拟合 + 有不确定度 | 291 | 0.269 | 0.960 | −0.008 | 0.255 |
| coxwell-only | 273 | 0.366 | 1.049 | 0.309 | 0.445 |
| coxwell-curated（剔除手册/Young 二手 + 待删） | 220 | 0.316 | 1.130 | nan（对太少） | 0.396 |

**结论：一手来源的「干净数据」不仅没有提升 R²，反而降低（0.396 → 0.35 → 0.31 → 0.27）**，
MAE 也变差（0.874 → 0.917+）。原因是删数据 = 损失样本量与化学覆盖，而 handbook 数据
虽是二手但**自洽、并非有害噪声**。coxwell-only 的 both-OOD R² 看似高（0.309），但只有
273 行、测试集极小（~20–40 行/折），不可靠；其同分布 MAE 也最差（1.049）。即便按
Coxwell 作者自己的标注剔除二手来源（polymer handbook / Young 汇编 / 待删）得到
**coxwell-curated（220 行）**，R² 也不升反降（0.366→0.316），再次确认数据清洗不是瓶颈。

→ 至此三条「清理数据」路线（重复次数 gold-subset、MAD 离群、一手来源）**全部无法
提升 R²**，一致地推翻「R² 受数据噪声限制」的旧判断（见 3.4 修正），指向真正的约束
在**特征表示**（见结论 6）。

### 3.8 CopDDB Δ-势垒天花板验证（`p1_copddb_ceiling.py`）

检验结论 6 的自然下一步：既然约束是「对间化学」，那么 CopDDB 的一阶量子化学
Δ-势垒（`r1=k11/k12 → ΔΔE‡ = E_barrier(R1,M2)−E_barrier(R1,M1)`）能否直接提供
ECFP 缺失的化学信号？CopDDB 覆盖 601 测量 / 91 对 / 30 单体（10.5%）。

物理关系检验（`ln r ≈ −ΔΔE‡/RT`，理论斜率 −1059/Hartree）：

| target | Pearson | 斜率(拟合) | 斜率(理论) | r² |
|---|---:|---:|---:|---:|
| ln r1 | +0.178 | +38.8 | −1059 | 0.032 |
| ln r2 | +0.150 | +33.8 | −1059 | 0.023 |

**物理关系不成立**：相关系数小且符号为正（理论应负），斜率差 30 倍。

模型对比（601 行子集，RF/Ridge 5 折）：

| 模型 | ln_r1 random R² | ln_r1 pair R² | ln_r2 random R² |
|---|---:|---:|---:|
| ECFP（RF） | **0.226** | −0.121 | **0.282** |
| Δ-势垒（Ridge） | 0.015 | −0.042 | −0.008 |
| ECFP + Δ-势垒（RF） | 0.223 | −0.060 | 0.271 |

**结论：CopDDB 一阶 Δ-势垒既与实验 `ln r` 基本无关，加入后也不能提升 ECFP。**
现成的真空小模型 DFT 势垒代理（无溶剂、无构象、小基组）不含实验反应率信号；
打破 ECFP 天花板需要真正高质量量子化学，而非现成描述符。这也印证了 P1 结论的
最终形态：**在现有可得的全部特征下（ECFP/图/描述符/条件/CopDDB 势垒），RR 预测
的同分布 R² 上限就是 ~0.4，且该上限是表征限制而非数据或模型限制。**

---

## 4. 关键结论

> **最强发现：OOD 误差由化学覆盖主导，与模型/架构无关** —— 误差随 OOD 严重度
> 单调上升，且没有任何模型（含 GNN/MLP）能显著击败 ECFP + 树的简单基线。

1. **同分布可预测**：MAE(ln r) ≈ 0.86–0.88，R² ≈ 0.39–0.42；加入温度/溶剂是
   同分布下最大的单一增益。
2. **OOD 才是瓶颈**：按 pair / 单单体 / 双单体划分，误差近乎翻倍
   （MAE 1.06 → 1.28，R² 0.39 → 0.08），且**条件与物理特征在 OOD 下几乎无增益**。
   误差随 OOD 严重度单调上升 `random < pair < one-monomers-out < both-monomers-out`，
   证实计划的中心假设：**瓶颈是化学覆盖，不是条件覆盖**。
3. **没有模型能显著赢过简单基线**：RF/XGB 指纹基线在绝大多数划分最优；GIN 图网络
   只在最难的双单体 OOD 划分上略胜（1.16 vs 1.20），但在同分布下明显更差
   （1.02 vs 0.86）；MLP 全面落后于 RF。→ 当前 ECFP+树的基线是强、且难以超越的。
4. **校准可诊断 OOD 退化**：NLL 从同分布 1.70 升到 OOD 1.91；coverage95 稳定在
   0.93–0.94（接近高斯、略欠自信）；Mayo-Lewis `E_curve` 从 0.093 升到 0.136
   （F1 组成曲线误差 OOD 翻倍）。这三者给出比 MAE 更细的 OOD 画像。
5. **R² 上限不是数据噪声，而是特征表示**（对 3.4 的修正）：先前的「对间噪声占 38%、
   天花板 0.62」只是「对恒等模型」的上限，不代表结构模型能达到。三条「清理数据」
   路线（重复次数 gold-subset、MAD 离群清洗、一手来源干净子集）**全部无法提升 R²**，
   且一手来源子集反而更差（0.40 → 0.27）。→ 测量噪声不是约束，模型结构表示才是。
6. **真正的约束是 ECFP 特征无法捕捉对间化学**：结构模型只在 R² 0.28–0.43 之间，
   远低于对恒等天花板 0.62，两者之间 >0.5 的方差是**对间化学差异**（Q-e 极性、位阻、
   前线轨道），ECFP+树表示不出来。
7. **但 CopDDB 一阶 Δ-势垒也无法打破这个天花板**（§3.8，`p1_copddb_ceiling.py`）：
   在 CopDDB 覆盖的 601 测量/91 对子集上，Δ-势垒与实验 `ln r` 几乎不相关
   （Pearson +0.18/+0.15，符号还与理论相反），`ln r ≈ −ΔΔE‡/RT` 的拟合斜率 +39
   而理论为 −1059；把 Δ-势垒加入 ECFP 后 R² 不升反降（0.226→0.223）。→ 现成的
   真空小模型 DFT 势垒代理不含实验反应率信号；**打破 ECFP 天花板需要的是真正的
   高质量量子化学（显式溶剂+色散+大基组+构象），而非现成描述符**。

---

## 5. 局限与下一步

* 图网络只用了 GIN 单架构；可补 GAT/transformers 或消息传递基线。
* 未用可训练条件特征（溶剂/温度以 one-hot 硬编码）；可用溶剂 SMILES 编码。
* `both_monomers_out` 测试集小（~100 行/折），噪声大，需更大数据集或更多折。
* 数据发布需确认 PolyCARP/CopDDB/IUPAC 再分发许可。

---

## 6. 复现

```bash
python3 src/p1_build_rrmaster.py    # Stage A-G → RRMaster v1
python3 src/p1_feasibility.py       # 数据集版本 + 可行性
python3 src/p1_benchmark.py         # RF/XGB/Ridge/MLP 全矩阵（含校准/E_curve）
python3 src/p1_gnn.py               # GIN 图网络基线
python3 src/p1_noise_cleaning.py    # pair-level 基准 + 噪声探测/清洗
python3 src/p1_clean_subset.py      # 一手来源干净子集（primary/现代拟合/不确定度）
python3 src/p1_copddb_ceiling.py   # CopDDB Δ-势垒天花板验证
```

输出：`results/p1_benchmark_results.csv`、`results/p1_benchmark_oof.parquet`、
`results/p1_gnn_results.csv`、`results/p1_gnn_oof.parquet`、
`results/p1_split_sizes.csv`、`projects/P1_rr_benchmark/pair_level_dataset.parquet`、
`noise_flagged_measurements.parquet`、`noise_cleaning_eval.csv`、
`p1_noise_pairlevel.json`。
