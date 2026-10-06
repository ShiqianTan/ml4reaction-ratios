# P4 补充分析记录

## 审计 S0：Johnston Tg 指标的样本级拆分

对现有 `rr_sequence_tg_validation.csv` 按 `seq_type` 拆分原始绝对误差：31 条中 2 条标为 alternating、29 条标为 random。两条 alternating 的模型预测都直接等于输入的 `Tg_AB`，而该数据表中观测 Tg 也等于同一个 `Tg_AB`，所以这两条误差为 0 K。其余 29 条 random 样本仍使用研究体系的 `Tg_AB` 参数；这不是对新体系的留出验证。

| 子集 | n | Fox MAE (K) | dyad 模型 MAE (K) |
|---|---:|---:|---:|
| 全部 | 31 | 15.28 | 7.84 |
| 仅 random | 29 | 14.98 | 8.38 |
| alternating | 2 | 19.65 | 0.00 |

去掉直接零误差的两行后，样本级 MAE 改善仍存在，但因为 `Tg_AB` 由同一体系数据提供，仍不能视为独立预测。后续比较必须按体系留出，并只用训练体系估计序列 Tg 参数。

## 实验 S1：固定 PLGA 组成内的序列指标–Tg 关联

### 目的

检查一组不属于 Johnston 拟合集的公开实验数据，判断在固定单体组成后，实测序列指标是否仍与 Tg 有关联。数据来自 Patil 与 Won 的 PLGA 研究；该工作通过改变进料速率改变序列均匀性，并报告了 `13C NMR` 序列统计和 DSC Tg（*Macromolecules* 2024, DOI: [10.1021/acs.macromol.4c00106](https://pubs.acs.org/doi/10.1021/acs.macromol.4c00106)）。

### 分析

- 12 个样品，4 个固定 LA 质量分数组（0.4、0.5、0.6、0.7），每组 3 个进料速率。
- 对每个变量先减去所属组成组均值，再计算与 Tg 的 Pearson 相关；以组内置换 Tg 标签进行双侧精确置换检验（每组 3! 种，共 1,296 种联合排列）。
- 同时报告 6 项探索指标的 Holm 多重比较校正。此检验用于发现关联，不是因果识别。

### 结果

| 变量 | 组成内相关 r | 原始精确 p | Holm 校正 p |
|---|---:|---:|---:|
| alternation | −0.829 | 0.0069 | 0.0417 |
| blockiness | +0.640 | 0.0224 | 0.0895 |
| LL | +0.824 | 0.0069 | 0.0417 |
| LG | +0.641 | 0.0224 | 0.1120 |
| feed rate | +0.755 | 0.0285 | 0.0895 |
| dispersity | −0.035 | 0.9259 | 1.0000 |

在四个组成组中，最低进料速率相对最高进料速率的 Tg 都较低，差值为 −1.84 至 −4.36 °C；同时 alternation 上升，blockiness 下降。组成内 alternation 与 Tg 的负相关在 Holm 校正后仍显著。该结果与“序列分布携带组成之外的信息”相符，但不能把 Tg 差异单独归因于序列：进料速率是共同变化的工艺变量，且当前摘录没有绝对分子量，只有 dispersity。

### 解释边界

1. 这是一个单体对、单篇研究的 12 个样本；不能证明跨化学体系可泛化。
2. 该 PLGA 体系由开环聚合制备，不能直接验证 P4 的 Mayo–Lewis 自由基共聚链路。
3. 该研究本身已经提出 PLGA 序列分布与 Tg 的关系，因此本分析是独立数据上的再分析和稳健性检查，不是新颖性主张。
4. 显著的组内关联仍可能受分子量、残余单体、热历史及其他随进料速率改变的变量影响。

### 当前结论

可把这组数据作为跨聚合机理的外部一致性检查，支持将论文问题收窄为“在组成固定时，实测序列统计是否与 Tg 残差有关”。它不能替代对 RR 推 dyad 的验证，也不能弥补 Johnston 结果的参数泄漏问题。

### 复现

```bash
python3 src/p4_supplementary_plga.py
```

输入、脚本和结果：

- `projects/P4_sequence_tg/plga_sequence_tg_input.csv`
- `src/p4_supplementary_plga.py`
- `projects/P4_sequence_tg/p4_supplementary_plga.json`
- `projects/P4_sequence_tg/p4_supplementary_plga_tests.csv`
- `projects/P4_sequence_tg/p4_supplementary_plga_contrasts.csv`

## 实验 S2：进料组成假设与单体对重采样

### 目的与方法

对 210 个真实共聚物 Tg 测量、40 个单体对做同一组留一单体对交叉验证。保留实测共聚物摩尔组成作为组成特征；序列特征分别按两种输入计算：

1. 现有近似：把实测共聚物组成当作进料组成。
2. 敏感性版本：对每个样品用 Mayo–Lewis 方程反解一个能产生实测共聚物组成的瞬时进料组成，再据此计算序列描述符。

比较 Fox 重量分数基线、组成/反应率模型、组成/序列模型。对每个预测模型与 Fox 的 MAE 差，按单体对做 5,000 次配对 cluster bootstrap。

### 结果

- 反解进料与实测共聚物摩尔分数的绝对差中位数为 **0.089**，因此把共聚物组成当进料组成确实改变了不少样品的序列输入。
- Fox 重量分数基线 MAE 为 **9.28 K**；组成 + 序列的 RF MAE 为 **10.36 K**（旧近似）和 **10.17 K**（反解进料）。
- 组成 + 序列的 Ridge MAE 从 **10.33 K** 降到 **10.08 K**，仍高于 Fox 基线。
- RF 序列模型相对 Fox 的 MAE 差，其单体对 bootstrap 95% 区间为 **−1.12 至 +3.81 K**（旧近似）和 **−0.95 至 +3.07 K**（反解进料），都跨过零。

这说明进料/聚合物组成不一致是重要的输入不确定性，但在这个数据集上用单个反解进料替代后，并没有得到稳健的序列预测增益。该反解假设低转化或单一瞬时状态；对存在明显组成漂移的样品，真正的序列统计需要随转化积分。

### 解释边界

- 当前 40 对样品主要是近无规体系，强非理想序列样本少；阴性结果不排除强序列体系中的 Tg 效应。
- 留一单体对仍可能让同一个单体同时出现在训练与测试集。
- Bootstrap 区间描述这批单体对的抽样波动，不是实验测量误差或跨化学空间置信区间。

### 复现

```bash
python3 src/p4_supplementary_feed_sensitivity.py
```

产物：`p4_supplementary_feed_sensitivity.json`、`p4_supplementary_feed_sensitivity.csv` 和带反解进料特征的 parquet，均位于 `projects/P4_sequence_tg/`。

## 实验 S3：RR 跨来源离散度向 dyad 预测的传播

### 方法

对 40 个 Tg 单体对，从 RRMaster 取正值、标记为 primary 的记录。旧版计算把来源 RR 中位数和每个样品不同的组成一起放入同一分布。它计算正确，但不能被解释成“来源造成的 dyad 不确定性”。新版固定 `f=0.5` 计算来源间 RR 散布，并另用每个 pair 的中位 RR 计算样品组成变化带来的 dyad 散布；来源 ID 优先用文献题名/出处，不能把数据库导入桶当作原始文献来源。

### 结果

- 40/40 对都有 RR 记录；其中 378/495 条 primary RR 记录有可辨认的引用来源，共 185 个来源标签。24/40 个 pair 有多个可辨认来源。
- 固定等摩尔组成后，跨文献来源的 heterodyad 90% 宽度中位数为 **0.155**，四分位数为 **0.262**。
- 固定 pair-median RR、只改变 Tg 样品的报告组成时，heterodyad 90% 宽度中位数为 **0.278**，四分位数为 **0.365**。
- 旧版混合宽度复现为中位数 **0.320**、四分位数 **0.443**。这是 RR 来源与样品组成共同造成的宽度，不能称为纯来源差异，也不能将前两项简单相加去“解释”它。

这说明原先的 **0.32 主要不是一个干净的 RR 来源不确定性统计量**；组成变化本身就足以形成接近此量级的 dyad 分散。跨来源差异仍存在，但大小依赖可识别引用来源的覆盖。

### 限制

来源间离散度仍混合温度、溶剂、聚合方式、拟合方法和数据质量差异，**不能称为统计置信区间或纯测量误差**。

### 复现

```bash
python3 src/p4_supplementary_rr_spread.py
```

输出：`p4_supplementary_rr_spread.json` 与 `p4_supplementary_rr_spread.csv`。

## 实验 S4：RR 来源/条件分层与来源留出

### RR 条件分层

在同一单体对内，以固定 `f=0.5` 的预测 heterodyad 为响应，先减去 pair 均值，再比较重复条件单元的组均值差异。单独分层时，温度区间、溶剂族、拟合类别、聚合方法分别解释约 **7.9%、5.1%、12.6%、6.5%** 的可分析组内平方和；进一步限定在相同 pair、相同引用来源内，解释量降至 **0%、0.15%、0.01%、3.9%**。比例是条件类别的描述性组内解释量，不是因果份额，彼此不可相加。

条件元数据完整度有限：温度 71%、溶剂 75%、拟合方法 54%、聚合方法 64%；来源级控制下许多 pair-文献-条件组合没有重复值。当前数据因此不能可靠回答 0.155 的来源间宽度有多少由某一个条件造成。更严谨的结果是：本批数据没有足够的同文献同 pair 条件对照来拆分这些贡献。

### 留一来源预测

当前 Tg 数据只有 `tgvalue` 和 `styrenic` 两个来源，且 styrenic 的 7 个 pair 中有 5 个在 tgvalue 中重复。故同时报告原始留源结果和剔除跨来源重复 pair 的严格结果：

| 测试来源（训练另一来源） | 严格测试样本 | Fox MAE | RF 组成校正 MAE | RF 组成+序列 MAE |
|---|---:|---:|---:|---:|
| tgvalue | 137 行 / 33 对 | 7.66 K | 12.13 K | 12.61 K |
| styrenic | 18 行 / 2 对 | 23.57 K | 23.79 K | 23.53 K |

序列特征在第一折使 RF 组成模型误差增加 0.48 K；第二折只改善 0.26 K，且测试仅 2 对。因此未观察到可重复的来源外推增益。此前留一单体对的 RF 组成+序列 MAE 为 10.36 K，Fox 为 9.28 K；配对单体对 bootstrap 区间跨零。当前数据对跨来源泛化的支持不足。

复现：`PYTHONPATH=src python3 src/p4_condition_and_source_cv.py`。结果位于 `p4_condition_source_cv.json`。

## 实验 S5：文献中一组可逐行核验的 RR—NMR dyad—Tg 案例

对 Sanghvi 等 2000 年的苯乙烯–丙烯腈微乳液数据录入 6 行；其中 5 行报告 DSC Tg。论文同时给出系统级 `rS=1.49, rA=0.022`、13C-NMR triad/diad 统计和 Tg，且报告 Barton 二元组 Tg 与 DSC 的一致性。论文在 70 °C、苯乙烯/丙烯腈微乳液条件下进行实验。[原文及表格](https://electronicsandbooks.com/edt/manual/Magazine/E/European%20Polymer%20Journal/2000%20Volume%2036/2000%20Volume%2036/10/2275-2283.pdf)。

把文中聚合位点的 `M0S/M0A` 转为进料苯乙烯摩尔分数，以论文 RR 计算 Mayo–Lewis dyad；其 heterodyad 对 NMR 值的 MAE 为 **0.082**（5 个有 Tg 样本中为 0.077）。用五个 Tg 点做逐点留一、从其余四点拟合 dyad–Tg 交叉项：Fox MAE **3.22 °C**，NMR dyad MAE **0.78 °C**，RR 推 dyad MAE **1.07 °C**。这是一条“RR 推序列可能超过 Fox”的正向先例。

但 n=5、单一单体对/文献；RR 由同一聚合系列估计，留出的是 Tg 而不是 RR；每折也用其余四个 Tg 点拟合 dyad–Tg 系数。因此只能称为**单体系内部留一案例**，不是独立化学体系的端到端证据。更重要的是，NMR dyad 与 Barton Tg 的关联早已有文献证明，P4 的新颖性必须收窄到“RR 序列预测何时能转移并稳定改善 Fox”。

Brar 与 Sunita 关于苯乙烯–乙基丙烯酸酯（1993）和丙烯腈–丙烯酸甲酯（1992）的文章也报告 RR、13C-NMR 序列和 Tg 三类结果；目前还没有取得可复核的完整逐样本表。[论文清单和提取状态](sequence_tg_literature_case_audit.md)。

输入、脚本和结果：`san_microemulsion_matched_case.csv`、`src/p4_literature_SAN_case.py`、`san_microemulsion_matched_case_results.csv`、`san_microemulsion_tg_loo_predictions.csv`、`p4_san_matched_case.json`。

## 实验 S6：S/MA 文献案例复核——RR 序列何时会低于 Fox

### 数据与方法

从 Brar 与 Sunita 的苯乙烯–丙烯酸甲酯文章 Table 3–4 转录 6 个样品的
13C-NMR dyads、Kelen–Tudos `rS=1.20, rM=0.18` 和 Barton Tg；只有 3 个样品有
DSC Tg。条件为 70 °C、本体聚合、BPO 0.5 wt%、5–10% 转化率。[原始论文及表格](https://nopr.niscpr.res.in/bitstream/123456789/45983/1/IJCA%2031A%2812%29%20903-908.pdf)

对同一组 RR 做了两种计算：①直接使用表中 feed 摩尔分数；②由 NMR dyads 得到
聚合物苯乙烯摩尔分数 `xS=FSS+FSA/2`，再反解 Mayo–Lewis 进料，使 RR 预测的瞬时
共聚物组成与该值匹配。第二种只是隔离序列的敏感性分析，不是独立测量的 feed。
Tg dyad 模型固定均聚物端点 100 °C 和 6 °C，每次留一个 DSC Tg，仅用其他两点拟合
一个 `Tg_SM` 交叉项。

### 结果

| 分析 | n | 结果 |
|---|---:|---:|
| feed RR heterodyad 对 NMR | 6 | MAE 0.266 |
| 组成约束 RR heterodyad 对 NMR | 6 | MAE 0.233 |
| Fox（NMR 推得组成）Tg | 3 | MAE 2.78 K |
| NMR dyad Tg 模型，样品留一 | 3 | MAE 1.37 K |
| feed RR dyad Tg 模型，样品留一 | 3 | MAE 3.70 K |
| 组成约束 RR dyad Tg 模型，样品留一 | 3 | MAE 1.36 K |

按报告 feed 计算时，terminal RR 明显高估 heterodyad，且 Tg 误差高于 Fox。用目标样品
的 NMR 组成反解 feed 后，dyad-Tg 误差降至 Fox 以下，但该步骤使用了目标样品的实测
组成，不能视为独立 RR 预测。每折又只有两个 Tg 点拟合交叉项，`n=3` 的 MAE 仅作
诊断，不是可靠性能估计。

### 对主线的含义

与 S5 的 SAN 正向案例合看，现有证据指向可检验的条件命题：RR 序列只有在 terminal
模型及其聚合位点组成输入能恢复样品真实序列时，才可能改善 Fox；若条件不匹配，序列
修正会使 Tg 预测退化。最值得投入的下一项工作是再找到可逐行核验的独立来源匹配案例，
然后做留一来源和留一单体对验证。当前两个案例一个正、一个负，不能声称存在跨体系增益。

复现命令：`PYTHONPATH=src python3 src/p4_literature_SMA_case.py`。数据和输出见
`styrene_methyl_acrylate_matched_case.csv`、`styrene_methyl_acrylate_matched_case_results.csv`、
`styrene_methyl_acrylate_tg_loo_predictions.csv`、`p4_sma_matched_case.json`。

## 文献筛查 S7：2026 年 E/V 序列–热性质研究

检索到 Cipullo 等 2026 年开放获取的乙烯–乙烯基环己烷研究。文章以天然丰度
2D 13C–13C INADEQUATE 完成谱峰归属，给出 3 个样品的 NMR dyad/triad 分布，以
定量序列统计确定两类催化剂体系的反应率；作者报告组成和分子量相近、序列分布
不同的样品出现约 **10 °C 的 Tg 差异**。[文章 DOI](https://doi.org/10.1016/j.polymer.2026.130510)

Table 2 中样品 1 与 3 的 dyad 可直接核对：按 `xV=[VV]+[EV]/2` 计算，使用表中
已四舍五入的 dyad 得到 V 含量约为 18.05 与 18.25 mol%；原文报告样品 1 的 13C NMR
组成是 17.9±0.1%，样品 3 正文报告为 18.1 mol%。这些差异可能来自表格舍入/独立定量，
不应把 dyad 反算值当作更精确组成。两者 `[EV]` 为 35.1% 与 26.5%，`[VV]` 为 0.5%
与 5.0%。原文报告样品 1/3 的 Tg 约为 −23/−14 °C，Mn 约为 268/401 kDa；差异约
9 °C，与作者概述的约 10 °C 一致。组成很接近，但分子量仅属可比而非相同，所以不能
将全部 Tg 差异单独归因于序列。

Table 4 的样品级 RR/序列参数也已抄录：Sample 1/C1 `rE=19, rV=0.003, nE=4.68,
nV=1.03, B=1.2`；Sample 2/C2 `14, 0.14, 2.64, 2.12, 0.8`；Sample 3/C2
`16, 0.11, 6.17, 1.37, 0.9`。论文称 Sample 2 与 3 独立导出的 C2 反应率在不确定度内
一致。现已提取 SI Table S1，并完成用 Sample 2 的 RR 预测 Sample 3 序列的条件组成留出；
完整逐样品数据另存于
[`ev_2026_matched_case.csv`](ev_2026_matched_case.csv)。

它是重要的近期近邻工作，直接排除了“首次证明序列改变 Tg”作为 P4 的论文主张。
不过该研究的 RR 来自同组样品的 NMR 序列统计，并非预先给定、独立测得后再预测
留出样品序列。SI 条件已转录，但没有总液体/溶剂体积或乙烯吸收历程。因此暂列为
**近邻先例及条件组成下的跨样品序列检查，而非进料盲预测案例**。应将“拟合 RR”与
“应用 RR”分开，避免同一批 NMR 数据既估计 RR 又验证 RR。

## 文献筛查 S8：AN–EA 1992 原文数据可得性

Brar 与 Sunita 的 AN–EA 论文摘要确认：本体聚合、BPO 引发；用 nonlinear EVM
估计反应率；由 13C NMR 三元序列计算 dyad 和平均序列长度；并称 NMR 三元分布与
RR 计算值吻合、DSC Tg 与 dyad 估计吻合。[期刊记录](https://doi.org/10.1016/0014-3057(92)90087-I)

这是目前最像第三个匹配案例的旧文献，但检索到的出版商页面只提供摘要，原文 PDF
受限，未找到可复核的逐样品组成、dyad、RR 和 Tg 表格。因此只列为待获取原文的高优先级
候选，暂不录入分析数据或当作验证结果。后续应优先通过图书馆/作者版本取得原文，再检查
是否能以不泄漏目标样品序列的方式重建 RR→dyad 预测。

## 实验 S9：E/V 补充信息转录与 Sample 2→3 序列留出

从本地 SI DOCX 的 Table S1 转录温度、乙烯压力、催化剂量、VCH 体积、反应时间和产率，
整理在 `ev_2026_s1_conditions.csv`。三组条件均为 80 °C、乙烯 80 psi、MAO/BHT=1:0.4；
Sample 1/2/3 分别使用 1.00/2.00/0.75 mL VCH，产率为 149/100/82 mg。

将 Sample 2 的 C2 RR（`rE=14,rV=0.14`）用于 Sample 3，固定 Sample 3 报告的
共聚物组成 `xV=18.1 mol%`，通过 terminal Mayo–Lewis/Markov 关系预测序列。结果：

| 统计量 | dyad | triad |
|---|---:|---:|
| MAE（百分点） | 0.159 | 0.725 |
| RMSE（百分点） | 0.176 | 0.790 |

EV 预测 26.32%、实测 26.5%；VV 预测 4.94%、实测 5.0%。该预测没有用 Sample 3 NMR
重拟合 RR，是目前清晰的跨样品序列检验；但它把报告的目标聚合物组成作为输入，且 RR
本身由同体系 NMR 序列推导，因此不是完全独立的进料到序列盲预测。

Table S1 没有给出甲苯/总液体体积或乙烯吸收量随时间变化。产率与最终组成可以事后估计
已聚合单体量，但不能替代初始进料输入和半间歇过程中的乙烯浓度轨迹。故 SI 让案例条件
更完整，却没有解除端到端预测的关键数据缺口。后续如获得独立的 feed/uptake 记录，才适合
增加由初始进料驱动的全盲预测；当前保留条件组成下序列准确度的解释。

复现：`PYTHONPATH=src python3 src/p4_ev_c2_sequence_holdout.py`。输入、逐项误差和汇总见
`ev_2026_matched_case.csv`、`ev_2026_c2_sample2_to_3_sequence_holdout.csv`、
`p4_ev_c2_sequence_holdout.json`。

## 实验 S10：序列偏离与 dyad–Tg 交互的二维作用图

### 目的与分解

为把“RR 序列何时可能改善 Fox”与“序列效应为何因体系而异”连起来，在固定聚合物
组成下，把 Barton 二元组模型相对于 Bernoulli 无规链的差写成：

`ΔTg_seq = (H_seq − 2xA(1−xA)) × (Tg_AB − (Tg_AA+Tg_BB)/2)`

其中 `H_seq` 是 RR 或报告序列类型给出的总异二元组分数；`2xA(1−xA)` 是同组成
Bernoulli 序列的异二元组分数；第二项是异二元组的 Tg 交互偏移。该乘积分解说明：
只有序列明显偏离无规、且异二元组 Tg 与均聚二元组均值不同，才有较大的序列杠杆。
交互项为正/负决定修正方向。

### Johnston 31 样本的诊断结果

对 5 个体系做了 fixed-composition 分解：random 样本由系统 RR 反解一个使瞬时组成匹配
报告共聚物组成的 feed，alternating 样本用 `H_seq=1`。这是敏感性/机理诊断，不能当成
盲预测。

| 体系 | n | `Tg_AB−均聚均值` (K) | 中位数 `|ΔH|` | 中位数 `|ΔTg_seq|` (K) | 计算序列位移范围 (K) | Fox MAE (K) |
|---|---:|---:|---:|---:|---:|---:|
| aMS–AN | 6 | −19.0 | 0.315 | 5.98 | −9.50 至 −4.73 | 25.00 |
| AN–MMA | 5 | −26.0 | 0.148 | 3.85 | −5.11 至 −0.99 | 13.60 |
| S–AN | 7 | +9.0 | 0.262 | 2.36 | +0.56 至 +4.50 | 6.64 |
| MMA–VCl | 7 | −41.4 | 0.058 | 2.42 | −3.52 至 −0.13 | 20.35 |
| BMA–VCl | 6 | −35.9 | n/a | n/a | n/a | 11.13 |

MMA–VCl 的 dyad 交互绝对值最大，但序列预测偏离理想无规较小，故计算得到的序列杠杆
小于 aMS–AN。S–AN 的交互项为正，说明偏向交替可提高该 Barton 序列项；其余四个已知
体系交互为负。BMA–VCl 缺少可用 RR，本表没有从 Bernoulli fallback 推断“无序列效应”。

Fox MAE 明显大于上述序列位移范围，说明在 Johnston 表中，简单二元组序列校正并不能
解释全部 Fox 偏差。Fox 残差还包括组成加和规则、经验参数和其他化学效应；因此该图用于
界定可能的序列贡献大小与方向，而非把整个 Fox 残差都归于序列。

### 解释边界与含义

- `Tg_AB` 是来源级系统参数，和汇总的 Tg 样本同源；这里没有按体系外推，因此不是
  留出体系预测。
- RR 有来源级或数据库估计；random 行的 feed 是由报告共聚物组成反解的诊断值，不是
  实测投料，且无法处理高转化时的组成漂移。
- 31 行只有 5 个单体对，范围不可视为普适阈值或显著性检验。
- 该分解给出后续文献采集优先级：优先找 `|ΔH|` 与 `|Tg_AB−均聚均值|` 都大，且 RR、
  组成、NMR dyad、Tg 能对齐的样品；再检验它们是否真的有 leave-source / leave-pair
  的 Fox 误差改善。

复现：`PYTHONPATH=src python3 src/p4_sequence_regime_map.py`。逐样本分解和汇总分别见
`p4_sequence_regime_map.csv` 与 `p4_sequence_regime_map.json`。

## 实验 S11：两阶段证据矩阵——序列准确度与 Tg 增益分开判断

为避免把“RR dyad 对 NMR 准”误写成“RR dyad 改善 Tg”，将当前结果分成两个连续但独立
的检验：

1. **序列重建：**报告 RR/条件输入是否预测了目标样品的 NMR dyad/triad；必须标明是否
   用目标样品组成、目标 NMR 或同一系列拟合 RR。
2. **热学意义：**仅在有实测 DSC Tg 的样本上比较 Fox 与 RR-sequence Tg 模型，并说明
   RR、序列交互项和样品来源是否在测试折之外。

证据矩阵见 `p4_two_stage_evidence_matrix.csv`。当前可量化图景：

- 六个选取的 RR→NMR 序列检查合并 dyad MAE 为 0.055、heterodyad MAE 为 0.043；这批
  样品不是与 Tg 同行的独立预测基准。
- SAN 有较好的序列重建（heterodyad MAE 0.082），样品内 leave-one-out 的 RR dyad Tg
  MAE 1.07 °C，优于 Fox 3.22 °C；但 RR 和 Tg 交互参数没有按来源/单体对独立留出。
- S/MA 按报告 feed 推 sequence 的 heterodyad MAE 为 0.266，RR-dyad Tg LOO MAE
  3.70 K，高于 Fox 2.78 K。把 feed 反解到目标样品组成后，sequence MAE 仍为 0.233；
  Tg MAE 虽降到 1.36 K，但目标组成参与了输入，不能当独立预测。
- E/V 的 Sample 2→3 序列留出 dyad MAE 为 0.159 个百分点、triad MAE 为 0.725 个
  百分点，但固定了 Sample 3 最终组成。Sample 1/3 的 Tg 约差 9 °C，可作为热学关联，
  不能隔离分子量和催化剂/聚合条件的影响。
- PLGA 2024 提供了一个热学侧的直接支持案例：12 个样本分属 4 个固定组成组，实测
  13C NMR 序列描述符与 Tg 在组内有关联；交替度的精确置换检验原始 `p=0.0069`，多重
  比较 Holm 校正后 `p=0.0417`。但该体系是开环聚合，feed rate 同时改变序列及其它
  材料属性，不检验乙烯基 RR→序列链路。

这组证据目前支持一个必要条件、但不支持充分条件：RR 序列重建差时，Tg 序列修正可能
恶化；重建准确后，仍要证明序列–Tg 关系可在独立样品/来源上迁移。SAN 目前是正向
候选，S/MA 是失败边界，E/V 是较强的条件序列留出；还没有跨单体对、无目标信息泄漏的
端到端成功案例。

### 新找到的匹配文献候选

Brar 与 Sunita 的 styrene–ethyl acrylate 1993 论文摘要明确报告 EVM 反应率、NMR triad
序列与 Barton Tg，并称 RR 计算 triads 与观测 NMR 相符、Barton Tg 与 DSC 相符。它是
新的高优先级候选，但本轮没有取得可核验的逐样本表格，因此未将摘要判断量化。来源：
[primary article DOI](https://doi.org/10.1002/macp.1993.021940618)。AN–EA 1992 也仍是
匹配候选：[primary article DOI](https://doi.org/10.1016/0014-3057(92)90087-I)。

另找到 Nguyen 等 2001 年的 MMA–EA 研究。摘要报告在低转化样品中，Barton/Johnston Tg
计算使用 NMR dyad 概率比使用经典 RR 方法更精确；高转化时差异更明显，加入 dyad 级
tacticity 又略微改善 Tg。它直接指出“RR 序列误差如何传入 Tg”的条件效应，因此是重要
反例/机制候选；目前只有摘要层面的定量叙述，尚未录入逐样本数值。[primary article DOI](https://doi.org/10.1002/pi.694)。

下一阶段优先从这些论文取得逐行表格，按同一协议录入组成、feed/转化率、
RR 来源与拟合方法、NMR dyad/triad、DSC Tg 和分子量；若新增案例不足以做来源/单体对
留出，则应将论文定位为多案例的有界验证，而非泛化性能论文。

## S12：主线新颖性压力测试与可证伪边界

### 文献先例改变了什么

主线不能再表述成“首次检验 RR 推导序列能否解释 Tg”。Barton 1970 已提出用 RR 得到
二元组序列分布、代入 dyad Tg 模型，并报告 11 个共聚体系与文献数据相符；SAN 2000、
S/EA 1993、AN/EA 1992 也明确报告了 RR/NMR 序列与 DSC/Barton Tg 的比较。MMA/EA 2001
进一步指出，高转化样品中 NMR dyad 对 Tg 计算优于经典 RR dyad，加入 dyad tacticity
略有改善。这些先例意味着“序列影响 Tg”与“RR sequence 能做 Tg 修正”都不是新发现。

仍可成立的差异化问题是：**当 RR、进料/瞬时组成、转化率与序列模型来源都被明确审计，
RR 序列相对 NMR 序列损失多少 Tg 精度；这种损失是否随组成漂移、模型失配与 dyad 热学
杠杆系统变化；其 Fox 增益能否迁移到未见来源/单体对？** 这是有边界的验证基准，而非
新建一种 dyad Tg 方程。Barton 对 11 个历史系统的总体 agreement 不等于今日要求的
样本级盲测：发表摘要没有给统一 Fox 对照、序列重建误差与按来源/单体对留出的误差。

### 两个阶段的可证伪预期

| 阶段 | 预期 | 反证/边界信号 | 当前观测 |
|---|---|---|---|
| RR→NMR 序列 | 低转化、组成和 RR 适配同一聚合条件时，dyad 误差相对小 | 高转化、组成漂移、相分配/半间歇有效 feed 未知、penultimate 机制存在时误差升高 | SAN 较好；S/MA 按报告 feed 明显失败；E/V 只在固定终点组成时跨样品准确；MMA/EA 摘要报告高转化时 NMR 优势更明显 |
| NMR→Tg | 序列偏离 Bernoulli 足够大且 dyad Tg 交互明显时，测量序列有可观修正杠杆 | 序列接近 Bernoulli 或交互接近零时，正确重建序列也不会显著优于 Fox；Tg 增益也可能被分子量/立构规整性等混杂 | S10 计算范围 −9.5 至 +4.5 K；PLGA 组成内有关联但多属性共同变化；SAN 小样本 LOO 有利、S/MA 报告 feed 下有害 |
| 端到端外推 | RR 预测应在目标组成/feed 可预先获得、且 Tg 关系不重拟合目标测试组时仍胜 Fox | 任何一阶段失败或输入用到目标 NMR/目标组成，就不能算盲预测增益 | 目前没有跨单体对、无目标序列信息泄漏的完整成功案例 |

这里“低转化较好、高转化较差”目前是待验证的方向性预期，不是已由本项目多案例回归
证实的定律。MMA/EA 是唯一直接点出 conversion 的文献候选；必须取得各行数值并核对
其 RR 的估计样本与 Tg 样本是否相互独立，才能把摘要结论变成定量证据。

### 论文是否值得继续的最低数据门槛

1. 逐样本匹配至少覆盖 3 个独立来源/聚合条件和 3 个单体对；每行同时有 RR、RR 拟合
   方法/来源、feed 或可重建瞬时 feed、转化率、共聚物组成、NMR dyad/triad、Tg 和可用
   混杂量（至少 Mn、tacticity 或热历史其中之一）。
2. 先报告 RR→NMR 的 dyad/triad MAE 与误差区间；区分目标组成已知、反应进料已知、
   目标 NMR 参与 RR 拟合三种信息等级，不把 conditional prediction 称作 blind。
3. 再在相同样本上比较 Fox、直接 NMR-sequence、RR-sequence；用 leave-source 与
   leave-pair 两种外推拆分。只在样本内重新拟合 dyad Tg interaction 的结果留作诊断。
4. 把预测增益对照两个预先定义的热学杠杆：序列相对 Bernoulli 的偏离 `|ΔH|` 和
   dyad interaction `|ΔT_AB|`。这可以解释“序列重建准但 Tg 无增益”并非矛盾。
5. 若提取后仍只有 SAN、S/MA 两个可量化端到端案例，文章应改为严谨的 matched-case
   synthesis/方法学边界，不报告 leave-pair 泛化结论，也不把结果包装成普适 ML 模型。

本轮文献核查来源：Barton 1970 [primary DOI](https://doi.org/10.1002/polc.5070300161)、
MMA/EA 2001 [primary DOI](https://doi.org/10.1002/pi.694)、S/EA 1993
[primary DOI](https://doi.org/10.1002/macp.1993.021940618)、AN/EA 1992
[primary DOI](https://doi.org/10.1016/0014-3057(92)90087-I)。
