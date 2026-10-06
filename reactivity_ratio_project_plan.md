# Reactivity Ratio 研究项目规划（可执行版：含数据来源与处理流程）

> 目标：围绕自由基共聚反应活性比（reactivity ratios, \(r_1,r_2\)）规划 **4 个可形成独立论文的项目**。优先选择：数据可获得、无需立即开展新实验、可共享数据基础设施、能回答明确科学问题的方向。
>
> 本版已删除此前可行性较低的 **Terminal / Penultimate / Depropagation Kinetic Model Selection** 项目，并将“数据从哪里来、如何获取、如何清洗、如何去重、如何形成训练/验证集”写入每个项目。

---

# 0. 项目组合与推荐顺序

| 优先级 | 项目 | 可行性 | 新颖性 | 数据准备难度 | 是否需新实验 | 推荐程度 |
|---|---|---:|---:|---:|---:|---:|
| **P1** | **ReactivityRatioBench：Condition / Quality / Uncertainty / Physics-aware RR Prediction** | ★★★★★ | ★★★★☆ | 中 | 否 | **★★★★★** |
| **P2** | **Temperature-Dependent RR Thermodynamics：\(\Delta\Delta H^\ddagger,\Delta\Delta S^\ddagger\)** | ★★★★☆ | ★★★★★ | 中高 | 否 | **★★★★☆** |
| **P3** | **ML Prior + Bayesian / Active Experimental Design for RR Estimation** | ★★★★☆ | ★★★★☆ | 中 | 否；可后续实验验证 | **★★★★☆** |
| **P4** | **RR → Sequence → Tg/Property：物理约束的多尺度 Δ-learning** | ★★★★☆ | ★★★★★ | 中高 | 否 | **★★★★☆** |

## 推荐执行顺序

```text
P1: 统一 RR 数据库 + Benchmark
        │
        ├──> P2: 温度依赖 kinetic thermodynamics
        │
        ├──> P3: ML prior + adaptive RR experiment design
        │
        └──> P4: RR → sequence → Tg/property
```

**P1 是所有项目的公共数据底座，必须优先做。**

---

# 1. 公共数据底座：所有项目共享的数据源

## 1.1 核心公开数据源

| 数据源 | 可获得内容 | 当前已知规模/特点 | 在本项目中的用途 | 链接 |
|---|---|---|---|---|
| **CoPolDB** | monomer、SMILES/InChI、copolymer pair、\(r_1,r_2\)、文献 DOI | 论文发表时：864 monomers、1954 copolymers、2991 copolymerization reactions；725 unique DOI 可对应 2449 条反应 | P1 主 RR 数据入口；P2 温度文献追溯入口 | https://www.copoldb.jp/ ; https://doi.org/10.1039/D3PY01372C |
| **PolyCARP / copolymer-reactivity** | VLM 文献抽取数据、reaction conditions、结构、RR、processed dataset、splits、model code | 3791 unique reactions、1206 publications、4969 measurement rows；包含 solvent、temperature、polymerization mechanism 等条件字段 | P1 condition-aware 数据；P2 多温度筛选；作为自动 PDF 抽取流程模板 | https://github.com/lamalab-org/copolymer-reactivity |
| **Chen & Jackson ReactivityRatios** | RR dataset、Polymer Handbook OCR extraction、DFT descriptors、barrier results、ML code | repo 含 `dataset/`、`dataset_gen/`、`descriptor_gen/`、`DFT_barrier_result/`；README 报告约 1500 experimental RR measurements | P1 benchmark/descriptor/数据血缘；Δ-learning baseline | https://github.com/TheJacksonLab/ReactivityRatios |
| **Coxwell et al. 2026 Structure–Reactivity** | >450 RR primary-literature dataset、SMILES、literature source、computed parameters、descriptors | SI ZIP 明确包含 RR dataset + SMILES + literature source + model descriptors | P1 高质量 primary-literature core；数据质量标注 | https://doi.org/10.1021/acs.macromol.6c00233 |
| **Coxwell et al. 2026 Better Models or Better Data** | 重复测量/数据噪声分析、模型外推分析 | 专门讨论 RR 预测性能与实验数据噪声上限 | P1 uncertainty/replicate 分析；quality-aware learning | https://doi.org/10.1021/acs.macromol.6c01785 |
| **CopDDB** | radical–monomer quantum descriptors、reaction energy、barrier、orbital/steric descriptors | 50 monomers × 50 radical–monomer combinations = 2500 pairs；24 类 reaction/electronic/geometric/conventional descriptors | P1 physics-informed descriptor；P2 辅助解释 | https://github.com/hatanaka-lab/CopDDB ; https://doi.org/10.1039/D4DD00266K |
| **Alexander et al. 2026 robotic dataset** | 标准化机器人 RR 测定；NMR/实验数据；6 个 monomer pairs | 6 组 IUPAC-compatible standardized RR systems | P3 最重要 retrospective active-design benchmark；P1 外部高质量验证集 | https://doi.org/10.1039/D6SC02232D |
| **IUPAC Reactivity-Ratios repository** | IUPAC \(f_0-X-F\) Python/Excel 代码、test set、Contour 兼容格式 | 统一输入为 \(f_0,X,F,\Delta F\)；含开放 Python 实现和测试集 | P1 quality checking；P3 likelihood/fitting/实验设计 | https://github.com/IUPAC-Polymer-Division/Reactivity-Ratios |
| **Polymer Handbook / Greenley RR tables** | 经典大规模自由基共聚 RR 汇编 | 多篇 ML 工作的数据祖先 | P1 legacy supplement；必须进行 provenance 去重 | 通过图书馆/订阅获取 |
| **你已下载的全文 PDF 文献集** | 原始表格、RR 条件、温度、溶剂、feed、conversion、fitting method、uncertainty | 当前最重要的“条件补全 + 数据质量标注”来源 | P1–P4 的人工/多模态补充数据 | 本地文献库 |

---

# 2. 公共数据处理流水线

所有项目不要分别清洗数据。先构建一个统一的 `RRMaster` 数据库。

## 2.1 建议目录

```text
reactivity-ratio-research/
│
├── data/
│   ├── raw/
│   │   ├── copoldb/
│   │   ├── polycarp/
│   │   ├── chen_jackson/
│   │   ├── coxwell/
│   │   ├── copddb/
│   │   ├── alexander_robotic/
│   │   ├── iupac/
│   │   └── literature_pdf_extracted/
│   │
│   ├── staging/
│   ├── standardized/
│   ├── provenance/
│   └── benchmark_splits/
│
├── metadata/
│   ├── references.csv
│   ├── monomer_registry.csv
│   ├── solvent_registry.csv
│   └── source_lineage.csv
│
└── projects/
    ├── P1_rr_benchmark/
    ├── P2_temperature_rr/
    ├── P3_active_rr_design/
    └── P4_sequence_tg/
```

---

## 2.2 Stage A — 原始数据入库

**原则：raw 数据永远不覆盖。**

每条记录至少增加：

```text
source_dataset
source_file
source_row_id
source_doi
source_reference
extraction_type
    = original_csv / SI_excel / database / OCR / VLM / manual
retrieval_date
```

输出：

```text
data/staging/all_raw_measurements.parquet
```

---

## 2.3 Stage B — 单体身份标准化

对所有来源中的 monomer：

1. 名称映射：trivial name / abbreviation / IUPAC name
2. SMILES 校验
3. RDKit sanitize
4. canonical SMILES
5. InChIKey
6. 去 counterion / solvent / salt
7. 保留原始 SMILES，不覆盖

建立：

```text
monomer_id
monomer_name_raw
monomer_name_canonical
smiles_raw
smiles_canonical
inchikey
source
manual_review_flag
```

### 特别注意

- MMA / methyl methacrylate 等简称必须统一。
- 同一 monomer 的不同 tautomer / protonation state 应人工确认。
- 对 vinyl monomer，双键反应位点必须一致。
- Radical representation 与 monomer representation 分开保存。

---

## 2.4 Stage C — RR 方向标准化

RR 有方向性：

\[
r_1=\frac{k_{11}}{k_{12}},\qquad
r_2=\frac{k_{22}}{k_{21}}
\]

所以不能简单把 monomer pair 排序后不处理 label。

### 推荐规则

创建 canonical monomer ordering：

```text
canonical_pair = sorted(InChIKey_1, InChIKey_2)
```

若原始顺序与 canonical order 相反，则同时交换：

```text
monomer_1 <-> monomer_2
r1 <-> r2
f1 <-> f2
F1 <-> F2
```

保留：

```text
original_orientation
canonical_orientation
orientation_swapped = True/False
```

---

## 2.5 Stage D — Reaction condition 标准化

统一字段：

```text
temperature_K
solvent_name
solvent_smiles
solvent_family
initiator
polymerization_method
polymerization_mechanism
pressure
reaction_time
```

### Temperature

- °C 全部转 K：

\[
T_K=T_{^\circ C}+273.15
\]

- 若文献给范围，则保留：

```text
temperature_min_K
temperature_max_K
temperature_nominal_K
```

### Solvent

建立 solvent registry：

```text
solvent_id
name
canonical_smiles
 dielectric_constant
ET30 (optional)
Hansen_deltaD/P/H (optional)
```

对于 bulk polymerization：

```text
solvent = BULK
```

不要填作 missing。

---

## 2.6 Stage E — RR 数值与派生 target

保留原值：

```text
r1_raw
r2_raw
```

经过确认后：

```text
r1
r2
ln_r1
ln_r2
log10_r1
log10_r2
```

如有温度，则计算：

\[
\Delta\Delta G_i^\ddagger=-RT\ln r_i
\]

对应：

```text
ddg1_kJ_mol
ddg2_kJ_mol
```

### QC

- \(r\le0\)：flag，不直接进入标准训练集。
- 极端值如 \(r>100\) 或 \(r<0.01\)：保留，但 `extreme_rr_flag=True`。
- OCR 提取的科学计数法必须与 PDF 原表二次核对。

---

## 2.7 Stage F — 数据血缘与重复去除

这是 P1 最重要的数据工程部分。

同一实验可能出现在：

```text
Primary paper
  ↓
Polymer Handbook
  ↓
CoPolDB
  ↓
Safari/other ML dataset
```

不能把它们当作 4 个独立样本。

### 去重键按优先级匹配

**Level 1：精确 primary DOI + monomer pair + RR**

```text
primary_doi
canonical_pair_id
r1
r2
```

**Level 2：DOI 缺失时使用文献信息**

```text
author + year + pair + temperature
```

**Level 3：近似重复**

```text
pair + round(r1,2) + round(r2,2) + temperature
```

Level 3 只做候选，不自动合并。

### 建议字段

```text
measurement_id
primary_measurement_id
duplicate_group_id
primary_source
secondary_sources
is_primary_measurement
```

**训练/测试 split 必须以 `primary_measurement_id` 或 `duplicate_group_id` 为 grouping unit。**

---

## 2.8 Stage G — Experimental quality score

依据 IUPAC / Lynd / van Herk 的建议建立质量标签。

### Quality A

- primary literature 可追溯
- \(f_0,X,F\) 信息可恢复
- nonlinear/EVM/IUPAC-compatible fitting
- uncertainty 或 confidence contour 可获得
- temperature / solvent 完整

### Quality B

- primary literature 可追溯
- 条件完整
- RR estimation 方法明确
- 缺 uncertainty 或完整原始点

### Quality C

- primary literature 可追溯
- 旧式 Fineman–Ross / Kelen–Tüdős 等线性化估计
- 条件部分缺失

### Quality D

- secondary source/handbook only
- experimental protocol 不可恢复

额外保存连续评分，例如：

```text
quality_score = 0–100
```

评分组件：

```text
primary_source_traceable   20
conditions_complete        15
f0_X_F_available           20
modern_fitting_method      15
uncertainty_available      15
replicate_information      10
raw_data_available          5
```

---

# Project 1 — ReactivityRatioBench

## 3.1 研究目标

建立一个真正可公平比较 RR 预测模型的 benchmark，同时回答：

1. **condition 是否真正提高 RR prediction？**
2. **当前误差到底来自 model 还是 data quality？**
3. **OOD 场景下已有模型还能否工作？**
4. **\(\ln r\) 或 \(\Delta\Delta G^\ddagger\) 是否优于直接预测 \(r\)？**
5. **physics baseline + Δ-learning 是否比纯数据驱动模型更稳健？**

---

## 3.2 P1 数据来源与使用优先级

### Tier A — 主训练数据

#### 1. CoPolDB

提取：

```text
monomer identities
SMILES
r1/r2
reference DOI
```

用途：大规模基础 RR 数据。

处理：

- 不直接视为 primary measurement。
- 根据 DOI 回溯原始论文。
- 与 Polymer Handbook / Chen / Safari 做 provenance 去重。

---

#### 2. PolyCARP

优先读取 repo 内 released CSV，而不是重新跑 VLM extraction：

```text
src/copolextractor/extracted_reactions.csv
processed_data.csv / paper_dataset
```

重点字段：

```text
monomer structures
r1/r2
temperature
solvent
polymerization mechanism
source paper
```

用途：**condition-aware 模型的核心数据源。**

同时可以复用其 extraction pipeline 处理你手动下载的订阅 PDF。

---

#### 3. Coxwell 2026 Structure–Reactivity SI

直接下载 SI ZIP。

包含：

```text
reactivity ratio dataset
SMILES
literature source
computational parameters
descriptors
```

用途：

- primary-literature-oriented 高质量 subset
- 与 CoPolDB/PolyCARP 交叉验证
- OOD external validation

---

### Tier B — 物理 descriptor / baseline

#### 4. CopDDB

从 GitHub clone CSV 数据。

重点字段：

```text
DE_barrier
DE_reaction
E_Rad_SOMO
E_Rad_LUMO
E_Mon_HOMO
E_Mon_LUMO
orbital gaps
%Vbur
volume
reaction geometry
```

处理：

CopDDB 基础单元是：

```text
(radical M1*, monomer M2)
```

而 RR \(r_1\) 需要比较：

```text
(M1*, M1)
(M1*, M2)
```

所以构建：

\[
\Delta x_1=x(M_1^*,M_1)-x(M_1^*,M_2)
\]

对 \(r_2\) 类似。

用途：

- physics descriptor baseline
- Δ-learning baseline
- SHAP/chemical interpretation

---

#### 5. Chen & Jackson

优先复用：

```text
dataset/
DFT_barrier_result/
descriptor_gen/
```

不要一开始重跑全部 DFT。

第一阶段只需读取作者已经提供的：

- experimental labels
- descriptors
- barrier results

并建立与 `RRMaster` 的 monomer/pair 映射。

---

### Tier C — 高质量外部验证

#### 6. Alexander 2026

6 个标准化机器人实验 monomer pairs 建议**全部从训练集排除**，作为 high-quality external test。

这样可以检验：

> literature-trained model 对标准化现代实验的真实泛化能力。

---

## 3.3 P1 最终数据集版本

建议发布 4 个版本：

### Dataset-U：Unique Measurements

每个 independent experiment 一行。

### Dataset-HQ：High Quality

仅 Quality A/B。

### Dataset-C：Condition Complete

要求：

```text
r1/r2
T
solvent
monomer structures
```

全部可用。

### Dataset-Legacy

允许缺 condition，用于和过去论文公平比较。

---

## 3.4 特征生成

### Structure-only

- ECFP4 2048 bit
- RDKit descriptors
- molecular GNN embedding

### Condition

- temperature standardized
- solvent RDKit descriptors / learned embedding
- BULK 单独 categorical token
- mechanism one-hot / embedding

### Physics

- Q–e（若可获得）
- CopDDB descriptors
- Chen/Jackson FMO descriptors

---

## 3.5 模型与实验矩阵

### Baseline

```text
ECFP + Random Forest
ECFP + XGBoost
RDKit descriptors + GBDT
MLP
GNN pair encoder
```

### Ablation

```text
Structure only
+ Temperature
+ Solvent
+ Temperature + Solvent
+ Physics descriptors
```

### Targets

```text
r
ln(r)
DeltaDeltaG‡
```

### Δ-learning

\[
\hat y=y_{Q-e}+\Delta_{ML}
\]

以及：

\[
\hat y=y_{physics\ descriptor\ baseline}+\Delta_{ML}
\]

---

## 3.6 Benchmark split

必须按 group split，防止 pair leakage。

1. Random measurement split
2. Pair split
3. One-monomer-out
4. **Both-monomers-out**
5. Chemical family / scaffold split
6. Publication/source split
7. Temporal split
8. Condition split

其中最重要的是：

\[
\boxed{\text{Both-monomers-out}}
\]

---

## 3.7 评价指标

### RR 参数

- MAE(ln r)
- RMSE(ln r)
- MAE(\(\Delta\Delta G^\ddagger\))
- \(R^2\)

### Copolymerization behavior

Mayo–Lewis：

\[
F_1=
\frac{r_1f_1^2+f_1f_2}
{r_1f_1^2+2f_1f_2+r_2f_2^2}
\]

定义：

\[
E_{curve}=\int_0^1
|F_1^{pred}(f)-F_1^{true}(f)|df
\]

### uncertainty

- NLL
- calibration curve
- 80/90/95% interval coverage

---

## 3.8 P1 Go/No-Go

**Go：**

- provenance 去重后仍有 >1500–2000 条 usable RR measurements；
- condition-complete subset >800；
- both-monomers-out split 可构造。

即使 condition-complete 少于预期，P1 仍可退化为：

> **Data-quality + provenance + OOD benchmark**

因此 P1 风险最低。

---

# Project 2 — Temperature-Dependent RR Thermodynamics

## 4.1 核心科学问题

不是直接学习：

\[
(M_1,M_2,T)\rightarrow r
\]

而是利用：

\[
\ln r=
-\frac{\Delta\Delta H^\ddagger}{RT}
+\frac{\Delta\Delta S^\ddagger}{R}
\]

学习：

\[
(M_1,M_2,\mathrm{solvent})
\rightarrow
(\Delta\Delta H^\ddagger,\Delta\Delta S^\ddagger)
\]

再恢复：

\[
r(T)=
\exp\left[-\frac{\Delta\Delta H^\ddagger}{RT}
+\frac{\Delta\Delta S^\ddagger}{R}\right]
\]

---

## 4.2 P2 数据来源

主要来自 P1 `RRMaster`，但重点使用以下来源做 temperature enrichment：

### 1. PolyCARP

因为包含 temperature + solvent + RR，优先统计：

```text
canonical_pair_id
solvent_id
temperature_K
r1
r2
source_doi
```

---

### 2. CoPolDB → DOI → 原文 PDF

CoPolDB 主要用于发现具有多个文献测量的 monomer pair。

对同一 pair：

1. 收集所有 DOI
2. 搜索不同温度测量
3. 从 PDF 表格/正文提取：

```text
T
solvent
r1/r2
method
uncertainty
```

---

### 3. Polymer Handbook

用于补 legacy multi-temperature series。

但必须追溯 primary reference，至少标记：

```text
secondary_only = True/False
```

---

### 4. 你下载的全文 PDF

优先检索关键词：

```text
temperature dependence
reactivity ratio at
°C
Kelen-Tudos
Fineman-Ross
Arrhenius
activation enthalpy
activation entropy
```

---

## 4.3 P2 数据处理

### Step 1：定义 temperature series

同一 series 至少要求：

```text
canonical_pair_id 一致
solvent 一致或均为 bulk
polymerization mechanism 一致
```

推荐附加：

```text
source/lab 尽量一致
fitting method 尽量一致
```

---

### Step 2：去除伪 multi-temperature series

不能简单把不同文献不同温度拼起来。

分两类：

#### Strict Series

```text
same pair
same solvent
same study/lab
same fitting method
>=3 temperatures
```

#### Extended Series

```text
same pair
same solvent
multiple studies allowed
>=3 temperatures
```

主结果用 Strict，Extended 做敏感性分析。

---

### Step 3：拟合 thermodynamic parameters

对每一个方向分别拟合：

\[
y=\ln r
\]

\[
x=1/T
\]

\[
y=a+bx
\]

得到：

\[
\Delta\Delta H^\ddagger=-bR
\]

\[
\Delta\Delta S^\ddagger=aR
\]

若有 RR uncertainty，进行 weighted least squares。

保存：

```text
ddh_kJ_mol
dds_J_mol_K
fit_r2
fit_rmse
n_temperature_points
temperature_span_K
parameter_uncertainty
```

---

### Step 4：质量控制

要求：

- 至少 3 个 distinct T；优先 ≥4
- 建议 temperature span ≥30 K
- \(\ln r\) vs \(1/T\) 明显非线性者单独标记
- 单个异常点不能自动删除；保留 robust regression 对比

---

## 4.4 P2 Go/No-Go：必须先统计再决定

### 强 Go

- Strict Series ≥150 个 monomer-pair-direction records
- 或 ≥75 个完整 binary pairs（两方向均有）

可做独立 ML paper。

### 条件 Go

- 50–150 个 series

适合：

- multi-task / physics-informed small-data ML
- Gaussian process
- descriptor regression
- family-level mechanistic analysis

### No-Go / 转型

- <50 个可靠 multi-temperature series

不要硬做深度学习。

转为：

> temperature-dependence meta-analysis + mechanistic descriptor study

或者暂缓 P2。

---

## 4.5 推荐实验

比较：

### Direct

\[
(M_1,M_2,T,S)\rightarrow\ln r
\]

### Physics-constrained

\[
(M_1,M_2,S)\rightarrow
(\Delta\Delta H^\ddagger,\Delta\Delta S^\ddagger)
\rightarrow\ln r(T)
\]

评价：

- within-range interpolation
- temperature extrapolation
- unseen-pair extrapolation

---

# Project 3 — ML Prior + Bayesian / Active Experimental Design

## 5.1 核心问题

如果已有 P1 模型提供：

\[
p(r_1,r_2|M_1,M_2,T,S)
\]

能否减少确定 RR 所需的实验数量？

目标不是重新发明 RR fitting，而是研究：

> **machine-learned prior 是否能改善 sequential experimental design？**

---

## 5.2 P3 核心数据源

### 1. Alexander et al. 2026 robotic RR dataset

这是最重要的数据源。

包含标准化测量的 6 个 monomer pairs：

```text
BMA–BA
BMA–St
St–BA
St–MMA
GMA–BA
GMA–St
```

用途：

- retrospective sequential-design benchmark
- 不需要重新做机器人实验
- 可以逐个“隐藏”实验点，模拟主动实验选择

---

### 2. IUPAC Reactivity-Ratios repo

直接复用：

```text
f0XFtestset.csv
IUPAC f0XF method in PythonV15.py
Converter f to F.xlsx
```

统一数据格式：

```text
f0
X
F
DeltaF
```

用途：

- 标准 likelihood / RR estimation
- 模拟实验数据
- 验证拟合代码

---

### 3. Lynd et al. / IUPAC recommendation examples

若 SI 中有原始或数字化 dataset，可作为额外 validation。

---

### 4. Wakiuchi et al. 2026

若 SI 有 gradient-flow curve / inline ATR-FTIR 数字数据：

可作为“dense trajectory”外部 case study。

如果只有图，不作为主训练源。

---

## 5.3 P3 数据处理

### Step 1：统一实验观测格式

```text
experiment_id
pair_id
f0
X
F
sigma_F
temperature
solvent
```

---

### Step 2：构造候选实验空间

例如：

```text
f0 = 0.05, 0.10, ..., 0.95
```

实验输出通过 terminal-model forward model 或真实机器人测量 lookup 获得。

---

### Step 3：构建 prior

对同一个 pair，P1 给出：

\[
\log r\sim\mathcal N(\mu,\sigma^2)
\]

比较：

1. noninformative prior
2. literature generic prior
3. P1 ML prior
4. oracle-like prior（上限参考）

---

### Step 4：Sequential experiment simulation

初始：

```text
D0 = empty 或一个随机实验
```

每轮：

1. 计算 posterior \(p(r_1,r_2|D_t)\)
2. 对所有候选 \(f_0\) 计算 expected information gain
3. 选择 \(f_0^*\)
4. 加入对应 observation
5. 更新 posterior

停止条件可为：

```text
JCI area < threshold
posterior SD(r1), SD(r2) < threshold
composition-curve uncertainty < threshold
```

---

## 5.4 Baseline 策略

必须比较：

- equally spaced feed compositions
- random
- classical optimal design / sensitivity-based
- Bayesian design without ML prior
- **Bayesian design with P1 ML prior**

---

## 5.5 评价指标

```text
number_of_experiments_to_threshold
posterior_error_r1
posterior_error_r2
JCI_area
composition_curve_error
failure_rate
```

尤其推荐报告：

\[
\mathrm{Efficiency}=
\frac{N_{baseline}-N_{method}}{N_{baseline}}
\]

---

## 5.6 P3 Go/No-Go

P3 对数据量要求不高，因此非常可行。

**Go 条件：**

- Alexander 数据能恢复逐实验点 \(f_0,X,F\)
- IUPAC fitting 可完全复现

即使只有 6 个 pairs，也足够做 methodological proof-of-concept；关键在于每个 pair 有多个实验点。

如果原始点不足：

- 使用 IUPAC synthetic simulator 生成不同 ground-truth \(r_1,r_2\)
- Alexander 数据作为真实 case study

---

# Project 4 — RR → Sequence → Tg / Property

## 6.1 核心科学问题

把两个孤立任务：

\[
\mathrm{structure}\rightarrow r
\]

和：

\[
\mathrm{structure/composition}\rightarrow T_g
\]

连接为：

\[
\boxed{
\mathrm{chemistry}
\rightarrow
\mathrm{reactivity}
\rightarrow
\mathrm{sequence}
\rightarrow
T_g
}
\]

重点回答：

> **由 RR 推导的 sequence information，能否解释 composition-only / Fox-type Tg 模型的 residual？**

---

## 6.2 P4 数据源

### Dataset A — P1 RRMaster

提供：

```text
monomer pair
r1/r2
T
solvent
uncertainty
```

优先使用 measured RR；无 measured RR 时可使用 P1 model prediction，但必须带 uncertainty flag。

---

### Dataset B — 共聚物 Tg 数据

优先顺序：

1. 你现有的 copolymer Tg 数据集
2. 已整理的 sequence–Tg 文献数据
3. Polymer Genome / literature-derived copolymer Tg（若能明确组成和 monomer pair）
4. 从当前全文 PDF 中补充同一 monomer pair 的 Tg / composition 数据

最低必需字段：

```text
monomer_1
monomer_2
composition
Tg
```

最好额外有：

```text
polymerization_temperature
solvent
conversion
molecular_weight
Tg_measurement_method
```

---

### Dataset C — Homopolymer Tg

每个 monomer 对应：

```text
Tg_A
Tg_B
```

用于 Fox / Gordon–Taylor baseline。

---

### Dataset D — 序列/组成理论验证

- Smith et al. 2018/2019：heteropolymer composition/drift
- Zhang et al. 2025：RR fingerprint / sequence regulation

主要用于公式与 qualitative validation，不一定直接作为训练标签。

---

## 6.3 P4 数据链接与处理

### Step 1：统一 monomer registry

必须使用与 P1 完全相同的：

```text
monomer_id
canonical_smiles
canonical_pair_id
```

不要重新建立一套 monomer naming。

---

### Step 2：匹配 Tg sample 与 RR

优先级：

#### Match-A：严格条件匹配

```text
same pair
polymerization T close
same solvent / bulk
```

#### Match-B：pair-level measured RR

```text
same pair
condition mismatch
```

#### Match-C：P1 predicted RR

无实验 RR 时使用 condition-aware P1 prediction。

增加：

```text
rr_source = measured_matched / measured_unmatched / predicted
rr_uncertainty
```

---

## 6.4 从 RR 生成 sequence descriptors

在 terminal model 下：

\[
P(A|A)=\frac{r_A f_A}{r_A f_A+f_B}
\]

\[
P(B|A)=\frac{f_B}{r_A f_A+f_B}
\]

\[
P(B|B)=\frac{r_B f_B}{r_B f_B+f_A}
\]

\[
P(A|B)=\frac{f_A}{r_B f_B+f_A}
\]

由此构造：

```text
P_AA
P_AB
P_BA
P_BB
mean_run_length_A
mean_run_length_B
blockiness
alternation_index
sequence_entropy
```

---

## 6.5 Conversion-dependent sequence

如果有 conversion，则不仅用初始 feed。

利用 Mayo–Lewis differential equation 模拟：

```text
f_A(X)
F_A(X)
```

在多个 conversion steps 计算 sequence statistics，然后形成：

```text
sequence_descriptor_initial
sequence_descriptor_mean
sequence_descriptor_final
sequence_descriptor_std_over_conversion
```

这比单点 \(r_1r_2\) 更有物理意义。

---

## 6.6 Tg baseline 与 Δ-learning

### Baseline 1 — Linear mixing

\[
T_g^{lin}=w_AT_{g,A}+w_BT_{g,B}
\]

### Baseline 2 — Fox

\[
\frac{1}{T_g^{Fox}}=
\frac{w_A}{T_{g,A}}+
\frac{w_B}{T_{g,B}}
\]

定义：

\[
\Delta T_g=T_g^{exp}-T_g^{Fox}
\]

训练：

\[
\Delta T_g=
f(\mathrm{sequence\ descriptors},
\mathrm{monomer\ features})
\]

最终：

\[
T_g^{pred}=T_g^{Fox}+\Delta T_g^{ML}
\]

---

## 6.7 最重要的消融

比较：

```text
A. composition + homopolymer Tg
B. A + r1,r2
C. A + r1*r2 / log(r1*r2)
D. A + derived sequence descriptors
E. A + conversion-dependent sequence descriptors
```

如果 D/E 显著优于 B/C，可以得到很强的科学结论：

> RR 的价值不是作为两个额外数字，而是通过 sequence statistics 将 polymerization kinetics 与 property prediction 连接起来。

---

## 6.8 Uncertainty propagation

若 RR 来自 P1 prediction：

\[
r_i\sim p(r_i)
\]

Monte Carlo sample：

```text
r1^(k), r2^(k)
    ↓
sequence descriptors^(k)
    ↓
Tg prediction^(k)
```

得到：

```text
Tg_mean
Tg_uncertainty_from_RR
```

可以进一步回答：

> RR uncertainty 对下游 Tg uncertainty 到底影响多大？

---

## 6.9 P4 Go/No-Go

先做 dataset overlap 统计。

### 强 Go

- ≥300 copolymer Tg samples 能匹配 measured RR
- 且覆盖 ≥50 monomer pairs

### 条件 Go

- measured RR overlap 100–300 samples
- 但 P1 model 对对应 pair 有较高 confidence

可做 measured + predicted RR 的 uncertainty-aware study。

### No-Go

- <100 Tg records 能可靠映射 monomer pair/composition

此时不要强行做实验数据论文；先扩充 Tg 文献数据。

---

# 7. PDF 全文抽取：统一的数据补全任务

你已经手动下载了核心全文 PDF，因此建议做一个统一的 structured extraction，而不是每个项目单独读论文。

## 7.1 每篇 RR 文献抽取 JSON

```json
{
  "reference": {
    "doi": "",
    "year": null,
    "title": ""
  },
  "systems": [
    {
      "monomer_1": "",
      "monomer_2": "",
      "smiles_1": "",
      "smiles_2": "",
      "temperature_C": null,
      "solvent": "",
      "initiator": "",
      "polymerization_method": "",
      "measurements": [
        {
          "r1": null,
          "r2": null,
          "r1_uncertainty": null,
          "r2_uncertainty": null,
          "fitting_method": "",
          "kinetic_model": "terminal",
          "f0_points": [],
          "conversion_points": [],
          "copolymer_composition_points": [],
          "table_or_figure": "",
          "page": null
        }
      ]
    }
  ]
}
```

---

## 7.2 抽取后的人工审核优先级

不是每条都人工审核，按价值分层：

### 必审

- multi-temperature series
- duplicate/conflicting RR
- Quality A/B 候选
- extreme RR
- OCR/VLM confidence 低的数据

### 抽样审

- 普通 legacy single RR
- 已与数据库完全一致的数据

---

# 8. 四个项目之间的数据复用关系

| 数据资产 | P1 | P2 | P3 | P4 |
|---|---:|---:|---:|---:|
| Monomer registry | ✓ | ✓ | ✓ | ✓ |
| RRMaster | ✓ | ✓ | prior | ✓ |
| Provenance / duplicate groups | ✓ | ✓ | — | ✓ |
| Temperature / solvent | ✓ | **核心** | ✓ | 可选 |
| Uncertainty / quality score | **核心** | ✓ | **核心** | ✓ |
| CopDDB/DFT descriptors | ✓ | ✓ | prior | 可选 |
| IUPAC \(f_0-X-F\) | quality check | — | **核心** | conversion 模拟 |
| Alexander robotic data | external test | — | **核心** | — |
| Copolymer Tg | — | — | — | **核心** |

---

# 9. 实际执行顺序

## Phase 1 — 先做数据底座

1. 下载/冻结各公开 dataset 版本
2. 建 monomer registry
3. 合并 CoPolDB / PolyCARP / Coxwell / Chen-Jackson
4. DOI/provenance 去重
5. 补全文 PDF 中的 temperature / solvent / fitting method
6. 建 quality score
7. 生成 RRMaster v1

输出：

```text
rrmaster_v1.parquet
monomer_registry.csv
source_lineage.csv
duplicate_groups.csv
```

---

## Phase 2 — 立刻并行做三个 feasibility statistics

### P1

```text
N_unique_measurements
N_unique_pairs
N_unique_monomers
N_condition_complete
N_quality_A_B
```

### P2

```text
N_pair_solvent_groups_with_>=3_temperatures
N_pair_solvent_groups_with_>=4_temperatures
median_temperature_span
```

### P4

```text
N_Tg_samples_with_measured_RR
N_unique_pairs_overlap
N_Tg_samples_with_condition_matched_RR
```

这三个数字决定 P2/P4 是否继续。

---

# 10. 最推荐的论文组合

## 第一篇：P1

**最稳妥，建议必须完成。**

核心贡献：

```text
cleaned provenance-aware RR dataset
+ rigorous OOD benchmark
+ condition-aware prediction
+ data-quality / uncertainty analysis
+ physics-informed target / delta-learning
```

---

## 第二篇：P2 或 P3

### 如果多温度数据足够

优先 **P2**。

### 如果多温度数据不足

立即切换 **P3**。

P3 不依赖大规模 literature data，因此是很好的备用第二篇。

---

## 第三篇：P4

前提是 RR–Tg overlap 足够。

P4 与现有 polymer Tg / sequence / delta-learning 工作协同最大，也最容易形成一个更大的 AI4Polymer 研究主线。

---

# 11. 最终推荐结论

如果按照“可行性 × 科学意义 × 数据成本 × 论文独立性”综合排序：

\[
\boxed{P1 > P3 \approx P2 > P4}
\]

但从长期研究主线看：

\[
\boxed{P1 \rightarrow P2 \rightarrow P4}
\]

故事最完整：

```text
monomer structure + condition
        ↓
reactivity ratio
        ↓
kinetic thermodynamics
        ↓
sequence distribution
        ↓
polymer property
```

而 **P3** 是一条相对独立但很实用的“ML → experimental design”支线，可以在 P2 数据量不足时优先推进。

---

# 12. 关键来源

- CoPolDB: https://www.copoldb.jp/
- Takahashi et al. 2024: https://doi.org/10.1039/D3PY01372C
- PolyCARP / copolymer-reactivity: https://github.com/lamalab-org/copolymer-reactivity
- Chen & Jackson repository: https://github.com/TheJacksonLab/ReactivityRatios
- Chen & Jackson 2025: https://doi.org/10.1021/acs.chemmater.5c01608
- Coxwell et al. 2026 Structure–Reactivity: https://doi.org/10.1021/acs.macromol.6c00233
- Coxwell et al. 2026 Better Models or Better Data: https://doi.org/10.1021/acs.macromol.6c01785
- CopDDB repository: https://github.com/hatanaka-lab/CopDDB
- CopDDB paper: https://doi.org/10.1039/D4DD00266K
- Alexander et al. 2026: https://doi.org/10.1039/D6SC02232D
- IUPAC Reactivity-Ratios code: https://github.com/IUPAC-Polymer-Division/Reactivity-Ratios
- IUPAC Recommendations 2025: https://doi.org/10.1515/pac-2024-0235
- Lynd et al. 2019: https://doi.org/10.1021/acs.macromol.8b01752
- van Herk & Liu 2024: https://doi.org/10.1002/mats.202400043
- Wakiuchi et al. 2026: https://doi.org/10.1021/acsomega.6c06915
