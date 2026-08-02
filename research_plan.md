# 耗散 SIDM 速度依赖截面 + 多模型排除图：研究计划

> **论文基础**：Schmidt, Fischer & Garny 2026, A&A, arXiv:2606.19428v1
> **目标期刊**：PRD（保底）/ PRL（若多模型排除图判别力足够）
> **竞争窗口**：Zhang & Yu 2026 已在 arXiv 抢弹性 SIDM 解释；Vegetti 团队后续会有更多致密天体样本。窗口在收窄，需尽快出初步结果。

---

## 0. 科学目标

在 Schmidt et al. 2026 的耗散 fSIDM 框架基础上，把 σ_T/m_χ 与 r_diss 从速度无关常数升级为**由具体暗 sector 拉氏量统一导出的速度依赖函数**，并用多模型对比 + 现有观测约束生成排除图，回答：

> **哪类暗 sector 模型（暗光子 brems / 原子 DM bound state / mirror DM）能同时解释 JVAS B1938+666 致密天体，且不被 Bullet Cluster / 星系团核 / 矮星系核约束排除？**

核心判别逻辑：速度依赖打破 Schmidt et al. 2026 Appendix G 的 (λ, μ) 重标度对称 → 不同模型在 σ_T(v)–r_diss(v) 平面上给出不同形状曲线 → 对 B1938+666 给出模型相关的可证伪参数区间。

---

## 1. 阶段划分

| 阶段 | 内容 | 执行方 | 计算资源 | 产出 |
|------|------|--------|----------|------|
| **P1** | 三类模型的微分截面文献整理 + σ_T(v)、r_diss(v) 解析/数值计算 | 我（本地） | 无 | `cross_sections/` 模块 + 函数曲线图 |
| **P2** | 流体模型（GravothermalSIDM 公开版）接入速度依赖 + 三类模型孤立晕演化 | 我（本地） | 普通工作站 | 三类模型的 t_coll、密度剖面、γ_2D 演化轨道 |
| **P3** | B1938+666 重标定（速度依赖打破对称，需逐模型做） | 我（本地） | 普通工作站 | 三类模型在 r_2D/rs 平面上的可解释区域 |
| **P4** | 观测约束叠加（Bullet Cluster、星系团核、矮星系核）→ 多模型排除图 | 我（本地） | 无 | **核心论文图：σ_T–r_diss 多模型排除图** |
| **P5** | N-body 验证（OpenGadget3 / 自实现 fSIDM 耗散模块） | **你/课题组** | **HPC** | 验证流体模型结论的稳健性 |
| **P6** | 撰写 + 投稿 | 共同 | 无 | 论文初稿 |

**当前执行**：P1 → P4，全部可在本地完成。P5 标记为"等待资源确认"。

---

## 2. P1：三类模型的微分截面与速度依赖参数

### 2.1 模型选择与文献依据

| 模型 | 暗 sector 内容 | 耗散通道 | 主要文献 |
|------|----------------|----------|----------|
| **M1: 暗光子 bremsstrahlung** | DM χ + 暗光子 A' (质量 m_{A'}) | χχ → χχ + A' | Lankester-Broche & Pradler 2026 (JCAP 2026, 034); Essig et al. 2019 (PRL 123, 121102) |
| **M2: 原子 DM (bound state)** | DM "原子" 含束缚态 | 碰撞激发 + 辐射衰变 | Kaplan et al. 2010 (JCAP 2010, 021); Wise & Zhang 2014; Gurian et al. 2022 |
| **M3: Mirror DM** | 镜像重子 + 镜像光子 | 镜像 brems + 镜像原子线冷却 | Foot 2004 (IJMPD 13, 2161); Foot & Silagadze 2013 |

### 2.2 待执行任务

- [ ] **M1**：从 Lankester-Broche & Pradler 2026 提取 dσ/dΩ dk0 的解析形式（暗光子软辐射近似），数值积分得 σ_T(v; m_{A'}, α_D, m_χ) 与 r_diss(v; 同参数)
- [ ] **M2**：从 Kaplan et al. 2010 / Gurian et al. 2022 提取原子 DM 的激发截面 + 衰变能损，构造 r_diss(v; E_b, α_D)
- [ ] **M3**：从 Foot 2004 提取镜像 DM 的 brems + 线冷却，整理成同一参数化形式
- [ ] 三类模型在同一 (v, σ_T/m_χ, r_diss) 空间画叠加图，识别判别区域
- [ ] **检查清单**：每个模型的 σ_T(v)、r_diss(v) 在 v → 0 和 v → ∞ 的渐近行为；是否有特征速度尺度（m_{A'} 或 E_b 引入的 v_* = m_{A'}/m_χ 等）

### 2.3 产出

- `src/cross_sections/` — 三个模型各自的 Python 模块
- `figures/P1_sigma_r_diss.png` — 三类模型在 (σ_T/m_χ, r_diss) 平面，颜色编码 v
- `notes/P1_derivations.md` — 关键推导笔记

---

## 3. P2：流体模型 + 速度依赖孤立晕演化

### 3.1 方法

使用公开的 **GravothermalSIDM**（https://github.com/kboddy/GravothermalSIDM，Nishikawa et al. 2020; Outmezguine et al. 2023; Gad-Nasr et al. 2024），按 Schmidt et al. 2026 Appendix D 接入冷却率

$$C(\rho, \nu) = \frac{8\sqrt{\pi}}{?} \frac{\sigma_T}{m_\chi}(r_{\rm diss} - 1) \rho^2 \nu^3$$

（精确系数见原论文 Eq. 19）。

**速度依赖扩展**：原论文冷却率假设 σ_T、r_diss 是常数。速度依赖下需修改流体方程：
- 局部相对速度分布取 Maxwell-Boltzmann（温度 ν²）
- σ_T、r_diss 替换为热平均 ⟨σ_T v⟩、⟨r_diss σ_T v⟩/⟨σ_T v⟩
- 热传导率 κ 同样改为热平均形式

### 3.2 待执行任务

- [ ] 克隆 GravothermalSIDM，确认其冷却率接口可注入任意 C(ρ, ν)
- [ ] 实现 P1 三个模型的 ⟨σ_T(v)⟩、⟨r_diss(v) σ_T(v)⟩ 热平均（Maxwell-Boltzmann 权重）
- [ ] 跑 NFW 初始条件（同 Schmidt et al. 2026: r_s=3.6 kpc, ρ_0=7.09×10⁻⁶ M_⊙/kpc³）三类模型 + 弹性对照
- [ ] 提取 t_core、t_coll、密度剖面、γ_2D(r=0.2/0.8/3.2 kpc) 演化轨道
- [ ] 与 Schmidt et al. 2026 Fig. 2/4/9 对比验证（弹性情形应复现）

### 3.3 产出

- `src/fluid_runner/` — 接入速度依赖的流体模型封装
- `figures/P2_evolution.png` — 三类模型的密度/ν/γ_2D 演化
- `data/P2_collapse_times.csv` — t_core、t_coll vs 模型参数

---

## 4. P3：B1938+666 重标定（速度依赖打破对称）

### 4.1 关键变化

Schmidt et al. 2026 Appendix G 的 (λ, μ) 重标度对称**只在速度无关时成立**。速度依赖引入特征速度 v_* = m_{A'}/m_χ（M1）、v_* = √(E_b/m_χ)（M2）后：

- λ 仍可匹配空间尺度，但 μ 不再能自由调——σ_T(v) 的形状被锁定
- Table 2 的重标定需逐模型重新求解：固定 (m_χ, m_{A'}, α_D)，扫描初始 NFW 参数 (r_s, ρ_0) 找到能匹配 M(r<20 pc) 与 M(r<90 pc) 的解

### 4.2 待执行任务

- [ ] 对每个模型，在 (r_s, ρ_0, 模型参数) 三维空间扫描
- [ ] 匹配条件：M(r<20 pc) = (4.25±0.21)×10⁵ M_⊙，M(r<90 pc) = (1.167±0.039)×10⁶ M_⊙
- [ ] 演化时间约束：t_evo ≤ t(z_obs=0.881) = 6.37 Gyr
- [ ] 输出每个模型可解释 B1938+666 的参数区域

### 4.3 产出

- `figures/P3_B1938_regions.png` — 三类模型在参数空间的可解释区域
- `data/P3_rescaled_params.csv` — 重标定后的物理参数

---

## 5. P4：多模型排除图（核心论文图）

### 5.1 观测约束收集

| 约束 | 物理量 | 限制对象 | 文献 |
|------|--------|----------|------|
| Bullet Cluster (1E 0657-56) | 子晕质量比偏移 | σ/m at v~3000 km/s | Randall et al. 2008 (ApJ 679, 1177); σ/m < 1.25 cm²/g |
| 星系团核 (Aquila, A3827 等) | 核心形态 | σ/m at v~1000 km/s | Sagunski et al. 2021 (JCAP 2021, 024); σ/m < ~1-3 cm²/g |
| 矮星系核 (Draco, Fornax 等) | 核密度 | σ/m at v~50 km/s | Read et al. 2018 类；σ/m ~ 0.5-50 cm²/g 区间 |
| B1938+666 致密天体 | 质量分布 | 耗散 + 高 σ/m at v~50 km/s | Vegetti et al. 2026; Powell et al. 2025 |

### 5.2 待执行任务

- [ ] 把每个约束转化为 (v, σ_T/m_χ, r_diss) 空间的不等式
- [ ] 在三类模型的 (m_χ, m_{A'}, α_D)（或对应参数）平面上叠加约束
- [ ] 输出**主图**：以 σ_T(v=50 km/s) 为横轴、r_diss(v=50 km/s) 为纵轴，三类模型曲线 + 约束带 + B1938+666 可解释区
- [ ] 判别力评估：三类模型在主图上是否占据不重叠区域？是否有模型被完全排除？

### 5.3 产出

- `figures/P4_exclusion_main.png` — **论文核心图**
- `figures/P4_exclusion_param_space.png` — (m_{A'}, α_D) 平面补充图
- `notes/P4_discrimination_analysis.md` — 判别力评估与 PRL 可行性判断

---

## 6. P5：N-body 验证（需要 HPC，等待资源）

### 6.1 待执行任务

- [ ] 确认是否有 OpenGadget3 访问权（联系 Schmidt/Fischer/Garny 团队）
- [ ] 若无：在 AREPO/GADGET-4 上实现耗散 fSIDM 模块（参考 Fischer et al. 2021, 2024, 2026 的算法描述）
- [ ] 选 P4 中 1-2 个判别力最强的模型参数点，跑 5×10⁶ 粒子孤立晕
- [ ] 验证流体模型预测的 t_coll、密度剖面、γ_2D 轨道
- [ ] 若流体模型定性结论被推翻，回 P2 调整

### 6.2 资源需求

- 单次 5×10⁶ 粒子孤立晕模拟：~10⁴-10⁵ CPU 核时（参考 Schmidt et al. 2026）
- 3-5 个参数点验证：~10⁵ 核时总量
- 存储：每快照 ~GB 级，需 TB 级工作盘

### 6.3 待确认

> **请告知：**
> 1. 你/课题组是否有 OpenGadget3 或等价 fSIDM 代码访问权？
> 2. 是否有 HPC 配额（SLURM/PBS 集群，~100 核×1 周量级）？
> 3. 若没有代码但有 HPC，是否愿意投入工程量在 AREPO 上重实现？

---

## 7. P6：撰写与投稿

### 7.1 论文结构（草案）

1. Introduction — 引 Schmidt et al. 2026 的速度无关结果 + 重标度对称，指出其局限
2. Velocity-dependent cross sections for three dark sector models — P1
3. Gravothermal evolution with velocity dependence — P2 + P5
4. Re-interpreting JVAS B1938+666 — P3
5. Multi-model exclusion — P4
6. Discussion & conclusions

### 7.2 PRL 与 PRD 的判别标准

- **PRL**：P4 主图上三类模型有清晰不重叠区域，且至少一类被现有观测完全排除 → "first discrimination among dark sector models via compact-object lensing"
- **PRD（保底）**：P4 主图有模型相关差异但判别力不绝对 → 仍是有价值的模型约束工作

P4 完成后由 P4 判别力评估决定投稿目标。

---

## 8. 时间表（估计）

| 阶段 | 内容 | 预估时长 | 依赖 |
|------|------|----------|------|
| P1 | 微分截面整理 + σ_T(v)、r_diss(v) 计算 | 1-2 周 | 无 |
| P2 | 流体模型接入 + 三类模型演化 | 1-2 周 | P1 |
| P3 | B1938+666 重标定 | 1 周 | P2 |
| P4 | 多模型排除图 + 判别力评估 | 1-2 周 | P1-P3 |
| P5 | N-body 验证 | 4-8 周（含代码准备） | P4，等待资源 |
| P6 | 撰写投稿 | 2-3 周 | P5 |

**P1-P4 流体模型路径**：5-7 周可出初步结果，足以判断 PLR 可行性。

---

## 9. 当前进度

- [x] 计划制定（本文件）
- [x] **P1 完成**（初步）：六场景文献整理（LB2026 统一框架）+ σ_T(v)、r_diss(v) 参数化实现 + 三模型对比图
- [x] **P2 完成**：GravothermalSIDM 接入 + 热平均 + 耗散 Halo 类 + 四模型演化（弹性 t_core=4.35 Gyr 匹配 Schmidt）
- [x] **P3 完成**：B1938+666 重标定（18148 匹配点，紧凑解区域，对称性破缺在低 v 被 kinematically suppressed）
- [x] **P4 完成**：多模型排除图（三模型斜率 σ_high/σ_low 不同 → 可观测判别力存在，所有模型均 viable）
- [ ] P5（等待资源确认）
- [ ] P6

### P1 关键发现（2026-07-25）

1. **LB2026 是统一框架**：Lankester-Broche & Pradler 2026 (arXiv:2509.12317)
   不是单一模型，而是六类耗散 SIDM 场景的统一 Born 极限推导，给出
   ready-to-use 闭式表达式。比原计划的"三模型分散提取"更干净。

2. **对称性破缺的物理来源已明确**：
   - 无质量发射：r_diss = const（恢复 Schmidt et al. 重标度对称）
   - 有质量发射：r_diss(v) 引入特征速度 v* = √(2m_{φ,V}/μ)
   - v < v* 时 Boltzmann 抑制，v > v* 时恢复效率

3. **天体物理相关参数空间**：m_χ = 10 GeV 时，v* 落在 100-1000 km/s
   需要 m_{φ,V} ~ 0.3-30 keV（light mediator SIDM 典型参数）

4. **三模型判别力已验证**（figures/P1_sigma_r_diss.png）：
   - M1（暗光子 m_V=1.1 keV, v*≈200）：矮星系无耗散，星系团强耗散
   - M2（标量 m_φ=6.95 keV, v*≈500）：MW 尺度才开始耗散
   - M3（无质量对照）：恒定 r_diss=1.05（Schmidt et al. 极限）

### P1 局限（待 P2 改进）

1. r_diss 振幅 C=0.05 是 fiducial；需用 LB2026 Sec 5.1.2 完整四极矩表达式校准
2. 仅四极矩通道（全同粒子）；偶极（可区分两分量 DM）未实现
3. Boltzmann 近似；需数值积分 Eq (3.28) + C_φ/C_V 因子

### P2 完成（2026-07-25）

**流体模型实现**：
- 获取 GravothermalSIDM 源码（WebFetch from raw.githubusercontent.com，master 分支）
- 实现热平均模块 (`src/fluid_runner/thermal_avg.py`)：Maxwell-Boltzmann ⟨σ_T v⟩、⟨r_diss σ_T v⟩
- 实现耗散 Halo 包装类 (`src/fluid_runner/dissipative_halo.py`)：注入速度依赖 σ_m(T)，加冷却率 C_cool

**弹性验证**（与 Schmidt et al. 2026 Fig 2 对比）：
- NFW r_s=3.6 kpc, ρ_0=7.09e-3, σ_m=50, w=100, n_shells=100, t_epsilon=1e-2
- **t_core = 4.35 Gyr**（Schmidt ~4-5 Gyr）✓
- 流体模型在深坍缩阶段时间步趋零（数值奇点），用 rho_factor_end=1000 截断

**四模型结果**（data/P2_collapse_times.csv）：
| 模型 | t_core [Gyr] | ρ_min | t_final [Gyr] |
|------|-------------|-------|---------------|
| Elastic | 4.35 | 1.73e-2 | 6.02 |
| Const r_diss=1.05 | 0.90 | 3.40e-2 | 3.09 |
| M1 dark photon (v*=200) | 6.40 | 3.14e-2 | 7.49 |
| M2 scalar (v*=500) | 4.23 | 3.13e-2 | 7.49 |

**关键发现**：
1. 耗散加速坍缩（const r_diss=1.05 的 t_core 比弹性快 5x）— 与 Schmidt Fig 3 一致
2. M1 比弹性慢（v*=200 km/s 在矮星系尺度 r_diss≈1）— **打破重标度对称**
3. M2 与弹性相近（v*=500 km/s 高于晕速度，耗散被抑制）
4. M3（无质量）因 Rutherford 发散数值不稳定，已用 const_rdiss_1.05 替代

**产出文件**：
- `figures/P2_evolution.png` — 四模型密度/ν/γ_2D 演化对比
- `data/P2_collapse_times.csv` — t_core 等观测量
- `notes/P2_validation.md` — 验证笔记

---

## 10. 风险与备选

| 风险 | 备选 |
|------|------|
| GravothermalSIDM 流体模型在速度依赖下不稳定 | 自实现简化流体方程（参考 Balberg et al. 2002） |
| 某模型微分截面文献不完整（如原子 DM 激发截面） | 退化为两模型对比（M1 + M3），仍可发 |
| P4 判别力不足 | 加 P2 中的速度依赖 smoking gun（γ_2D 轨道差异）作为补充判据 |
| OpenGadget3 不可得 | AREPO 重实现，或仅用流体模型 + caveat 投 PRD |
