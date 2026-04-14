# 计算机博士后攻克 KRAS G12D 靶向药：从 AI 设计到实物合成的完整路线图

> 调研时间：2025–2026年  
> 目标：不只是计算设计——**真正造出 KRAS G12D 靶向药候选分子**  
> 关键词：KRAS G12D、AI 药物设计、化学合成、生物验证、临床前开发  
> 定位：计算机博士 → 生物医学博士后，发挥计算优势 + 补齐实验短板

---

## 目录

1. [背景与机遇](#1-背景与机遇)
2. [KRAS G12D 靶点全景](#2-kras-g12d-靶点全景)
3. [第一阶段：知识体系构建（0–3 个月）](#3-第一阶段知识体系构建03-个月)
4. [第二阶段：确定差异化研究方向（3–6 个月）](#4-第二阶段确定差异化研究方向36-个月)
5. [第三阶段：计算设计与虚拟筛选（6–12 个月）](#5-第三阶段计算设计与虚拟筛选612-个月)
6. [第四阶段：化学合成——把分子造出来（9–18 个月）](#6-第四阶段化学合成把分子造出来918-个月)
7. [第五阶段：体外生物学验证（12–20 个月）](#7-第五阶段体外生物学验证1220-个月)
8. [第六阶段：体内药效与药代动力学（18–28 个月）](#8-第六阶段体内药效与药代动力学1828-个月)
9. [第七阶段：先导化合物优化闭环（贯穿全程）](#9-第七阶段先导化合物优化闭环贯穿全程)
10. [第八阶段：成果产出、专利与产业化（24+ 个月）](#10-第八阶段成果产出专利与产业化24-个月)
11. [实验室搭建与资源规划](#11-实验室搭建与资源规划)
12. [合作网络构建策略](#12-合作网络构建策略)
13. [风险管理与 Plan B](#13-风险管理与-plan-b)
14. [关键成功因素](#14-关键成功因素)
15. [参考资料](#15-参考资料)

---

## 1. 背景与机遇

### 1.1 为什么 KRAS G12D 是绝佳切入点

KRAS G12D 曾被称为"不可成药"靶点（undruggable），但 2024–2026 年发生了历史性突破：

| 药物 | 机制 | 最新进展（2025–2026） | 关键数据 |
|------|------|----------------------|---------|
| **Zoldonrasib**（RMC-9805） | RAS(ON) 三复合物共价抑制剂 | FDA 突破性疗法认定（2026） | NSCLC ORR 61%，DCR 89% |
| **Setidegrasib**（ASP3082） | 首个 KRAS G12D 蛋白降解剂 | Phase I 全球数据发表（NEJM, 2026） | NSCLC ORR 36%，PFS 8.3个月 |
| **GFH375**（VS-7375） | 口服 G12D ON/OFF 双态抑制剂 | 全球首个 G12D Phase III（2026） | PDAC ORR 52%，DCR 100% |
| **Daraxonrasib**（RMC-6236） | 多靶 RAS(ON) 抑制剂 | FDA 突破性疗法认定（2025） | PDAC mPFS 9个月 |
| **ARV-806** | PROTAC 降解剂 | 临床前数据（AACR 2025） | >25× 优于现有抑制剂 |
| **HRS-4642** | 非共价 G12D 抑制剂 | ESMO 2025 数据 | 多线治疗后仍有效 |

**但目前仍无药物获批上市**，且耐药机制未解——这正是计算背景的人最大的机会窗口。

### 1.2 计算机博士的独特优势

传统药化团队的短板恰好是你的长板：

| 传统药化团队的瓶颈 | 你的计算优势 |
|-------------------|------------|
| 化学空间探索依赖经验 | 生成式 AI 可遍历数十亿分子 |
| 耐药突变响应滞后 | ML 可提前预测耐药并设计规避策略 |
| 构效关系（SAR）分析手动 | GNN/Transformer 可自动学习 SAR |
| 多目标优化困难 | 强化学习/贝叶斯优化可同时优化活性+成药性+选择性 |

### 1.3 核心理念：计算-实验闭环

本路线图的核心不是"做完计算交给别人"，而是：

```
AI 设计 → 化学合成 → 生物测试 → 数据反馈 → 模型更新 → 新一轮设计
     ↑                                                        ↓
     └────────────────── 闭环迭代 ──────────────────────────────┘
```

**你必须深度参与每一个环节**，而不仅仅是计算部分。

---

## 2. KRAS G12D 靶点全景

### 2.1 生物学背景

**KRAS（Kirsten Rat Sarcoma Viral Oncogene Homolog）** 是人类癌症中最常见的致癌基因之一。

**正常状态**：KRAS 在 GDP（失活态）和 GTP（激活态）之间循环切换：
```
                SOS1（GEF）
    KRAS-GDP  ─────────────→  KRAS-GTP
    （OFF）    ←─────────────  （ON）
                 GAP（NF1）
                内源 GTPase
```

**G12D 突变的本质**：
- 第 12 位甘氨酸（Gly, G）→ 天冬氨酸（Asp, D）
- Asp 的负电荷侧链阻碍 GAP 催化的 GTP 水解
- KRAS 被锁死在 GTP 结合的 **持续激活态**
- 下游 MAPK（RAF→MEK→ERK）和 PI3K/AKT/mTOR 通路持续激活
- 导致细胞不受控增殖

**临床分布**：
| 癌种 | KRAS G12D 突变频率 |
|------|-------------------|
| 胰腺导管腺癌（PDAC） | ~36%（最高频） |
| 结直肠癌（CRC） | ~12% |
| 非小细胞肺癌（NSCLC） | ~4% |
| 子宫内膜癌 | ~5% |

### 2.2 蛋白质结构关键特征

**可用晶体结构（PDB）**：

| PDB ID | 描述 | 分辨率 | 年份 | 特点 |
|--------|------|--------|------|------|
| **7RPZ** | G12D + MRTX1133 | 高分辨率 | 2022 | 非共价抑制剂复合物，必读结构 |
| **7F0W** | G12D + GDP（Switch I 开放） | 高分辨率 | 2022 | 关键构象参考 |
| **9L6F** | G12D + ASP3082 | 3.18 Å | 2025 | 首个降解剂复合物结构 |
| **9P45** | G12D + BBO-11818 | 1.35 Å | 2025 | 超高分辨率，非共价抑制剂 |
| **9C3K** | G12D/M67R + GDP | 1.70 Å | 2025 | 抑制突变体，机制研究 |
| **6GJ7** | G12D + BI-2852 | — | 2019 | 早期抑制剂复合物 |
| **4EPR** | G12D + GDP | — | 2012 | 经典参考结构 |

**关键结构区域**（必须理解）：
- **P-loop（残基 10–17）**：G12D 突变所在位置
- **Switch I（残基 30–38）**：效应蛋白结合界面
- **Switch II（残基 60–76）**：GAP 结合界面，allosteric pocket 所在
- **Asp12 侧链**：负电荷，可形成离子相互作用——设计抑制剂的关键把手

### 2.3 当前靶向策略分类

```
KRAS G12D 靶向策略
├── 直接抑制
│   ├── 共价抑制剂（利用 Asp12 附近亲核位点）
│   ├── 非共价抑制剂（MRTX1133 类，离子+氢键网络）
│   └── 三复合物共价（Zoldonrasib，锁定 OFF 态）
├── 蛋白降解
│   ├── 分子胶水/降解剂（Setidegrasib/ASP3082）
│   └── PROTAC（ARV-806，E3 连接酶介导降解）
├── 间接策略
│   ├── SOS1 抑制（阻断上游激活）
│   ├── SHP2 抑制（阻断上游激活）
│   └── 下游联合（MEK/ERK 抑制剂联合）
└── 新兴策略
    ├── 多靶 RAS 抑制（Daraxonrasib，覆盖 G12X/G13X/Q61X）
    ├── RNA 干扰/ASO（基因水平沉默）
    └── 免疫联合（G12D 新抗原疫苗 + 检查点抑制剂）
```

---

## 3. 第一阶段：知识体系构建（0–3 个月）

### 3.1 生物学基础（必须补的课）

**核心知识模块**：

| 模块 | 内容 | 推荐资源 | 优先级 |
|------|------|---------|--------|
| RAS 信号通路 | MAPK 和 PI3K/AKT 通路 | Molecular Biology of the Cell (Alberts) 第 15 章 | ⭐⭐⭐⭐⭐ |
| 蛋白质结构基础 | 氨基酸、二级/三级/四级结构 | Lehninger Biochemistry 第 4 章 | ⭐⭐⭐⭐⭐ |
| 癌症生物学 | 致癌基因、肿瘤抑制基因 | Cancer Biology (Weinberg) | ⭐⭐⭐⭐ |
| 药物化学基础 | 成药性、ADMET、PK/PD | Medicinal Chemistry (Patrick) | ⭐⭐⭐⭐ |
| 有机化学基础 | 反应机理、官能团 | Organic Chemistry (Clayden) 精选章节 | ⭐⭐⭐ |
| 细胞生物学实验 | Western blot、细胞培养、MTT | YouTube + 实验室 hands-on | ⭐⭐⭐⭐⭐ |

**实验技能快速入门**（必须亲手做）：

| 实验技术 | 学习方式 | 时间投入 | 为什么必须学 |
|---------|---------|---------|------------|
| 细胞培养 | 跟实验室技术员学 | 2 周 | 后续所有细胞实验的基础 |
| Western Blot | 导师/博后演示 + 自己做 | 1 周 | 验证通路抑制的核心手段 |
| MTT/CCK-8 增殖实验 | 1 天学会操作 | 2 天 | 测细胞毒性/IC50 |
| 基础分子克隆 | 实验手册 + 实操 | 2 周 | 理解质粒、突变体构建 |
| PyMOL / ChimeraX | 在线教程 | 3 天 | 蛋白结构可视化与分析 |

### 3.2 文献精读（Top Priority）

**必读论文清单**：

**靶点与结构**：
1. "Discovery of MRTX1133, a Noncovalent, Potent, and Selective KRAS G12D Inhibitor" — *J. Med. Chem.*
2. "Targeting the untargetable: accelerated discovery of KRAS G12D inhibitors" — *ScienceDirect*, 2025
3. "New exploration of KRAS G12D inhibitors and the mechanisms of resistance" — *PMC*, 2025

**降解剂策略**：
4. "Discovery of KRAS (G12D) selective degrader ASP3082" — *Nature Chemical Biology*, 2025
5. ARV-806 PROTAC 临床前数据 — *AACR-NCI-EORTC*, 2025

**临床数据**：
6. Zoldonrasib Phase I 结果 — *AACR 2025*
7. Setidegrasib Phase I 结果 — *NEJM*, 2026
8. GFH375 Phase I/II 数据 — *ASCO/ESMO 2025*

**AI 药物设计方法**：
9. DiffSBDD: Structure-based Drug Design with Equivariant Diffusion Models — *GitHub + 论文*
10. "Quantum biological convergence: quantum computing accelerates KRAS drug discovery" — *Nature Signal Transduction*, 2025
11. "Machine Learning-Driven Drug Repurposing for KRAS G12C and KRAS G12D" — *bioRxiv*, 2025
12. "Optimizing drug design by merging generative AI with a physics-based approach" — *PMC*, 2025

### 3.3 同时做的事：熟悉实验室环境

> **关键原则**：不要等计算做完再进实验室。从第一天起就泡在实验室。

- **第 1 周**：与实验室负责人确认你可以使用的仪器、试剂和经费
- **第 2–4 周**：跟一个有经验的博后/PhD 学习细胞培养（传代、冻存、复苏）
- **第 5–8 周**：独立完成 MIA PaCa-2（KRAS G12D）细胞系的培养和增殖实验
- **第 9–12 周**：掌握 Western Blot，能独立检测 p-ERK / total ERK

---

## 4. 第二阶段：确定差异化研究方向（3–6 个月）

### 4.1 五条差异化路径

#### 路径 A：AI 驱动的全新分子生成（De Novo Design）⭐⭐⭐⭐⭐

**核心思路**：用生成式 AI（扩散模型/GAN/VAE/Transformer）直接设计全新的 G12D 抑制剂分子。

- **工具链**：DiffSBDD、Pocket2Mol、TargetDiff、DrugGPT、REINVENT
- **你的优势**：深度学习模型架构设计、训练优化、采样策略
- **差异化点**：结合 G12D 特有的 Asp12 残基电荷特性，设计 mutation-specific 生成条件
- **到实物的路径**：生成分子 → 合成可行性评估 → CRO 合成 → 你做生物测试

#### 路径 B：耐药机制预测与克服策略 ⭐⭐⭐⭐⭐

**核心思路**：当前 G12D 抑制剂最大问题是耐药，用 ML 预测耐药突变并提前设计克服策略。

- **方法**：分子动力学（MD）模拟 + 深度突变扫描（DMS）+ 图神经网络（GNN）
- **数据源**：COSMIC 数据库、临床试验耐药突变谱
- **到实物的路径**：预测耐药突变 → 构建突变体细胞系 → 测试你设计的分子是否抗耐药

#### 路径 C：PROTAC/分子胶水的计算设计 ⭐⭐⭐⭐

**核心思路**：蛋白降解剂是最前沿策略，三元复合物设计极其复杂，正是计算的用武之地。

- **工具**：AlphaFold 3、PROTAC-Model、DeepPROTACs
- **到实物的路径**：设计三元复合物 → 合成 PROTAC → Western Blot 验证 KRAS 蛋白降解

#### 路径 D：多组学数据整合与联合靶点发现 ⭐⭐⭐⭐

**核心思路**：用 AI 发现与 G12D 合成致死的基因，指导联合用药。

- **数据源**：TCGA、DepMap、CCLE
- **到实物的路径**：发现联合靶点 → 购买已有抑制剂 → 联合用药实验

#### 路径 E：端到端自动化药物发现平台 ⭐⭐⭐

**核心思路**：构建全自动 AI 管线，以 KRAS G12D 为示范。

- **到实物的路径**：平台输出候选分子 → 合成验证 → 发表方法+验证论文

### 4.2 推荐组合策略

> **路径 A（生成式设计）+ 路径 B（耐药穿透）**
> 
> 用生成式 AI 设计新分子，同时预测并规避耐药，形成 **"AI 设计 + 耐药穿透"** 的独特定位。

**关键：每条路径都必须有明确的"到实物的路径"**，不能止步于计算。

---

## 5. 第三阶段：计算设计与虚拟筛选（6–12 个月）

### 5.1 计算环境与工具

| 类别 | 工具 | 用途 |
|------|------|------|
| 蛋白结构预测 | AlphaFold 3, RoseTTAFold | 预测 KRAS G12D-药物复合物结构 |
| 分子对接 | AutoDock Vina, Gnina, Glide | 虚拟筛选与打分 |
| 分子动力学 | GROMACS, OpenMM, AMBER | 验证结合稳定性与动态行为 |
| 分子生成 | DiffSBDD, Pocket2Mol, REINVENT | AI 生成新化学实体 |
| 化学信息学 | RDKit, DeepChem | 分子描述符、ADMET 预测 |
| ML 框架 | PyTorch, PyTorch Geometric | 图神经网络、生成模型 |
| 蛋白质语言模型 | ESM-2/3, ProtTrans | 突变效应预测、蛋白表示学习 |
| 合成可行性 | AiZynthFinder, ASKCOS (MIT) | **评估生成分子的可合成性** |
| 数据库 | PDB, ChEMBL, ZINC, Enamine REAL | 结构数据与化合物库 |

### 5.2 关键计算实验路线

#### Step 1：靶点结构准备（2 周）

```
操作流程：
1. 从 PDB 下载多个 G12D 结构（7RPZ, 9P45, 7F0W, 9L6F）
2. 用 PyMOL 对齐所有结构，分析构象差异
3. 选择最佳结构用于对接（推荐 9P45，分辨率 1.35Å）
4. 用 GROMACS 对 apo 蛋白做 100ns MD 模拟
5. 从 MD 轨迹中聚类提取 5–10 个代表性构象（ensemble docking 用）
6. 用 AlphaFold 3 补充预测缺失的 loop 区域
```

#### Step 2：已知抑制剂 Benchmark（2 周）

```
操作流程：
1. 从 ChEMBL 收集所有已知 KRAS G12D 抑制剂的结构和 IC50 数据
2. 包括：MRTX1133, RMC-9805, GFH375, HRS-4642, BI-2852 等
3. 进行 retrospective docking：已知活性分子 vs. decoy 分子
4. 验证你的 docking 流程能正确 rank 已知活性分子（AUC > 0.8）
5. 如果效果不好 → 调整 scoring function、尝试 Gnina（ML-based scoring）
6. 建立可靠的 scoring baseline → 这是后续所有筛选的基准
```

#### Step 3：AI 生成新分子（4 周）

```
操作流程：
1. 用 DiffSBDD 在 G12D Switch II pocket 中生成 ~10,000 分子
2. 设计自定义条件约束：
   - 必须：与 Asp12 形成离子/氢键相互作用
   - 必须：满足 Lipinski Rule of Five
   - 优选：分子量 300–500，LogP 1–4
3. 生成分子 → AutoDock Vina/Gnina 重新打分
4. Top 1000 → RDKit ADMET 过滤（溶解度、代谢稳定性、hERG 安全性）
5. Top 200 → 分子动力学验证（短程 10ns MD）
6. Top 50 → 深入分析

⚠️ 关键增加步骤（为了造出实物）：
7. Top 50 → AiZynthFinder / ASKCOS 合成可行性评估
8. 剔除无法合成或合成成本过高的分子
9. 最终得到 Top 20 "可合成候选分子"
```

#### Step 4：深度计算验证（4 周）

```
操作流程：
1. 对 Top 20 分子进行 200ns MD 模拟（每个分子）
2. 分析：
   - RMSD 稳定性（<2Å 为优）
   - 结合自由能（MM-PBSA/GBSA，越负越好）
   - 关键残基相互作用持续时间（Asp12, Thr35, Asp69 等）
   - 水分子桥接网络
3. 选出 Top 10 最稳定候选

4. 对 Top 10 进行耐药突变扫描：
   - 模拟 G13D, Q61H, A59G, Y64H 等已知/预测耐药突变
   - 用 FoldX / Rosetta / ESM-2 预测突变对结合力的影响
   - 选择对 ≥3 种突变仍保持 >50% 结合力的分子

5. 最终输出：Top 5 "可合成 + 高亲和力 + 耐药穿透" 候选分子
```

---

## 6. 第四阶段：化学合成——把分子造出来（9–18 个月）

> **这是从"计算论文"到"真正药物"的关键跨越。**

### 6.1 三条合成路径（按可行性排序）

#### 路径 1：CRO 外包合成（最快，推荐先行）⭐⭐⭐⭐⭐

**CRO（Contract Research Organization）= 合同研究组织**，专业做化学合成的公司。

| CRO 公司 | 特点 | 单个分子合成费用 | 周期 |
|---------|------|----------------|------|
| **WuXi AppTec（药明康德）** | 全球最大 CRO 之一，合成能力极强 | ¥5,000–30,000 | 2–6 周 |
| **Pharmaron（康龙化成）** | 中国头部 CRO | ¥5,000–25,000 | 2–6 周 |
| **Enamine** | 乌克兰/美国，价格有竞争力 | $1,000–5,000 | 3–8 周 |
| **Mcule** | 在线下单，小批量快速合成 | $500–3,000 | 2–4 周 |
| **Molport** | 化合物采购+定制合成 | 变化大 | 2–8 周 |

**操作流程**：
```
1. 将 Top 5 候选分子的 SMILES/结构图发给 CRO
2. CRO 评估合成可行性 → 报价（通常 1 周内）
3. 确认后支付定金 → CRO 开始合成
4. 2–6 周后收到产品（通常 5–50mg）
5. CRO 提供 NMR、MS、HPLC 纯度报告

💡 Pro Tip：
- 先合成 2–3 个最有信心的分子，不要一次全做
- 首次合成 5–10mg 即可，够做初步生物测试
- 保留所有分析报告，发论文需要
```

**预算估算**：5 个分子 × ¥15,000/个 = ¥75,000（~$10,000）

#### 路径 2：合作实验室合成（最省钱）⭐⭐⭐⭐

找有机化学/药物化学实验室合作：

**如何找合作者**：
- 你的博后导师推荐（最佳渠道）
- 同校/同城的药化课题组
- 学术会议上 networking
- 在线平台（如 ResearchGate）

**合作模式**：
```
你提供：
  ✓ 计算设计的分子结构
  ✓ 合成路线建议（AiZynthFinder 输出）
  ✓ 预测的活性数据作为合成优先级参考
  ✓ 共同署名

对方提供：
  ✓ 有机合成实验操作
  ✓ 纯化与表征（NMR, MS, HPLC）
  ✓ 合成经验（可能优化你的合成路线）
```

#### 路径 3：自己学合成（最慢但最有价值）⭐⭐⭐

如果所在实验室有合成条件：

**快速入门路线**（3 个月）：
```
第 1 月：有机化学实验基础
  - 简单反应操作（回流、萃取、旋蒸、柱层析）
  - 跟着有经验的人做 2–3 个已知反应

第 2 月：学习常见药化反应
  - Suzuki 偶联（最常用的 C-C 键形成反应）
  - 酰胺化反应（药物分子中最常见的键）
  - 还原胺化

第 3 月：尝试合成一个简单的候选分子
  - 选最简单的那个（步骤 ≤4 步）
  - 有人指导的情况下完成
```

### 6.2 ⚠️ 不要忽略的关键步骤

**合成可行性前置评估**（在计算阶段就要做）：

```python
# 使用 RDKit + AiZynthFinder 评估合成可行性的伪代码示意
from rdkit import Chem
from rdkit.Chem import Descriptors

def assess_synthesizability(smiles):
    mol = Chem.MolFromSmiles(smiles)
    
    # 1. Synthetic Accessibility Score (SA Score, 1-10, 越低越好)
    sa_score = calculate_sa_score(mol)  # 目标: < 4
    
    # 2. 分子量和复杂度
    mw = Descriptors.MolWt(mol)         # 目标: < 500
    
    # 3. 立体中心数量（手性中心越多越难合成）
    chiral_centers = count_chiral(mol)   # 目标: ≤ 2
    
    # 4. 是否含有已知的难合成基团
    problematic = check_problematic_groups(mol)
    
    # 5. AiZynthFinder 逆合成分析
    routes = aizynthfinder_analyze(smiles)  # 目标: 有 ≥1 条可行路线
    
    return {
        'sa_score': sa_score,
        'mw': mw,
        'chiral_centers': chiral_centers,
        'problematic_groups': problematic,
        'synthesis_routes': routes,
        'verdict': 'SYNTHESIZABLE' if sa_score < 4 and routes else 'DIFFICULT'
    }
```

### 6.3 候选分子的结构表征要求

合成出来后，必须确认分子是你设计的那个：

| 表征方法 | 确认内容 | 谁来做 |
|---------|---------|--------|
| ¹H NMR | 氢原子环境，确认结构 | CRO 提供 / 校内 NMR 中心 |
| ¹³C NMR | 碳骨架确认 | CRO 提供 / 校内 NMR 中心 |
| LC-MS | 分子量确认（±1 Da） | CRO 提供 / 校内质谱平台 |
| HPLC | 纯度（>95%） | CRO 提供 / 自己做 |
| HRMS | 高分辨质谱，精确分子式 | 发论文必需 |

---

## 7. 第五阶段：体外生物学验证（12–20 个月）

> **拿到合成的分子后，必须证明它真的能抑制 KRAS G12D。**

### 7.1 实验体系概览

```
体外验证金字塔（从简到繁）：

                    ┌─────────┐
                    │ 机制验证 │ ← Western Blot (p-ERK, KRAS降解)
                   ┌┴─────────┴┐
                   │ 细胞活性   │ ← IC50 (MTT/CCK-8)
                  ┌┴───────────┴┐
                  │ 蛋白结合实验 │ ← SPR / DSF / MST
                 ┌┴─────────────┴┐
                 │ 无细胞生化活性 │ ← KRAS-SOS1 结合抑制 / GTPase 活性
                ┌┴───────────────┴┐
                │   化合物溶解性    │ ← 动力学溶解度测定
                └─────────────────┘
```

### 7.2 具体实验操作

#### 实验 1：化合物溶解性测定（第 1 周）

```
目的：确认化合物在水溶液中可溶，否则后续实验都做不了

操作：
1. 配制 10mM DMSO 母液（最终 DMSO 浓度 < 0.5%）
2. 在 PBS buffer（pH 7.4）中梯度稀释
3. 37°C 孵育 24h
4. 离心/过滤后 UV 或 LC-MS 测定浓度
5. 溶解度 > 10 μM 即可继续

⚠️ 如果不溶：考虑加成盐、改侧链、或使用辅溶剂
```

#### 实验 2：蛋白水平结合实验（第 2–4 周）

**方法选择**：

| 方法 | 原理 | 需要的蛋白量 | 优先级 |
|------|------|------------|--------|
| **DSF（差示扫描荧光法）** | 配体结合使蛋白热稳定性提高 | 微量 | ⭐⭐⭐⭐⭐ 首选 |
| **MST（微量热泳动）** | 配体结合改变蛋白热泳动行为 | 微量 | ⭐⭐⭐⭐ |
| **SPR（表面等离子共振）** | 实时测量结合动力学 ka/kd/KD | 较多 | ⭐⭐⭐⭐⭐ 金标准 |
| **ITC（等温滴定量热）** | 热力学参数（ΔH, ΔS, KD） | 较多 | ⭐⭐⭐ |

**推荐：先做 DSF（最简单最快）→ 有结合信号再做 SPR（定量 KD）**

```
DSF 实验操作：
1. 购买/表达 KRAS G12D 重组蛋白（购买推荐：Abcam, R&D Systems, Addgene 质粒自己表达）
2. 蛋白 + SYPRO Orange 荧光染料 + 化合物（梯度浓度）
3. 用 qPCR 仪做热变性曲线（25→95°C，1°C/min）
4. 看 Tm shift：
   - Tm shift > 2°C → 有结合，继续
   - Tm shift < 1°C → 弱结合或无结合
5. 阳性对照：MRTX1133（可从 MedChemExpress 购买，~$200/mg）
```

```
SPR 实验操作（Biacore 或 OpenSPR）：
1. 将 KRAS G12D 蛋白固定在芯片上
2. 流过梯度浓度的化合物（通常 0.1–100 μM）
3. 实时记录结合/解离曲线
4. 拟合获得 KD 值
   - KD < 100 nM → 优秀
   - KD 100 nM – 1 μM → 良好
   - KD > 10 μM → 需优化
```

#### 实验 3：细胞活性测定（第 4–8 周）

**必需细胞系**：

| 细胞系 | KRAS 状态 | 癌种 | 用途 | 来源 |
|--------|----------|------|------|------|
| **MIA PaCa-2** | G12C | 胰腺癌 | 选择性对照 | ATCC |
| **PANC-1** | G12D | 胰腺癌 | **主要测试细胞** | ATCC |
| **AsPC-1** | G12D | 胰腺癌 | 验证细胞系 | ATCC |
| **LS180** | G12D | 结直肠癌 | 验证细胞系 | ATCC |
| **SW480** | G12V | 结直肠癌 | 选择性对照 | ATCC |
| **H358** | G12C | NSCLC | 选择性对照 | ATCC |
| **HEK293** | Wild-type | 正常细胞 | 毒性对照 | ATCC |

```
细胞增殖抑制实验（IC50 测定）：

Day 0：
  1. 将细胞接种到 96 孔板（3000–5000 cells/well）
  2. 每个细胞系 3 个生物学重复

Day 1：
  3. 加入化合物（10 个浓度梯度，通常 0.001–100 μM）
  4. 阳性对照：MRTX1133（购买）
  5. 阴性对照：DMSO（等体积）

Day 4（72h 后）：
  6. 加入 CCK-8 试剂（10 μL/well）
  7. 37°C 孵育 1–4h
  8. 酶标仪读 OD450
  9. GraphPad Prism 拟合剂量-反应曲线 → 得到 IC50

目标：
  - 在 G12D 细胞系（PANC-1, AsPC-1）IC50 < 1 μM
  - 在 non-G12D 细胞系 IC50 > 10 μM（证明选择性）
  - 在 HEK293 IC50 > 50 μM（低正常细胞毒性）
```

#### 实验 4：通路验证——Western Blot（第 6–10 周）

```
目的：证明你的分子真的抑制了 KRAS 下游信号

操作：
1. PANC-1 细胞接种于 6 孔板，生长至 70–80% 汇合度
2. 加化合物（IC50 浓度 + 2× IC50 + 5× IC50）处理 2h, 6h, 24h
3. 裂解细胞，收蛋白 → BCA 定量
4. SDS-PAGE 电泳 → 转膜（PVDF）
5. 抗体孵育：

   必检蛋白          抗体                   预期结果
   ──────────       ──────                 ──────────
   p-ERK1/2        Cell Signaling #4370    处理后 ↓↓↓
   total ERK1/2    Cell Signaling #4695    不变（loading 对照）
   p-AKT(S473)     Cell Signaling #4060    处理后 ↓↓
   total AKT       Cell Signaling #4691    不变
   KRAS            Abcam ab275876          如果是降解剂→↓↓↓
   β-actin         Cell Signaling #4967    不变（整体 loading 对照）

6. 化学发光检测 → ImageJ 定量
7. 与 MRTX1133 对照比较
```

#### 实验 5：KRAS G12D 选择性验证（第 8–12 周）

```
操作：
1. 在所有细胞系上重复 IC50 实验
2. 计算选择性指数：
   - SI = IC50(WT KRAS) / IC50(G12D)
   - 目标 SI > 10

3. 如果有条件，做 KRAS 过表达/敲除实验：
   - 在 G12D 细胞中 siRNA 敲低 KRAS → 化合物效果应消失
   - 在 WT 细胞中过表达 G12D → 化合物效果应出现
   - 这是证明"on-target"的最强证据
```

### 7.3 实验耗材预算

| 项目 | 预算 | 说明 |
|------|------|------|
| 细胞系购买（5-7 种） | ¥20,000–35,000 | ATCC/中科院细胞库 |
| 重组 KRAS G12D 蛋白 | ¥5,000–15,000 | 商购或自己表达 |
| 阳性对照化合物（MRTX1133 等） | ¥3,000–10,000 | MedChemExpress/Selleck |
| 抗体（6–8 种） | ¥15,000–25,000 | Cell Signaling Technology |
| 细胞培养耗材 | ¥10,000–20,000 | 培养基、血清、耗材 |
| CCK-8/MTT 试剂 | ¥2,000–5,000 | 多次实验 |
| Western Blot 耗材 | ¥5,000–10,000 | 凝胶、膜、ECL |
| SPR/DSF 耗材 | ¥5,000–15,000 | 芯片、试剂 |
| **体外实验总预算** | **¥65,000–135,000** | **约 $9,000–$18,000** |

---

## 8. 第六阶段：体内药效与药代动力学（18–28 个月）

> **从细胞到活体动物——证明你的分子在体内也有效。**

### 8.1 前提条件

开始体内实验前，候选分子必须满足：
- [x] 体外 IC50 < 1 μM（G12D 细胞系）
- [x] 选择性指数 > 10
- [x] 明确的 on-target 机制（p-ERK 抑制 / KRAS 降解）
- [x] 初步 ADMET 预测良好（计算 + 体外）
- [x] 所在机构已获得动物实验伦理批准（IACUC）

### 8.2 体外 ADMET 补充实验

在进动物之前，先做这些体外 ADMET：

| 实验 | 目的 | CRO 外包价格 | 自做可行性 |
|------|------|-------------|-----------|
| 肝微粒体代谢稳定性 | 预测肝脏清除率 | ¥3,000/化合物 | 需要微粒体 |
| Caco-2 渗透性 | 预测口服吸收 | ¥5,000/化合物 | 需要 Caco-2 细胞 |
| 血浆蛋白结合率 | 预测游离药物浓度 | ¥2,000/化合物 | 平衡透析法 |
| CYP 抑制 | 药物-药物相互作用 | ¥5,000/化合物 | CRO 推荐 |
| hERG 抑制 | 心脏安全性 | ¥8,000/化合物 | CRO 推荐 |
| 水溶性/LogD | 物化性质 | ¥1,000/化合物 | 可自做 |

**建议**：将体外 ADMET 打包给 CRO，通常有套餐价（¥15,000–30,000/化合物）。

### 8.3 药代动力学（PK）实验

```
目的：确认分子在动物体内能达到有效浓度

动物：CD-1 小鼠（非肿瘤鼠，只看 PK）
给药：
  - 口服（PO）：10–50 mg/kg，如果目标是口服药
  - 静脉注射（IV）：1–5 mg/kg，确定绝对生物利用度
  - 腹腔注射（IP）：10–30 mg/kg，备选给药途径

采样时间点：0.25, 0.5, 1, 2, 4, 8, 12, 24h（每时间点 3 只鼠）
检测：LC-MS/MS 测定血浆药物浓度

关键 PK 参数目标：
  - Cmax > 10× IC50（血浆峰浓度）
  - AUC：足够高（暴露量）
  - T1/2 > 2h（半衰期，越长越好）
  - F% > 20%（口服生物利用度，如果是口服药）
```

**PK 实验通常外包给动物 CRO**：
- 药明康德/康龙化成/昭衍新药：¥30,000–80,000/化合物（含分析）
- 学校动物中心 + 自行操作：¥10,000–30,000/化合物

### 8.4 体内药效实验（Xenograft 模型）

```
肿瘤异种移植模型：

动物：BALB/c nude 小鼠（裸鼠，免疫缺陷）
细胞：PANC-1 或 AsPC-1（KRAS G12D 胰腺癌细胞）

操作流程：
1. 皮下接种 5×10⁶ 肿瘤细胞（右侧腋下）
2. 等肿瘤长到 ~100 mm³ 后随机分组

实验分组（每组 6–8 只）：
  ┌──────────────────────┬──────────────────────┐
  │ 组别                  │ 处理                  │
  ├──────────────────────┼──────────────────────┤
  │ 溶剂对照组            │ Vehicle（溶剂）        │
  │ 低剂量组              │ 化合物 10 mg/kg QD    │
  │ 高剂量组              │ 化合物 30 mg/kg QD    │
  │ 阳性对照组            │ MRTX1133 30 mg/kg QD  │
  │ （可选）联合用药组     │ 化合物 + MEK 抑制剂   │
  └──────────────────────┴──────────────────────┘

测量指标：
  - 每 2–3 天用游标卡尺量肿瘤：V = (长 × 宽²) / 2
  - 每天称体重（安全性监测）
  - 给药 21–28 天

终点分析：
  - 肿瘤生长抑制率 TGI% = [1 - (T/C)] × 100%
    目标：TGI > 60% 为有效
  - 终末取肿瘤：
    a. 一半做 Western Blot（p-ERK, KRAS）
    b. 一半做 IHC（Ki67 增殖、Cleaved Caspase-3 凋亡）
  - 终末取血：测血药浓度（确认 PK/PD 关系）
```

### 8.5 体内实验预算

| 项目 | 预算 | 说明 |
|------|------|------|
| 裸鼠（40–50 只） | ¥15,000–25,000 | ~¥400/只 |
| 动物房费用 | ¥10,000–20,000 | 按天/笼计费 |
| PK 实验（2 个化合物） | ¥60,000–160,000 | 含动物+分析 |
| 药效实验耗材 | ¥10,000–20,000 | 注射器、量具等 |
| IHC 分析 | ¥5,000–10,000 | 外包或自做 |
| 化合物大量合成（100–500mg） | ¥20,000–50,000 | 体内实验需要更多量 |
| **体内实验总预算** | **¥120,000–285,000** | **约 $16,000–$39,000** |

### 8.6 如果没有动物实验条件

如果你的实验室没有动物房或 IACUC 批准：

1. **合作**：找有动物实验平台的合作课题组
2. **外包**：全流程外包给动物 CRO（药明康德体内药效评价套餐约 ¥150,000–300,000）
3. **3D 类器官替代**：用患者来源的胰腺癌类器官（organoid）做药效测试，不需要动物伦理
4. **推迟**：先把体外数据做扎实，体内作为后续合作的谈判筹码

---

## 9. 第七阶段：先导化合物优化闭环（贯穿全程）

> **这是整个项目的灵魂——计算与实验的迭代闭环。**

### 9.1 闭环迭代流程

```
                        ┌─────────────────────────┐
                        │    Round 1: 初始设计      │
                        │  AI 生成 → Top 5 候选分子  │
                        └────────────┬────────────┘
                                     ↓
                        ┌────────────────────────────┐
                        │    合成 → 生物测试           │
                        │  合成 3–5 个 → IC50, SPR    │
                        └────────────┬───────────────┘
                                     ↓
                      ┌──────────────────────────────┐
                      │    数据分析：什么有效什么没效？  │
                      │  活性分子的共同特征？           │
                      │  无活性分子失败的原因？          │
                      └──────────────┬───────────────┘
                                     ↓
                        ┌────────────────────────────┐
                        │    更新计算模型              │
                        │  · 将实验 IC50 加入训练集     │
                        │  · Fine-tune 生成模型        │
                        │  · 更新 QSAR/scoring 模型    │
                        └────────────┬───────────────┘
                                     ↓
                        ┌────────────────────────────┐
                        │    Round 2: 优化设计         │
                        │  基于 SAR 生成改进分子        │
                        │  聚焦活性骨架的衍生物         │
                        └────────────┬───────────────┘
                                     ↓
                              （重复 3–5 轮）
                                     ↓
                        ┌────────────────────────────┐
                        │    最终候选分子              │
                        │  IC50 < 100nM              │
                        │  选择性 > 50 倍             │
                        │  良好 ADMET + PK            │
                        └────────────────────────────┘
```

### 9.2 每轮迭代的具体操作

**Round 1（月 9–14）**：初始验证
- 合成 3–5 个计算设计的分子
- 测 DSF/SPR + IC50
- 获得首批实验数据

**Round 2（月 14–18）**：SAR 驱动优化
- 分析 Round 1 的 SAR（哪些基团重要？）
- AI 模型融入实验数据，生成新一批分子
- 合成 5–10 个优化分子
- 测全套体外活性

**Round 3（月 18–22）**：ADMET 优化
- 聚焦代谢稳定性和溶解度问题
- 通过引入代谢阻断基团（如氟取代、环丙基）改善 PK
- 合成 5–10 个 PK 优化分子
- 做体外 ADMET + 初步 PK

**Round 4（月 22–28）**：最终优化
- 综合活性、选择性、PK 数据
- 选择 1–2 个最佳候选分子进入体内
- 做 xenograft 药效实验

### 9.3 关键：让实验数据反哺计算模型

```python
# 主动学习闭环（Active Learning Loop）伪代码

class DrugDesignLoop:
    def __init__(self):
        self.model = GenerativeModel()       # DiffSBDD 或自定义
        self.scorer = AffinityPredictor()     # GNN-based
        self.experimental_data = []           # 实验结果存储
    
    def run_iteration(self, round_num):
        # 1. 生成候选分子
        candidates = self.model.generate(n=1000, pocket=KRAS_G12D_pocket)
        
        # 2. 计算打分 + 合成可行性过滤
        scored = self.scorer.predict(candidates)
        synthesizable = filter_by_sa_score(scored, threshold=4.0)
        
        # 3. 选择最有信息量的分子（不只是预测最好的）
        # 使用 uncertainty sampling：选预测不确定性最高的
        selected = self.acquisition_function(
            synthesizable, 
            strategy='uncertainty',  # 或 'expected_improvement'
            n_select=5
        )
        
        # 4. 合成 + 测试（实验室操作）
        experimental_results = wet_lab_testing(selected)
        # 返回: {smiles: IC50, binding_KD, selectivity, ...}
        
        # 5. 更新模型
        self.experimental_data.extend(experimental_results)
        self.scorer.fine_tune(self.experimental_data)
        self.model.update_reward(self.experimental_data)
        
        # 6. 分析本轮 SAR
        sar_insights = analyze_sar(self.experimental_data)
        
        return experimental_results, sar_insights
```

---

## 10. 第八阶段：成果产出、专利与产业化（24+ 个月）

### 10.1 论文发表策略

| 时间 | 产出 | 内容 | 目标期刊 |
|------|------|------|---------|
| 6–9 月 | 方法论文 | 新的生成模型 or 耐药预测框架 | J. Chem. Inf. Model. / Briefings in Bioinformatics |
| 12–15 月 | 计算+体外论文 | AI 设计 + DSF/SPR 结合验证 | J. Med. Chem. / ACS Med. Chem. Lett. |
| 18–22 月 | 完整体外论文 | 计算设计 + IC50 + WB + 选择性 + SAR | Nature Communications / Cell Chemical Biology |
| 24–28 月 | 体内验证论文 | 全部数据 + xenograft 药效 + PK | Nature Chemical Biology / Cancer Discovery |
| 30+ 月 | 临床前候选论文 | 优化分子 + 全面临床前数据 | Nature / Science |

### 10.2 专利策略

> **重要：在发表论文之前必须提交专利申请！**

```
专利时间线：

月 12：首批有活性的分子确认后 → 咨询专利律师
月 14：准备临时专利申请（Provisional Patent）→ 保护优先权
月 15：提交美国临时专利申请（成本 ~$2,000–5,000）
月 16：此时才能安全地发预印本或投稿
月 26：补充数据后转为 PCT 国际专利申请（成本 ~$10,000–20,000）

专利覆盖范围：
  ✓ 化合物结构（Markush 式覆盖系列化合物）
  ✓ 治疗 KRAS G12D 突变癌症的方法
  ✓ 化合物的合成方法
  ✓ AI 设计方法（如果有创新性）
```

### 10.3 预印本与开源

- 第一时间在 bioRxiv/chemRxiv 发布预印本（专利提交后！）
- 开源代码和模型（GitHub），增加引用和影响力
- 开源分子数据集供社区验证
- 开源计算工作流（增强方法论文的影响力）

### 10.4 产业化路径

```
学术成果 → 产业转化路径：

路径 1：License Deal（最常见）
  有活性数据 → 找药企谈 license
  里程碑付款 + 销售版税
  你继续做学术，药企负责后续开发

路径 2：创业
  组建团队 → Seed 融资（$1–3M）
  → 做完临床前全部数据
  → Series A（$10–30M）
  → 进入 IND 申报和 Phase I
  参考：Insilico Medicine 从 AI 到 Phase II 仅用约 3 年

路径 3：学术药物发现中心
  很多大学有 Drug Discovery Center
  → 申请内部 pilot grant（$50,000–200,000）
  → 做到 clinical candidate
  → License 或 spin-off
```

---

## 11. 实验室搭建与资源规划

### 11.1 你需要的硬件资源

**计算资源**：

| 资源 | 用途 | 获取方式 |
|------|------|---------|
| GPU 服务器（至少 4× A100） | 训练生成模型、MD 模拟 | 校内 HPC / 云计算（AWS/GCP） |
| CPU 集群 | 大规模对接、ADMET 筛选 | 校内 HPC |
| 存储（>10TB） | MD 轨迹、分子库 | 校内 / AWS S3 |

**实验资源**（确认你的实验室 / 合作者有）：

| 仪器 | 用途 | 是否必须 |
|------|------|---------|
| 细胞培养间（CO2 培养箱、生物安全柜） | 所有细胞实验 | ✅ 必须 |
| 酶标仪 | 读 IC50 板 | ✅ 必须 |
| Western Blot 全套 | 通路验证 | ✅ 必须 |
| qPCR 仪 | DSF 结合实验 | ⭐ 强烈推荐 |
| SPR（如 Biacore） | 结合亲和力定量 | ⭐ 推荐（可用校内共享平台） |
| LC-MS | 化合物分析 | ⭐ 推荐（可外送） |
| NMR | 结构确认 | 外送/CRO |
| 动物房 | 体内实验 | 后期需要（可合作/外包） |

### 11.2 总体经费预算

| 阶段 | 预算 | 说明 |
|------|------|------|
| 计算资源 | ¥50,000–100,000/年 | GPU 时间、云计算 |
| 化学合成（3 轮 × 5–10 个分子） | ¥150,000–400,000 | CRO + 合作 |
| 体外生物实验 | ¥65,000–135,000 | 细胞、抗体、试剂 |
| 体外 ADMET | ¥30,000–60,000 | CRO 外包 |
| PK 实验 | ¥60,000–160,000 | CRO 外包 |
| 体内药效 | ¥60,000–125,000 | 动物+耗材 |
| 杂项（差旅、发表费） | ¥30,000–50,000 | |
| **总预算（2–3 年）** | **¥445,000–1,030,000** | **约 $60,000–$140,000** |

### 11.3 经费来源

| 来源 | 金额 | 说明 |
|------|------|------|
| 博后导师课题经费 | 变化大 | 首选，最快获得 |
| NSFC 青年项目（中国） | ¥300,000 | 竞争激烈但值得申请 |
| NIH K99/R00（美国） | $250,000+ | 博后到独立 PI 的桥梁基金 |
| NIH R21（美国） | $275,000 | 探索性/高风险项目 |
| 校内 seed grant | $10,000–50,000 | 很多大学有，用于启动 |
| 产业合作（药企） | 变化大 | 用初步数据谈 |
| AI 药物发现专项 | 变化大 | 近年各国都有 |

---

## 12. 合作网络构建策略

> **一个人做不完所有事。构建合作网络是成功的关键。**

### 12.1 你需要的合作者

| 角色 | 为什么需要 | 如何找到 | 合作模式 |
|------|----------|---------|---------|
| **药物化学家** | 合成指导 + 分子优化 SAR 经验 | 同校药化系/会议 | 共同通讯作者 |
| **有机合成博士** | 实际合成操作 | 实验室内/同系 | 共同第一作者 |
| **癌症生物学家** | 细胞/动物模型 + 机制解读 | 导师推荐/同楼 | 共同通讯作者 |
| **结构生物学家** | 解析化合物-蛋白共晶结构 | 同校/共享平台 | 共同作者 |
| **临床医生（肿瘤科）** | 临床意义解读 + 未来临床试验 | 附属医院 | 顾问/共同作者 |
| **专利律师** | 知识产权保护 | 学校技术转让办公室 | 顾问 |

### 12.2 说服合作者的策略

```
你带来的价值（让合作者无法拒绝）：

1. 计算预测的可信度
   → 展示已知抑制剂的 retrospective 验证结果
   → "我的模型能正确预测 MRTX1133 的活性，AUC = 0.92"

2. 清晰的合成路线
   → AiZynthFinder 输出的逆合成路线
   → "这个分子只需要 3 步合成，起始原料都是市售的"

3. 明确的实验计划
   → 本文档就是你的武器
   → 展示你对实验细节的理解，不是"请帮我做实验"

4. 经费支持
   → 如果你能带来一部分经费（课题、grant），合作更容易

5. 共同署名 + 共同基金申请
   → 联合申请 NSFC/NIH，增加双方获批几率
```

### 12.3 网络建设时间线

```
月 1–3：内部网络
  → 认识本实验室/系的所有人
  → 找到 1-2 个可以教你实验的人
  → 与导师讨论合作策略

月 3–6：校内网络
  → 拜访药化系、肿瘤生物学系的 PI
  → 参加校内研讨会、报告你的计划
  → 锁定 1–2 个核心合作者

月 6–12：外部网络
  → 参加 1–2 个学术会议（AACR, ACS Medicinal Chemistry, CASP）
  → 在会议上做海报或口头报告
  → 联系领域内的通讯作者（通过邮件 cold email）

月 12+：产业网络
  → 参加 biotech 投资会议
  → 联系 CRO 建立长期合作
  → 与感兴趣的药企 BD 部门对接
```

---

## 13. 风险管理与 Plan B

### 13.1 主要风险与应对

| 风险 | 概率 | 影响 | 应对策略 |
|------|------|------|---------|
| AI 生成的分子全部无活性 | 中 | 高 | Plan B1: 转向已知活性骨架的 SAR 优化 |
| 合成困难/失败 | 中 | 中 | 选择高合成可行性分子；多个 CRO 备选 |
| 找不到实验合作者 | 低–中 | 高 | Plan B2: 全部外包给 CRO（更贵但可行） |
| 体外有活性但体内无效 | 中 | 高 | PK 优化（前药策略、代谢阻断基团） |
| 经费不足 | 中 | 高 | 分阶段申请；先用最小 budget 证明概念 |
| 选择性不够（毒性大） | 中 | 中 | 聚焦突变选择性设计；扩大选择性窗口 |
| 被竞争对手抢先 | 高 | 中 | 聚焦差异化（耐药穿透是你的独特角度） |

### 13.2 Plan B 详解

**Plan B1：如果 AI 生成分子全部失败**
```
1. 回退到 MRTX1133 的化学骨架
2. 用 AI 做 scaffold hopping（骨架跃迁）
3. 在已知活性骨架上做衍生物设计
4. 这种方式成功概率更高，但创新性降低
5. 仍然可以发论文（SAR 研究 + AI 辅助优化）
```

**Plan B2：如果找不到合作者**
```
全 CRO 外包路线：
  - 合成：药明康德/Enamine（¥10,000–30,000/分子）
  - 体外活性：Eurofins/Crown Bioscience（IC50 套餐 ¥20,000–50,000）
  - ADMET：Cyprotex/Eurofins（¥15,000–30,000/化合物）
  - PK：药明康德/昭衍（¥30,000–80,000/化合物）
  - 体内药效：Crown Bioscience/药明（¥100,000–200,000）
  
  总计：¥175,000–390,000/化合物（约 $24,000–$53,000）
  
  优点：你完全控制进度和数据所有权
  缺点：费用高，没有人面对面指导实验细节
```

**Plan B3：转向计算为主的方向**
```
如果实验确实完全不可行：
  1. 做纯计算 + 公开数据验证（用文献中的实验数据验证你的预测）
  2. 开发开源工具/数据库（如 KRAS 突变-药物数据库）
  3. 与有实验数据的团队合作做计算-实验联合分析
  4. 目标：方法论文 + 工具论文（影响力可能不如带实验数据的）
```

---

## 14. 关键成功因素

### 14.1 必须做到的 10 件事

1. **从第一天起就进实验室** —— 不要等计算做完再碰实验
2. **尽早锁定合作者** —— 第 3 个月前必须有至少一个实验合作方
3. **选择可合成的分子** —— 在计算阶段就用合成可行性过滤
4. **小步快跑** —— 每轮只合成 3–5 个分子，快速迭代
5. **让数据说话** —— 每个 claim 都要有实验证据支撑
6. **先专利后发表** —— 有活性分子后第一件事是联系专利律师
7. **聚焦差异化** —— "AI 设计 + 耐药穿透"是你的独特角度
8. **理解生物学** —— 计算再牛也要 make biological sense
9. **善用 CRO** —— 不会的实验外包，不要浪费时间自己摸索
10. **保持与临床的连接** —— 关注 Zoldonrasib/GFH375 的临床数据，耐药数据是你的下一个机会

### 14.2 时间线总览

```
月份    计算工作                        实验工作                        里程碑
────    ────────                       ────────                       ──────
0–3     文献调研，学习工具                学习细胞培养、Western Blot      ✓ 知识体系建立
3–6     确定方向，搭建流程                独立做细胞实验，找合作者         ✓ 方向确定
6–9     靶点准备，benchmark               继续练手，准备细胞系和试剂      ✓ 计算流程验证
9–12    AI 生成+筛选 Top 5               送合成，准备生物测试体系         ✓ 首批分子设计完成
12–15   分析实验数据，更新模型             DSF/SPR + IC50（Round 1）     📄 方法论文投稿
15–18   Round 2 设计优化                  Round 2 合成+测试             ✓ 活性优化中
18–22   Round 3 ADMET 优化               体外 ADMET + 初步 PK          📄 体外验证论文投稿
22–26   最终分子选择                      体内 xenograft 药效            ⚖️ 专利申请
26–30   总结分析                          补充实验                      📄 高影响力论文投稿
30+     扩展到其他突变                    临床前 IND-enabling            🎯 候选分子确定
```

### 14.3 最终目标

在 2–3 年博士后期间，实现以下具体产出：

- [x] **1 个以上** 具有体外活性（IC50 < 1μM）的 KRAS G12D 新化学实体
- [x] **完整的 SAR 数据** 覆盖 ≥20 个化合物
- [x] **体外机制验证** 证明 on-target（p-ERK 抑制/KRAS 降解）
- [x] **初步体内药效** 在 xenograft 模型中 TGI > 60%
- [x] **3–5 篇论文**（至少 1 篇 IF > 10 的杂志）
- [x] **≥1 项专利申请**
- [x] **开源代码和数据集**
- [x] **建立 "AI + 靶向药" 方向的学术声誉**

---

## 15. 参考资料

### 关键论文与综述

1. Wang X, et al. "Targeting the untargetable: accelerated discovery of KRAS G12D inhibitors." *ScienceDirect*, 2025. https://www.sciencedirect.com/science/article/pii/S3050787125000988
2. Li Y, et al. "New exploration of KRAS G12D inhibitors and the mechanisms of resistance." *PMC*, 2025. https://pmc.ncbi.nlm.nih.gov/articles/PMC11912584/
3. "Targeting KRAS G12D: Advances in Inhibitor Design." *PMC*, 2025. https://pmc.ncbi.nlm.nih.gov/articles/PMC12719395/
4. "Discovery of KRAS (G12D) selective degrader ASP3082." *Nature Chemical Biology*, 2025. https://www.nature.com/articles/s42004-025-01662-4
5. "Quantum biological convergence: quantum computing accelerates KRAS drug discovery." *Nature Signal Transduction and Targeted Therapy*, 2025. https://www.nature.com/articles/s41392-025-02239-2
6. "Machine Learning-Driven Drug Repurposing for KRAS G12C and KRAS G12D." *bioRxiv*, 2025. https://www.biorxiv.org/content/10.1101/2025.05.16.654410v1
7. "Optimizing drug design by merging generative AI with a physics-based approach." *PMC*, 2025. https://pmc.ncbi.nlm.nih.gov/articles/PMC12334747/
8. "Computational drug design in the artificial intelligence era." *Pharmacological Reviews*, 2025. https://pharmrev.aspetjournals.org/article/S0031-6997(25)07503-9/fulltext
9. "Computational Discovery of Selective KRAS G12D Inhibitors Using a Consensus Docking Pipeline." *ChemRxiv*. https://chemrxiv.org/doi/pdf/10.26434/chemrxiv.15000641
10. "Discovery of Novel Noncovalent KRAS G12D Inhibitors through Computational Methods." *Molecules*, 2024. https://www.mdpi.com/1420-3049/29/6/1229

### 临床试验数据

11. "Oral Investigational Agent Zoldonrasib Elicits Objective Responses in KRAS G12D-Mutated Lung Cancer." *AACR 2025*. https://www.aacr.org/about-the-aacr/newsroom/news-releases/oral-investigational-agent-zoldonrasib-elicits-objective-responses-in-patients-with-kras-g12d-mutated-lung-cancer/
12. "First-in-Class Investigational Drug KRAS G12D Degrader Shows Promise in Lung and Pancreatic Cancer." *MSK*, 2026. https://www.mskcc.org/news/first-in-class-investigational-drug-kras-g12d-degrader-shows-promise-in-lung-and-pancreatic
13. "GenFleet Therapeutics Announces Phase I Data from Clinical Study of GFH375." *GenFleet*, 2025. http://www.genfleet.com/en/press_release-91
14. "ESMO 2025 – Incyte impresses in KRAS G12D." *ApexOnco*, 2025. https://www.oncologypipeline.com/apexonco/esmo-2025-incyte-impresses-kras-g12d
15. "Arvinas Presents Preclinical Data for ARV-806." *BioSpace*, 2025. https://www.biospace.com/press-releases/arvinas-presents-preclinical-data-for-arv-806-demonstrating-robust-and-differentiated-activity-in-models-of-kras-g12d-mutated-cancer-at-the-2025-aacr-nci-eortc-international-conference-on-molecular-targets-and-cancer-therapeutics

### 计算工具与数据库

16. DiffSBDD: Structure-based Drug Design with Equivariant Diffusion Models. https://github.com/arneschneuing/DiffSBDD
17. AiZynthFinder (AstraZeneca): Retrosynthesis Planning. https://github.com/MolecularAI/aizynthfinder
18. ASKCOS (MIT): Computer-Aided Synthesis Planning. https://askcos.mit.edu
19. RDKit: Cheminformatics Toolkit. https://www.rdkit.org
20. GROMACS: Molecular Dynamics Simulation. https://www.gromacs.org
21. AutoDock Vina: Molecular Docking. https://vina.scripps.edu
22. Gnina: ML-based Molecular Docking. https://github.com/gnina/gnina
23. ESM-2/3 (Meta): Protein Language Models. https://github.com/facebookresearch/esm

### 蛋白质结构

24. PDB 7RPZ: KRAS G12D + MRTX1133. https://www.rcsb.org/structure/7RPZ
25. PDB 9P45: KRAS G12D + BBO-11818 (1.35Å). https://www.rcsb.org/structure/9P45
26. PDB 9L6F: KRAS G12D + ASP3082. https://www.rcsb.org/structure/9L6F
27. PDB 7F0W: KRAS G12D + GDP (Switch I open). https://www.rcsb.org/structure/7F0W
28. PDB 9C3K: KRAS G12D/M67R + GDP. https://www.rcsb.org/structure/9C3K

### 化合物采购与 CRO

29. MedChemExpress (MCE): https://www.medchemexpress.com （阳性对照化合物）
30. Enamine: https://enamine.net （化合物合成与采购）
31. WuXi AppTec: https://www.wuxiapptec.com （全流程 CRO）
32. Crown Bioscience: https://www.crownbio.com （体内药效评价）

---

> **最后的话**：KRAS G12D 正从"不可成药"转向"可成药但尚未完全解决"。计算方法在其中的角色越来越核心，而你的计算机背景正是最大的武器。但请记住——**目标不是发表另一篇计算论文，而是真正设计、合成、验证出一个有效的分子**。这需要你走出舒适区，拥抱实验，构建合作，迭代优化。Timing 极好，机会窗口正在打开。去做吧。
