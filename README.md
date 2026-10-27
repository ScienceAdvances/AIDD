# AIDD 教程配套代码与练习数据

这是 timholy.com AIDD 教程的配套练习仓库，按 **00–15 章**组织。示例采用 Debian 13 / RTX 3060 的实际计算，以及公开 CAR-T 数据的重新分析。先用保存的输入与输出练习指标复算，再按需要配置软件、执行新计算。

## 从这里开始

下载或克隆后，在仓库根目录执行。若已经有分析环境，直接使用它；否则可自行创建读者专用环境：

```bash
python3 -m venv .venv-analysis
source .venv-analysis/bin/activate
python -m pip install -r env/requirements-analysis.txt
python practice.py doctor
python practice.py list
python practice.py verify
```

先从 CAR-T 入手，可以从原始 Excel 一直复现到图：

```bash
python practice.py analyze 04 --output runs/04_first
```

准备好基础分析依赖后，可以一次复算所有有计算内容的章节：

```bash
python practice.py analyze all --output runs/first_pass
```

这个入口会重算已保存的几何、分数、过滤与描述统计；第 09 章会生成 RDKit 构象。它不启动结构设计模型、不安装软件或下载模型，不要求 GPU。`runs/` 中保留自己的输出，`examples/` 保留已核验的示例。输出名称已存在时换一个新名字。

## 按章节练习

| 章 | 内容与练习入口 | 练习方式 |
|---|---|---|
| 00 | [学习主线与练习导航](chapters/00_start/README.md) | 阅读与讲述 |
| 01 | [多肽结构与 HADDOCK3](chapters/01_haddock/README.md) | 离线复算；可选新计算 |
| 02 | [ProteinMPNN 序列设计](chapters/02_proteinmpnn/README.md) | 离线复算；可选新计算 |
| 03 | [GROMACS 分子动力学](chapters/03_gromacs/README.md) | 离线复算；可选新计算 |
| 04 | [CAR-T 原始数据与工程判断](chapters/04_cart/README.md) | 离线复算；可选新计算 |
| 05 | [跨案例证据与综合简报](chapters/05_integration/README.md) | 证据简报 |
| 06 | [受约束的序列优化](chapters/06_sequence_optimization/README.md) | 离线复算；可选新计算 |
| 07 | [GNINA 小分子重对接](chapters/07_gnina/README.md) | 离线复算；可选新计算 |
| 08 | [面试问答与研究讲述](chapters/08_interview/README.md) | 阅读与讲述 |
| 09 | [RDKit 化学身份与构象](chapters/09_rdkit/README.md) | 离线复算；可选新计算 |
| 10 | [PyRosetta 侧链精修与界面](chapters/10_pyrosetta/README.md) | 离线复算；可选新计算 |
| 11 | [RFdiffusion2 原子基序支架](chapters/11_rfdiffusion2/README.md) | 离线复算；可选新计算 |
| 12 | [RFantibody 抗体与 VHH](chapters/12_rfantibody/README.md) | 离线复算；可选新计算 |
| 13 | [AlphaFold 3 复合物预测](chapters/13_alphafold3/README.md) | 离线复算；可选新计算 |
| 14 | [Boltz-2 结构与亲和力输出](chapters/14_boltz2/README.md) | 离线复算；可选新计算 |
| 15 | [BoltzGen 结合分子设计](chapters/15_boltzgen/README.md) | 离线复算；可选新计算 |

建议先走 01→02→06→10 的多肽线，再读 03 的独立 MD 案例和 04 的 CAR-T 数据；09→07 是小分子线，11–15 为生成与预测进阶。每章的 README 给出输入、代码、命令、预期结果和独立练习，tutorial.md 提供文章式讲解。

## 目录怎样对应

```text
AIDD-practice/
├── practice.py             # list / doctor / verify / analyze / fetch / run
├── chapters/00_...15_.../   # 每章 README 和教程正文
├── scripts/                # 全章共享的运行、分析和验收代码
├── data/                   # 小型原始输入、工作簿、固定版 ProteinMPNN 源码
├── workflows/              # 对接任务、约束、模型输入 JSON/YAML/配置
├── examples/               # Debian 参考结果、结构、阶段表和小型轨迹
├── config/                 # 下载清单和工具路径示例
├── env/                    # 读者基础依赖与原 Debian 依赖记录
├── docs/                   # 环境、数据、复现、预期结果与上传说明
├── provenance/             # 来源、校验与分发清单
├── third_party/            # 第三方来源及许可证
└── runs/                   # 自己产生的结果，已忽略
```

## 需要 GPU 的新计算

先按[环境说明](docs/ENVIRONMENTS.md)安装或复用工具，然后复制工具配置：

```bash
cp config/tools.example.json config/tools.json
# 编辑 config/tools.json，填写自己的解释器、软件和模型资源位置
python practice.py run 02 --output runs/02_new
# 上一条只预览命令；准备好后执行
python practice.py fetch proteinmpnn_weights
python practice.py run 02 --output runs/02_new --execute
```

每章使用自己的工具环境，配置文件不进入 Git。执行器记录命令、工作目录、日志、返回码和预期产物；科学结论继续依据几何、分数、过滤标志及实验终点判断。

## 数据、版本与发布

基础仓库包含小型输入、阶段表、代表结构和可复算示例。模型权重、完整含水轨迹和大规模中间结构采用按需下载或可选 Release 数据包，参见[数据说明](docs/DATA.md)。数据来源与第三方许可见 [THIRD_PARTY_NOTICES](THIRD_PARTY_NOTICES.md)。

[预期结果](docs/EXPECTED_RESULTS.md)说明每项数字的含义；[复现说明](docs/REPRODUCIBILITY.md)说明版本、随机性和计算边界；[上传说明](docs/PUBLISHING.md)提供把这个目录独立上传到 GitHub 的步骤。ESMFold2 保持待办状态。
