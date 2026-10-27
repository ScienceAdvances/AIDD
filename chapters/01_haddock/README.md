# 01 · 多肽结构与 HADDOCK3

区分生成好姿态、给它高排名与聚类选择将其保留这三个环节。

[阅读教案](tutorial.md) · [回到总目录](../../README.md) · [结果解释](../../docs/EXPECTED_RESULTS.md)

## 输入与对应代码

data/peptide/HADDOCK3-protein-peptide/；examples/haddock_stages/；examples/debian13_20261023/haddock_analysis/。

所有命令均在仓库根目录执行。共享代码位于 `scripts/`，GNINA 的任务与辅助程序保存在 `workflows/gnina/`；本章无需复制另一份脚本。

## 先复算已保存结果

```bash
python practice.py analyze 01 --output runs/01_analysis
```

此入口读取仓库中的 Debian 示例并重新计算分析，不会启动 GPU 设计模型。

验收时应看到：smoke 的 EM 最佳 DockQ 约 0.569，最终集合最佳约 0.055；full 最终保留 32 个模型，最佳及评分第一 DockQ 均约 0.569。

## 再执行新的计算

先按[环境说明](../../docs/ENVIRONMENTS.md)配置工具。下面第一条只显示命令计划，确认路径与资源后再执行第二条；输出目录必须尚不存在。

```bash
python practice.py run 01 --output runs/01_new
python practice.py run 01 --output runs/01_new --execute
```

入口运行 12 模型 smoke。完成后可进入 `runs/01_new/practice/`，用自己的 HADDOCK3 执行 `haddock3 workflows/full.cfg`；full 的默认 4 个 CPU 核可按资源调整。本轮教材作者曾用 32 核，改变核数应保存实际配置。

## 独立练习与交付

从逐阶段 TSV 复算 DockQ，指出 full 柔性精修阶段评分第一 0.090 与集合最佳 0.566 的差异。再查看 AIR 中受体 active 与肽 passive 的作用。

保存实际输入、版本、命令、结果和一段解释。新计算的 `execution.json` 核验退出码与产物，候选是否达到研究目标还要继续读分数、几何和过滤标志。
