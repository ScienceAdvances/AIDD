# 06 · 受约束的序列优化

区分硬约束与优化目标，用共同八组解码顺序和 Pareto 选择解释提名。

[阅读教案](tutorial.md) · [回到总目录](../../README.md) · [结果解释](../../docs/EXPECTED_RESULTS.md)

## 输入与对应代码

examples/debian13_20261023/sequence_debian_only/ 的 candidates.csv、score_audit.npz 和 shortlist.csv。

所有命令均在仓库根目录执行。共享代码位于 `scripts/`，GNINA 的任务与辅助程序保存在 `workflows/gnina/`；本章无需复制另一份脚本。

## 先复算已保存结果

```bash
python practice.py analyze 06 --output runs/06_analysis
```

此入口读取仓库中的 Debian 示例并重新计算分析，不会启动 GPU 设计模型。

验收时应看到：197 条候选 → 176 条合格 → 5 条 Pareto → WT/D5E/S4E/L9K 四条提名；评分分母为十个可变位点。

## 再执行新的计算

先按[环境说明](../../docs/ENVIRONMENTS.md)配置工具。下面第一条只显示命令计划，确认路径与资源后再执行第二条；输出目录必须尚不存在。

```bash
python practice.py run 06 --output runs/06_new
python practice.py run 06 --output runs/06_new --execute
```

重新推理前先执行 `python practice.py fetch proteinmpnn_weights`。第 06 章默认对随仓库提供的六条 Debian 设计评分；若要换成自己的六条设计，直接调用 `scripts/optimize_sequence.py --design-run runs/02_new/design --output runs/06_own --device cuda`。

## 独立练习与交付

重新计算逐位 NLL 的平均、硬约束和 Pareto 集合。解释最低 NLL 的六突变样本为什么落选，以及 L9R 为何没有进入本轮四条提名。

保存实际输入、版本、命令、结果和一段解释。新计算的 `execution.json` 核验退出码与产物，候选是否达到研究目标还要继续读分数、几何和过滤标志。
