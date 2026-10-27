# 15 · BoltzGen 结合分子设计

沿当前 0.3.2 五阶段流程读取真实过滤标志，解释有排名却零通过的结果。

[阅读教案](tutorial.md) · [回到总目录](../../README.md) · [结果解释](../../docs/EXPECTED_RESULTS.md)

## 输入与对应代码

workflows/modern_tools/boltzgen/；examples/debian13_20261023/boltzgen_gpu/final_ranked_designs/。

所有命令均在仓库根目录执行。共享代码位于 `scripts/`，GNINA 的任务与辅助程序保存在 `workflows/gnina/`；本章无需复制另一份脚本。

## 先复算已保存结果

```bash
python practice.py analyze 15 --output runs/15_analysis
```

此入口读取仓库中的 Debian 示例并重新计算分析，不会启动 GPU 设计模型。

验收时应看到：2 条候选，0 条通过全部筛选；整体 filter RMSD 约 22.46–24.17 Å，自身 binder RMSD 约 1.33–1.39 Å。ESMFold2 未执行。

## 再执行新的计算

先按[环境说明](../../docs/ENVIRONMENTS.md)配置工具。下面第一条只显示命令计划，确认路径与资源后再执行第二条；输出目录必须尚不存在。

```bash
python practice.py run 15 --output runs/15_new
python practice.py run 15 --output runs/15_new --execute
```

入口按 BoltzGen 0.3.2 的 configure/execute 生成五阶段任务。它不执行 ESMFold2；不要用其他版本的阶段名称替代实际 steps.yaml。

## 独立练习与交付

对照 YAML 理解 exclude 和插入：CDR1 保留原环 5 位再加 1..5，CDR2 保留 4 位再加 1..5，CDR3 保留 8 位再加 1..12，得到长度范围 6..10、5..9、9..20。再解释排名第一为什么仍不合格。

保存实际输入、版本、命令、结果和一段解释。新计算的 `execution.json` 核验退出码与产物，候选是否达到研究目标还要继续读分数、几何和过滤标志。
