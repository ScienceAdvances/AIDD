# 07 · GNINA 小分子重对接

在共同受体坐标系中评价小分子姿势，并区分经验能量、CNNscore 与 CNNaffinity。

[阅读教案](tutorial.md) · [回到总目录](../../README.md) · [结果解释](../../docs/EXPECTED_RESULTS.md)

## 输入与对应代码

workflows/gnina/example_1iep/；workflows/gnina/jobs_1iep.tsv；examples/debian13_20261023/gnina_*_gpu/。

所有命令均在仓库根目录执行。共享代码位于 `scripts/`，GNINA 的任务与辅助程序保存在 `workflows/gnina/`；本章无需复制另一份脚本。

## 先复算已保存结果

```bash
python practice.py analyze 07 --output runs/07_analysis
```

此入口读取仓库中的 Debian 示例并重新计算分析，不会启动 GPU 设计模型。

验收时应看到：三种子 × 两模式共六任务，每次 20 个姿势；首位重原子 RMSD 为 0.786–0.793 Å，两种模式首位相同。

## 再执行新的计算

先按[环境说明](../../docs/ENVIRONMENTS.md)配置工具。下面第一条只显示命令计划，确认路径与资源后再执行第二条；输出目录必须尚不存在。

```bash
python practice.py run 07 --output runs/07_new
python practice.py run 07 --output runs/07_new --execute
```

## 独立练习与交付

从输出 SDF 重新提取属性并计算 RMSD；说明为什么不能先单独拟合配体，也不能由本例宣称 CNN 更优。

保存实际输入、版本、命令、结果和一段解释。新计算的 `execution.json` 核验退出码与产物，候选是否达到研究目标还要继续读分数、几何和过滤标志。
