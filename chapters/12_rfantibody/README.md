# 12 · RFantibody 抗体与 VHH

按框架与抗原两种基准对齐，检查 VHH 局部形状与整体结合摆放。

[阅读教案](tutorial.md) · [回到总目录](../../README.md) · [结果解释](../../docs/EXPECTED_RESULTS.md)

## 输入与对应代码

workflows/modern_tools/rfantibody/；examples/debian13_20261023/rfantibody_outputs/ 的三阶段 QV、评分及最终 PDB。

所有命令均在仓库根目录执行。共享代码位于 `scripts/`，GNINA 的任务与辅助程序保存在 `workflows/gnina/`；本章无需复制另一份脚本。

## 先复算已保存结果

```bash
python practice.py analyze 12 --output runs/12_analysis
```

此入口读取仓库中的 Debian 示例并重新计算分析，不会启动 GPU 设计模型。

验收时应看到：2 个主链 → 4 条序列 → 4 个 RF2 结构；两条相对抗原偏移超过 42 Å，另有一条约 1.15 Å、界面 PAE 4.68 Å。

## 再执行新的计算

先按[环境说明](../../docs/ENVIRONMENTS.md)配置工具。下面第一条只显示命令计划，确认路径与资源后再执行第二条；输出目录必须尚不存在。

```bash
python practice.py run 12 --output runs/12_new
python practice.py run 12 --output runs/12_new --execute
```

使用官方脚本和可写工作目录完成后续阶段，并核对每阶段实际文件。RF2 的 num_recycles=10 表示初次预测再加十次 recycle。

## 独立练习与交付

复算 QV 数量、目标和框架序列保持，以及靶标对齐后的抗体 RMSD。用叠合图解释高整体 pLDDT 为什么不能替代界面检查。

保存实际输入、版本、命令、结果和一段解释。新计算的 `execution.json` 核验退出码与产物，候选是否达到研究目标还要继续读分数、几何和过滤标志。
