# 11 · RFdiffusion2 原子基序支架

围绕十个功能原子生成支架，再检查序列设计与 Chai 复核后的配体几何。

[阅读教案](tutorial.md) · [回到总目录](../../README.md) · [结果解释](../../docs/EXPECTED_RESULTS.md)

## 输入与对应代码

workflows/modern_tools/rfdiffusion2/；examples/debian13_20261023/rfdiffusion2_outputs/compiled_metrics.csv。

所有命令均在仓库根目录执行。共享代码位于 `scripts/`，GNINA 的任务与辅助程序保存在 `workflows/gnina/`；本章无需复制另一份脚本。

## 先复算已保存结果

```bash
python practice.py analyze 11 --output runs/11_analysis
```

此入口读取仓库中的 Debian 示例并重新计算分析，不会启动 GPU 设计模型。

验收时应看到：1 个主链、8 条 LigandMPNN 序列、40 个 Chai 样本；每序列 best-of-five 的口袋对齐配体 RMSD 仍为 6.49–11.40 Å。

## 再执行新的计算

先按[环境说明](../../docs/ENVIRONMENTS.md)配置工具。下面第一条只显示命令计划，确认路径与资源后再执行第二条；输出目录必须尚不存在。

```bash
python practice.py run 11 --output runs/11_new
python practice.py run 11 --output runs/11_new --execute
```

原生适配器替换容器启动方式，保留官方科学阶段。需要先配置 RFdiffusion2、LigandMPNN、Chai 的解释器、源码和已获许可的模型资源；若只学习指标，可停留在上面的离线复算。完整 40 个复核结构见可选数据包。

## 独立练习与交付

从实际列名提取指标，区分生成约束几乎为零的基序 RMSD 与复核几何。说明 best-of-five 的选择偏好和关闭可选 ESM 嵌入的设置。

保存实际输入、版本、命令、结果和一段解释。新计算的 `execution.json` 核验退出码与产物，候选是否达到研究目标还要继续读分数、几何和过滤标志。
