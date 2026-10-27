# 04 · CAR-T 原始数据与工程判断

从公开工作簿逐单元格提取观测，比较表达、背景、抗原结合和功能之间的取舍。

[阅读教案](tutorial.md) · [回到总目录](../../README.md) · [结果解释](../../docs/EXPECTED_RESULTS.md)

## 输入与对应代码

data/cart/ 两份 Source Data 工作簿；examples/debian13_20261023/cart/。

所有命令均在仓库根目录执行。共享代码位于 `scripts/`，GNINA 的任务与辅助程序保存在 `workflows/gnina/`；本章无需复制另一份脚本。

## 先复算已保存结果

```bash
python practice.py analyze 04 --output runs/04_analysis
```

此入口会从原始工作簿重新提取和绘图。

验收时应看到：提取 204 个观测、49 个分组；L3 相对背景 30.372、相对存活 9.289，H2L3wt 分别为 12.657、80.131。

## 再执行新的计算

先按[环境说明](../../docs/ENVIRONMENTS.md)配置工具。下面第一条只显示命令计划，确认路径与资源后再执行第二条；输出目录必须尚不存在。

```bash
python practice.py run 04 --output runs/04_new
python practice.py run 04 --output runs/04_new --execute
```

## 独立练习与交付

从一条观测追溯工作簿、sheet 与单元格，重算均值、SD 和 SEM，再解释背景最低为何不一定最值得提名。

保存实际输入、版本、命令、结果和一段解释。新计算的 `execution.json` 核验退出码与产物，候选是否达到研究目标还要继续读分数、几何和过滤标志。
