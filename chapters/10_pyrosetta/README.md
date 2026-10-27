# 10 · PyRosetta 侧链精修与界面

固定主链与链间位置，比较四个候选的侧链适配、总分和界面指标。

[阅读教案](tutorial.md) · [回到总目录](../../README.md) · [结果解释](../../docs/EXPECTED_RESULTS.md)

## 输入与对应代码

data/peptide/1YCR_original.pdb；examples/debian13_20261023/pyrosetta/。

所有命令均在仓库根目录执行。共享代码位于 `scripts/`，GNINA 的任务与辅助程序保存在 `workflows/gnina/`；本章无需复制另一份脚本。

## 先复算已保存结果

```bash
python practice.py analyze 10 --output runs/10_analysis
```

此入口读取仓库中的 Debian 示例并重新计算分析，不会启动 GPU 设计模型。

验收时应看到：四份结构主链 CA 坐标保持；D5E 精修后总分约 −178.134 REU，L9K 界面列约 −27.100 REU。

## 再执行新的计算

先按[环境说明](../../docs/ENVIRONMENTS.md)配置工具。下面第一条只显示命令计划，确认路径与资源后再执行第二条；输出目录必须尚不存在。

```bash
python practice.py run 10 --output runs/10_new
python practice.py run 10 --output runs/10_new --execute
```

## 独立练习与交付

将 WT 与三个变体叠合，查看突变侧链和邻近原子。解释为什么总分最低、界面分数最低和面积最大不一定对应同一候选。

保存实际输入、版本、命令、结果和一段解释。新计算的 `execution.json` 核验退出码与产物，候选是否达到研究目标还要继续读分数、几何和过滤标志。
