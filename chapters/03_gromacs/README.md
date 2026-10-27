# 03 · GROMACS 分子动力学

先从一条已完成的 1 ns 轨迹学会读图，再准备自己的水、离子、拓扑和平衡过程。

[阅读教案](tutorial.md) · [回到总目录](../../README.md) · [结果解释](../../docs/EXPECTED_RESULTS.md)

## 输入与对应代码

data/md/md-intro-tutorial-main/data/input/；examples/debian13_20261023/md_analysis/。

所有命令均在仓库根目录执行。共享代码位于 `scripts/`，GNINA 的任务与辅助程序保存在 `workflows/gnina/`；本章无需复制另一份脚本。

## 先复算已保存结果

```bash
python practice.py analyze 03 --output runs/03_analysis
```

此入口读取仓库中的 Debian 示例并重新计算分析，不会启动 GPU 设计模型。

验收时应看到：本轮 201 帧、5 ps 间隔，Cα 对 EM 平均 RMSD 0.103076 nm，平均 Rg 1.889838 nm。

## 再执行新的计算

先按[环境说明](../../docs/ENVIRONMENTS.md)配置工具。下面第一条只显示命令计划，确认路径与资源后再执行第二条；输出目录必须尚不存在。

```bash
python practice.py run 03 --output runs/03_new
python practice.py run 03 --output runs/03_new --execute
```

入口运行 NVT/NPT 各 100 ps、生产 1 ns，并要求 CUDA 非键相互作用卸载。完整原始轨迹分析需要 GROMACS、MDTraj 和官方参考：先 `python practice.py fetch gromacs_reference`，再运行 `python scripts/analyze_md.py --gmx gmx --local-run runs/03_new/md --output runs/03_full_analysis`。

## 独立练习与交付

按 CSV 重画 RMSD，并说明周期边界、拟合对象和参考如何影响解释。包含的 protein_centered.xtc 只有蛋白坐标；完整含水轨迹见可选数据包。

保存实际输入、版本、命令、结果和一段解释。新计算的 `execution.json` 核验退出码与产物，候选是否达到研究目标还要继续读分数、几何和过滤标志。

## 可选：直接从蛋白轨迹复算

具备 MDTraj 后，运行 `python scripts/verify_md_subset.py --output runs/03_xtc_check`，重新读取 201 帧蛋白轨迹，计算 Cα 对完整 EM 参考的 RMSD，并与 CSV 比较。此练习不依赖 GROMACS，也不包含水或离子分析。
