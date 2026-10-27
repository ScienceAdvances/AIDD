# 13 · AlphaFold 3 复合物预测

从序列、空 MSA 和空模板开始，按完整精度排名选择样本，再与参考比较。

[阅读教案](tutorial.md) · [回到总目录](../../README.md) · [结果解释](../../docs/EXPECTED_RESULTS.md)

## 输入与对应代码

workflows/modern_tools/alphafold3/；examples/debian13_20261023/alphafold3_gpu_repaired/。

所有命令均在仓库根目录执行。共享代码位于 `scripts/`，GNINA 的任务与辅助程序保存在 `workflows/gnina/`；本章无需复制另一份脚本。

## 先复算已保存结果

```bash
python practice.py analyze 13 --output runs/13_analysis
```

此入口读取仓库中的 Debian 示例并重新计算分析，不会启动 GPU 设计模型。

验收时应看到：5 个样本；最高排名为 sample 3，iPTM 0.17、受体对齐后 p53 CA RMSD 20.801 Å。本轮未可靠恢复参考。

## 再执行新的计算

先按[环境说明](../../docs/ENVIRONMENTS.md)配置工具。下面第一条只显示命令计划，确认路径与资源后再执行第二条；输出目录必须尚不存在。

```bash
python practice.py run 13 --output runs/13_new
python practice.py run 13 --output runs/13_new --execute
```

AlphaFold 3 模型参数需自行按官方条件获取；本仓库不分发模型参数。空 MSA 输入跳过数据库搜索，仍需可工作的官方软件、CCD 数据和参数。

## 独立练习与交付

核对链长度与 PAE 区块，复算五份结构。说明这个负结果为何不能推翻已知结合，也不能独自证明 MSA 的因果作用。

保存实际输入、版本、命令、结果和一段解释。新计算的 `execution.json` 核验退出码与产物，候选是否达到研究目标还要继续读分数、几何和过滤标志。
