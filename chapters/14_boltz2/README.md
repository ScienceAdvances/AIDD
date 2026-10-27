# 14 · Boltz-2 结构与亲和力输出

把蛋白—肽结构预测与独立蛋白—小分子亲和力示例分开评价。

[阅读教案](tutorial.md) · [回到总目录](../../README.md) · [结果解释](../../docs/EXPECTED_RESULTS.md)

## 输入与对应代码

workflows/modern_tools/boltz/；examples/debian13_20261023/boltz_structure/ 与 boltz_affinity/。

所有命令均在仓库根目录执行。共享代码位于 `scripts/`，GNINA 的任务与辅助程序保存在 `workflows/gnina/`；本章无需复制另一份脚本。

## 先复算已保存结果

```bash
python practice.py analyze 14 --output runs/14_analysis
```

此入口读取仓库中的 Debian 示例并重新计算分析，不会启动 GPU 设计模型。

验收时应看到：结构最高置信度 model_0 的 iPTM≈0.956、p53 RMSD≈2.052 Å；独立小分子 y=2.2371197，预测 IC50≈172.63 µM，分类概率≈0.3534。

## 再执行新的计算

先按[环境说明](../../docs/ENVIRONMENTS.md)配置工具。下面第一条只显示命令计划，确认路径与资源后再执行第二条；输出目录必须尚不存在。

```bash
python practice.py run 14 --output runs/14_new
python practice.py run 14 --output runs/14_new --execute
```

第一个任务使用空 MSA；第二个官方亲和力任务启用 MSA server，会向配置的服务提交蛋白序列，并按需下载模型。

## 独立练习与交付

复算受体对齐后的肽偏差，并从原始 affinity JSON 换算单位。说明 y=−3 对应 1 nM，但模型 IC50 不能直接改名为 Kd。

保存实际输入、版本、命令、结果和一段解释。新计算的 `execution.json` 核验退出码与产物，候选是否达到研究目标还要继续读分数、几何和过滤标志。
