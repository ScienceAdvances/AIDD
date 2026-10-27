# 05 · 跨案例证据与综合简报

把不同体系的证据整理为简报，并为每个候选提出可检验的下一步。

[阅读教案](tutorial.md) · [回到总目录](../../README.md) · [结果解释](../../docs/EXPECTED_RESULTS.md)

## 输入与对应代码

examples/debian13_20261023/ 下的 HADDOCK、序列、MD、CAR-T 和 GNINA 记录。

所有命令均在仓库根目录执行。共享代码位于 `scripts/`，GNINA 的任务与辅助程序保存在 `workflows/gnina/`；本章无需复制另一份脚本。

## 先复算已保存结果

```bash
python practice.py analyze 05 --output runs/05_analysis
```

此入口读取仓库中的 Debian 示例并重新计算分析，不会启动 GPU 设计模型。

验收时应看到：生成 evidence.csv、candidate_validation_plan.csv、brief.md 与 build_summary.json；提名仍为 WT 与三个单点。

## 独立练习与交付

任选 L9K 写一句包含比较对象、NLL 差值、评价协议和下一步的结论。说明 Factor Xa MD 为何不能验证这条肽。

保存实际输入、版本、命令、结果和一段解释。新计算的 `execution.json` 核验退出码与产物，候选是否达到研究目标还要继续读分数、几何和过滤标志。
