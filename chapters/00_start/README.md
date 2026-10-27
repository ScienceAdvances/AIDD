# 00 · 学习主线与练习导航

从“正在研究什么分子”开始，把结构、模型分数和实验终点连接起来。

[阅读教案](tutorial.md) · [回到总目录](../../README.md) · [结果解释](../../docs/EXPECTED_RESULTS.md)

## 输入与对应代码

全仓库共用 data/、workflows/、scripts/ 和 examples/；新结果写入 runs/。

所有命令均在仓库根目录执行。共享代码位于 `scripts/`，GNINA 的任务与辅助程序保存在 `workflows/gnina/`；本章无需复制另一份脚本。

## 开始练习

先用 `python practice.py list` 找到章节，再执行 `python practice.py doctor` 检查当前解释器。

## 独立练习与交付

选择多肽或 CAR-T 主线，写出一个研究问题、一项已有证据和一项尚需验证的假设。

保存实际输入、版本、命令、结果和一段解释。新计算的 `execution.json` 核验退出码与产物，候选是否达到研究目标还要继续读分数、几何和过滤标志。
