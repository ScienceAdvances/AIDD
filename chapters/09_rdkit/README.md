# 09 · RDKit 化学身份与构象

把分子连接关系、化学状态和三维构象分别记录，为对接准备可追踪输入。

[阅读教案](tutorial.md) · [回到总目录](../../README.md) · [结果解释](../../docs/EXPECTED_RESULTS.md)

## 输入与对应代码

workflows/gnina/example_1iep/reference.sdf；scripts/run_rdkit_case.py。

所有命令均在仓库根目录执行。共享代码位于 `scripts/`，GNINA 的任务与辅助程序保存在 `workflows/gnina/`；本章无需复制另一份脚本。

## 先复算已保存结果

```bash
python practice.py analyze 09 --output runs/09_analysis
```

此入口会重新生成 RDKit 构象并计算性质。

验收时应看到：37 个重原子、形式电荷 +1；请求 20 个，保存 18 个且 MMFF 均收敛。不同 RDKit 版本可有数值或构象数量变化。

## 再执行新的计算

先按[环境说明](../../docs/ENVIRONMENTS.md)配置工具。下面第一条只显示命令计划，确认路径与资源后再执行第二条；输出目录必须尚不存在。

```bash
python practice.py run 09 --output runs/09_new
python practice.py run 09 --output runs/09_new --execute
```

增加构象数的独立练习：`python scripts/run_rdkit_case.py --num-confs 50 --output runs/09_50confs`。

## 独立练习与交付

先检查化学身份，再将 --num-confs 改成 50，使用新目录运行，比较构象覆盖和能量分布。解释 MMFF 能量为何不是结合自由能。

保存实际输入、版本、命令、结果和一段解释。新计算的 `execution.json` 核验退出码与产物，候选是否达到研究目标还要继续读分数、几何和过滤标志。
