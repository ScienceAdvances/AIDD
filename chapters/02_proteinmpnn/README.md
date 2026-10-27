# 02 · ProteinMPNN 序列设计

在已知 1YCR 骨架上保持受体和三个肽热点，核对六条设计的计分范围。

[阅读教案](tutorial.md) · [回到总目录](../../README.md) · [结果解释](../../docs/EXPECTED_RESULTS.md)

## 输入与对应代码

data/peptide/1YCR_original.pdb；examples/proteinmpnn/；data/aidd/ProteinMPNN/。

所有命令均在仓库根目录执行。共享代码位于 `scripts/`，GNINA 的任务与辅助程序保存在 `workflows/gnina/`；本章无需复制另一份脚本。

## 先复算已保存结果

```bash
python practice.py analyze 02 --output runs/02_analysis
```

此入口读取仓库中的 Debian 示例并重新计算分析，不会启动 GPU 设计模型。

验收时应看到：六条序列，肽链 13 位、固定 3/7/10、十个可变计分位点；重算 NLL 应与 FASTA 头部舍入值相符。

## 再执行新的计算

先按[环境说明](../../docs/ENVIRONMENTS.md)配置工具。下面第一条只显示命令计划，确认路径与资源后再执行第二条；输出目录必须尚不存在。

```bash
python practice.py run 02 --output runs/02_new
python practice.py run 02 --output runs/02_new --execute
```

重新推理前先执行 `python practice.py fetch proteinmpnn_weights`。第 06 章默认对随仓库提供的六条 Debian 设计评分；若要换成自己的六条设计，直接调用 `scripts/optimize_sequence.py --design-run runs/02_new/design --output runs/06_own --device cuda`。

## 独立练习与交付

把 PDB B19/B23/B26 映射成链内 3/7/10，从 NPZ 的 mask 和 chain_order 验证受体与热点保持。

保存实际输入、版本、命令、结果和一段解释。新计算的 `execution.json` 核验退出码与产物，候选是否达到研究目标还要继续读分数、几何和过滤标志。
