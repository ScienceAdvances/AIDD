# BoltzGen：Debian 五阶段设计、复核与筛选

> 本文对应网站教案的 Debian 版本。仓库练习请先读同目录 [README](README.md)：其中命令使用本仓库目录和可配置工具路径；正文的运行命令已按仓库路径调整；多工具任务统一通过可配置入口执行。


本章在 RTX 3060 上完成两个 VHH 候选的生成、逆折叠、Boltz2 复核、分析和过滤。最终排名有两行，**通过全部筛选的候选为零**。学习重点是沿着真实配置和过滤标志判断结果。

## 1. 本轮版本与阶段

复用 `/opt/micromamba/envs/boltzgen` 的 **BoltzGen 0.3.2、Python 3.12、Torch 2.6.0+cu126**。本环境的实际路线为：

```text
靶标 + VHH 框架 + CDR 定义
 → design
 → inverse_folding（boltzgen1_ifold.ckpt）
 → folding（Boltz2）
 → analysis
 → filtering
```

逆折叠在给定主链条件下生成序列；结构复核再用所得序列预测三维结构。两步的输入方向不同，不能因生成模型给出一个姿势就省略复核。

本轮没有 SolubleMPNN 阶段，也没有 ESMFold2。ESMFold2 按用户决定保留待办；其他源码提交即使显示相同版本号，也不能用来描述本环境已执行的阶段。

版本号在这里决定的不只是命令写法，也决定真正执行的科学阶段。当前逆折叠使用 `boltzgen1_ifold`，结构复核使用 Boltz2；如果直接套用其他主线版本的流程图，就会把未执行的模型写进结果。复现时应从本机 `steps.yaml` 和阶段日志重建路线，再解释各阶段作用。

## 2. 输入与设计范围

输入目录为 [workflows/modern_tools/boltzgen](../../workflows/modern_tools/boltzgen)。靶标取 `9bkq-assembly2.cif` B 链，VHH 使用 `7xl0.cif` 与官方 CDR 定义，顶层 `penguinpox.yaml` 将二者组合。

7XL0 的结构保留范围明确为 1..121，保留框架 92 个残基。三段 CDR 定义控制保留、设计、排除和插入：

CDR1 的原区域为 26..33，共 8 位；`exclude` 去掉 26..28 三位，保留 5 位，另插入 1..5 位，得到最终 6..10 位。CDR2 的原区域为 51..57，共 7 位；排除 51..53 三位后保留 4 位，再插入 1..5 位，得到 5..9 位。

CDR3 的原区域为 97..110，共 14 位；排除 97..102 六位后保留 8 位，再插入 1..12 位，得到 9..20 位。`design` 表示允许设计的区域，`exclude` 表示从原结构移除的部分，二者含义不同；可设计不等于从结构中删除。

最终 VHH 长度为 `92+CDR1+CDR2+CDR3`，范围 112..131。`structure_groups` 给框架和靶标提供内部几何条件，可设计环区则允许变化。检查这些定义能够防止把靶标或框架误设成可变位点。

来源、官方输入与派生 YAML 的哈希见[输入清单](../../workflows/modern_tools/input_manifest.json)。

设计范围可以理解为对结构中不同部分赋予不同任务：框架提供已有折叠的支撑，靶标提供待识别表面，CDR 承担主要变化。`structure_groups` 提供组内几何条件，不等于已经证明组与组之间会以生成姿势稳定结合。检查输入时，既要确认哪些坐标受到条件约束，也要确认哪些相对安排需要在后续重新预测中接受检验。

## 3. 配置、检查、执行

```bash
# 先在 config/tools.json 配置本章所用环境与模型资源
python practice.py run 15 --output runs/15_new --execute
```

`models` 对本案例执行 configure，复用已下载的资源。实际配置为 protocol=nanobody-anything、num_designs=2、budget=2、devices=1、workers=0、use_kernels=false，design 为 500 步。配置位于 `examples/debian13_20261023/boltzgen_gpu/config/`。

检查 `steps.yaml` 及各 YAML 的路径、阶段、模型和采样数，再执行已准备的任务：

```bash
# practice.py run 15 已依次执行 configure 与 execute。
# 查看 runs/15_new/boltzgen_gpu/steps.yaml 与 runs/15_new/execution.json。
```

同一目录已经完成，请先查看保存结果。做新实验应指定独立的输出目录，并记录实际配置；预算指定的是希望保留的排名规模，不承诺产生同样数量的合格设计。

`num_designs=2` 控制本轮尝试生成多少候选，`budget=2` 控制最终希望输出的排名规模。前者限定探索范围，后者限定保留数量，两者都不等于“保证获得两条合格结合分子”。若所有候选失败，报告仍可输出两行，并同时标出失败状态；因此应先读通过数，再读排名和结构。

## 4. 逐阶段验收

design 阶段应有两条初始 CIF 及元数据，inverse_folding 阶段应有两条最终序列及 CIF。folding 随后交付 Boltz2 复核结构与置信度，analysis 聚合各项指标并保持候选 ID 对应，filtering 最后交付排名、逐项过滤标志和报告。沿同一候选逐阶段追踪，能够定位几何在哪次转换发生变化。


本轮五阶段约 406 秒。四份 BoltzGen 模型资源的实际哈希见[资源指纹](../../examples/debian13_20261023/boltzgen_weight_hashes.json)。程序完成与科学条件合格分开验收。

## 5. 两候选的真实结果

排名第一的 `penguinpox_1`，design-to-target iPTM 为 0.41928，整体 filter RMSD 为 24.16721 Å，binder 自身 RMSD 为 1.33133 Å，没有通过全部筛选。排名第二的 `penguinpox_0` 对应为 0.25518、22.46054 Å 和 1.38864 Å，同样未通过全部筛选。

这两个结果共同指出，自身形状较接近设计，并没有确保复合物整体安排一致。排名只是在本轮候选中确定先后，是否达到筛选条件必须另外读取过滤标志。


两条 binder 自身折叠较接近设计，但复合物整体几何不一致。单独 binder 的较低 RMSD 不能替代结合摆放的检查。

`final_designs_metrics_2.csv` 仍有两行，因为此版本按过滤通过情况惩罚排名后，会继续输出预算数量。必须读取 `pass_filters` 和 `pass_filter_rmsd_filter`；不能只凭“最终目录中有两个结构”宣称通过筛选。

[逐候选表](../../examples/debian13_20261023/boltzgen_candidate_summary.csv)、[筛选摘要](../../examples/debian13_20261023/boltzgen_analysis_summary.json)、[原始全指标](../../examples/debian13_20261023/boltzgen_gpu/final_ranked_designs/all_designs_metrics.csv)。复算入口：

```bash
python scripts/analyze_debian_boltzgen.py --output runs/15_analysis
```

设想一个 VHH 保持自己的折叠形状，却相对抗原整体平移或转动。先把 VHH 自己叠合，就可以得到很小的 binder RMSD；把两者作为复合物检查，则会保留相对摆放变化，出现很大的整体误差。因此，两种 RMSD 并不矛盾。进一步定位时，应固定靶标再查看 VHH 的位置、朝向及原本设计的接触，而不能只看单体叠合。

`pass_filters` 回答是否满足整套条件，`pass_filter_rmsd_filter` 回答其中的几何条件是否通过。最终排名回答的是当前候选怎样排序，这三个问题需要分别读取。当前零通过应如实保留，让后续优化针对实际失败项展开。

## 6. 怎样解释负结果并设计下一轮

本轮只生成两个候选，用来检验流程和指标；零通过不能推导靶标不可设计。下一轮先叠合靶标，查看 binder 的相对位置变化，再检查表位、框架可及性和 CDR 几何。扩大采样时保持参数一致，再判断同类失败是否普遍。

如果只改变过滤阈值，应说明这是筛选敏感性分析，不能将放宽阈值视为生成质量改善。实验结合、特异性、表达和聚集仍需要独立验证。

若下一轮怀疑问题来自表位或框架朝向，可以固定采样规模，分别改变这些条件并记录复核后的相对摆放；若怀疑只是小样本没有探索到可行解，则保持输入与筛选规则，增加独立候选并比较失败模式。一次改变多个因素，虽然可能得到更好结果，却难以解释改善来自哪里。

放宽 RMSD 阈值适合检查结果对筛选规则是否敏感，但不会改变已经生成的结构。若新增宽松筛选，应保留原过滤结果和新规则，明确它只是候选晋级标准的变化。后续实验仍需检验表达、聚集、选择性和真实结合。

## 7. ESMFold2 待办

本机驱动与官方 CUDA 13 运行时要求不兼容，官方小例显存测量也超过 12GB。用户选择先交付现有结果，ESMFold2 保留待办，详见[兼容性评估](../../examples/debian13_20261023/ESMFold2_兼容性评估.md)。后续创建环境、驱动改动或使用其他 GPU 均需按既定边界确认。

## 8. 练习交付

从最终表逐条说明通过和失败的条件；把 binder 自身 RMSD 与靶标对齐后的相对位置画在同一份检查记录中；解释为何排名第一也不合格。交付内容应包括配置、两条序列、结构、完整过滤标志和下一轮可检验假设。
