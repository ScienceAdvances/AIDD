# PyRosetta：给 p53 候选换侧链、精修并分析界面

> 本文对应网站教案的 Debian 版本。仓库练习请先读同目录 [README](README.md)：其中命令使用本仓库目录和可配置工具路径；正文的运行命令已按仓库路径调整；多工具任务统一通过可配置入口执行。


第六章选出了 WT、D5E、S4E 和 L9K。ProteinMPNN 的分数反映这些序列与给定结构的相容性；现在要进一步看：把突变真实地放回三维结构，侧链能否合理排列，局部精修后得到什么能量和界面变化。

这一课使用**人源 MDM2—p53 的 1YCR**，与第二、六章共享同一序列体系。最终得到四份精修 PDB 和一张统一协议下的评分表。Rosetta 的 Python 接口名为 **PyRosetta**，下面直接从 Python 调用 Rosetta 的打分、侧链装填和最小化工具。

## 1. Debian 实际版本

本轮复用现有 **PyRosetta 2026.3+releasequarterly.5e498f1409、Python 3.11**，使用 CPU。直接检查已安装环境：

```bash
cd .
python -c 'import pyrosetta; print(pyrosetta.version())'
```

打分函数为 ref2015，固定主链和刚体 jump，允许所有侧链重排与最小化。具体协议和版本串见[运行摘要](../../examples/debian13_20261023/pyrosetta/run_summary.json)。许可与 API 见[官方说明](https://www.pyrosetta.org/downloads)。

## 2. 先确认输入和编号

原始结构为 [1YCR_original.pdb](../../data/peptide/1YCR_original.pdb)。本课只取蛋白 ATOM 记录，得到 A 链 85 个残基、B 链 13 个残基，总共 98 个。B 链为 `ETFSDLWKLLPEN`，PDB 编号从 17 到 29。

PyRosetta 用 **Pose** 保存结构、序列及能量信息。Pose 内部把所有残基按顺序编号；PDB 的链号和残基号由 `PDBInfo` 保存。两套编号之间要显式转换：

D5E 在肽内是第 5 位，对应 PDB B21 和 Pose 90；S4E 是第 4 位，对应 B20 和 Pose 89；L9K 是第 9 位，对应 B25 和 Pose 94。三个固定热点 F/W/L 的肽内位置是 3/7/10，对应 PDB B19/B23/B26 和 Pose 88/92/95。先逐项核对原残基身份，再执行突变，能够避免把一条语法正确的命令用到错误位置。


```python
import pyrosetta
pyrosetta.init('-mute all -constant_seed -jran 20261020 -ex1 -ex2aro')
pose = pyrosetta.pose_from_pdb(
    'examples/debian13_20261023/pyrosetta/1YCR_protein.pdb'
)
assert pose.total_residue() == 98
assert pose.chain_sequence(2) == 'ETFSDLWKLLPEN'
position = pose.pdb_info().pdb2pose('B', 21)
assert position == 90
print(pose.residue(position).name3())
```

这段检查使用教材已经保存的清理后 PDB；完整脚本会从原始文件重新制作它。`pdb2pose` 返回 0 表示没有找到，必须先排查链号、插入码或编号，再执行突变。

## 3. 这次精修具体允许什么移动

本课固定主链和两条链的相对位置，允许全部蛋白残基的侧链重新装填，再最小化侧链角度。这样四个候选共享同一个结合主链，可把注意力放在侧链适配上。

**Rotamer，侧链旋转异构体**，是常见侧链二面角组合。装填程序从这些组合中搜索相互相容的排列。**Repacking** 保持氨基酸身份，仅重新选择侧链；**Design** 则还允许改变残基种类。本课先指定单点突变，再使用 repacking。

**MoveMap** 指定最小化可调整的自由度。`bb` 对应主链角度，`chi` 对应侧链角度，`jump` 管理刚体连接。本例有一个 jump 连接两条链，固定它就保持当前复合物摆放。

```python
from pyrosetta import rosetta

scorefxn = pyrosetta.create_score_function('ref2015')
rosetta.protocols.simple_moves.MutateResidue(position, 'GLU').apply(pose)
task = pyrosetta.standard_packer_task(pose)
task.restrict_to_repacking()
task.or_include_current(True)
rosetta.protocols.minimization_packing.PackRotamersMover(
    scorefxn, task
).apply(pose)

move_map = rosetta.core.kinematics.MoveMap()
move_map.set_bb(False)
move_map.set_chi(True)
move_map.set_jump(False)
minimizer = rosetta.protocols.minimization_packing.MinMover()
minimizer.movemap(move_map)
minimizer.score_function(scorefxn)
minimizer.min_type('lbfgs_armijo_nonmonotone')
minimizer.tolerance(0.001)
minimizer.apply(pose)
```

`ref2015` 将范德华、溶剂化、氢键及构象偏好等项按权重组合，返回 Rosetta Energy Units，简称 **REU**。最小化在给定自由度内沿能量下降方向调整几何。其优化范围比整体重新折叠小，适合本课的局部比较。原理可参照 [Rosetta 打分函数文档](https://docs.rosettacommons.org/docs/latest/rosetta_basics/scoring/score-types)。

固定主链和 jump 同时带来可解释性与限制。好处是四条序列都在同一个结构假设上比较，差别主要来自侧链适配；限制是肽不能整体离开口袋，也不能通过显著主链重排寻找其他结合方式。因此，精修后仍贴在受体上并不构成结合稳定性的证据——这种贴合本来就是本次协议固定的条件。

侧链最小化与分子动力学也不同。它沿允许自由度降低当前评分，不以真实时间描述温度下的构象分布。若需要研究热运动，应建立相应力场体系；若需要比较更大构象调整，可新增受约束主链精修作为另一套实验，并保留它与本协议的区别。

## 4. 分析结合态与分离态

**InterfaceAnalyzer** 在既定界面划分下比较复合物与分离组分，并计算界面埋藏面积等指标。这里两条链恰好由 jump 1 分开：

```python
assert pose.num_jump() == 1
interface = rosetta.protocols.analysis.InterfaceAnalyzerMover(1)
interface.set_scorefunction(scorefxn)
interface.set_pack_input(False)
interface.set_pack_separated(True)
interface.apply(pose)
print(interface.get_interface_dG())
print(interface.get_interface_delta_sasa())
```

本版 Python 绑定采用 `InterfaceAnalyzerMover(1)` 这条实测入口。换成多链体系时，先检查 FoldTree 和 jump 划分；应让分组对应研究中的受体与结合分子。

`pack_input=False` 保留刚完成的结合态精修结果，`pack_separated=True` 允许分离后的组分重新装填侧链。由于分离态也会调整，这些选项会影响界面分数，比较候选时应保持一致。**ΔSASA** 是组分分离后增加的溶剂可及面积，单位 Å²，可帮助描述界面接触规模。[InterfaceAnalyzer 官方说明](https://docs.rosettacommons.org/docs/latest/application_documentation/analysis/interface-analyzer)解释了对应的能量与面积指标。

可以把界面分数差概括为 $E_{complex}-E_{separated}$。负值表示按当前函数和分离态准备方式，复合物比指定分离状态更受偏好。由于本例允许分离组分重新装填侧链，分离态不是简单把原坐标平移开后直接读数；这也是 `pack_separated` 必须随协议保存的原因。

ΔSASA 则是面积变化：两个组分接触后，有一部分原本能接触溶剂的表面被埋藏。面积较大意味着接触规模更大，却不保证接触化学合适；埋藏未满足的极性基团或不合理的疏水表面，都需要进一步检查。因此，面积、界面分数和具体氢键或碰撞应该一起读。

## 5. 运行整个候选组

完整脚本 [run_pyrosetta_case.py](../../scripts/run_pyrosetta_case.py)逐个从同一原始 Pose 克隆 WT 与三个变体，每个候选重置同一随机种子，然后执行上面的统一协议。

```bash
python scripts/run_pyrosetta_case.py \
  --output runs/10_my_pyrosetta
```

完成后检查 `run_summary.json` 的版本、原始文件哈希、`backbone_fixed` 和 `pack_separated`；再核对 `scores.csv` 的四条序列。`WT.pdb`、`D5E.pdb`、`S4E.pdb` 和 `L9K.pdb` 可以叠合查看。把主链固定以后，CA 坐标应保留输入位置，变化集中在侧链。

## 6. 实测结果，以及为什么会出现不同排序

本地记录见 [运行摘要](../../examples/debian13_20261023/pyrosetta/run_summary.json)和[完整评分表](../../examples/debian13_20261023/pyrosetta/scores.csv)。

WT 的总分从精修前 161.746 降至 −175.327 REU，界面分数差为 −26.914 REU，ΔSASA 为 1491.39 Å²。D5E 从 170.086 降至 −178.134 REU，界面分数差为 −26.719 REU，ΔSASA 为 1497.10 Å²。D5E 的精修后总分更低，但其界面列并没有同样领先。

S4E 从 439.011 降至 −175.903 REU，界面分数差为 −26.366 REU，ΔSASA 为 1539.51 Å²。L9K 从 186.518 降至 −172.323 REU，界面分数差为 −27.100 REU，ΔSASA 为 1492.46 Å²。L9K 的界面列最有利，S4E 的埋藏面积最大，这些观察各自描述不同方面，不能合并成一个没有定义的“最好”。


S4E 突变后、装填前的总分很高，经过侧链搜索明显下降。这提醒我们：直接把更大的侧链塞进原位置，会先产生不合适的接触；候选比较应使用统一精修后的状态。

D5E 的精修后总分最低，而界面分数最有利的是 L9K。两列优化的对象不同：总分覆盖整个复合物，界面分数比较结合态和指定分离态。把第六章的 ProteinMPNN 分数一并放进来，可以形成三项互补观察：序列适配、局部物理评分和界面接触。

这组结果来自单一主链、单一种子和一套侧链协议。当前可据此提出结构假设：优先检查 D5E 在局部结构中的适配，以及 L9K 的界面相互作用；候选提名继续结合溶解性、实验成本和已保留的热点。要建立稳定排序，下一轮应增加结构起点与随机种子，检查排序是否保持。

以 S4E 为例，较大的 Glu 侧链直接放到 Ser 的原位置时可能产生局部冲突，使精修前分数很高；装填后下降到与其他候选接近的范围，说明先前一部分高分来自未适配的几何。不能把这个大幅下降直接归因于突变提高了亲和力，因为 WT 和其他变体也需要经过同样准备。

D5E 与 WT 的界面分数只相差约 0.195 REU，L9K 相对 WT 也只改善约 0.186 REU。当前只有一个起点和一种随机种子，尚未量化协议波动；这些小差值更适合用于提出检查局部接触的假设，而非宣称已经建立可靠的活性顺序。

## 7. 从局部精修走向可靠比较

改变氨基酸会同时改变参考能量、侧链自由度、溶剂化和内部分子接触。跨序列比较总分时，应明确这些项都参与了变化。界面分析的分离态选择也会改变结果。面试里说明协议与比较范围，比只报一个“最低分”更有用。

如果研究问题需要主链重排，可另建一个受约束的 Relax 协议，并同时评估与参考结构的偏移。新协议应给全部候选使用相同约束、采样规模和评价标准；保留本课的固定主链结果，便于判断新增自由度带来了什么变化。

找不到安装包时核对官方 wheel 索引、Python 和系统架构；突变位置异常时，同时打印 PDB 编号、Pose 编号和原残基。自有配体或其他非标准残基需要对应 Rosetta 参数，应先单独验证读取和打分。

换成其他体系后若界面为空，先检查链、FoldTree 和 jump 分组。两次评分不一致则核对随机种子、打分权重、装填范围和最小化自由度。只有比较协议一致，数值差才具有明确含义。


## 8. 面试短答与练习

本课完成的是受约束的局部比较：把三个单点候选映射回 1YCR，让 WT 与变体共享主链、链间位置、装填和最小化协议，再分别读取总分与界面指标。D5E 总分最低，L9K 界面列更有利，下一步应检查具体接触及其对起点、随机种子的敏感性。

解释结果时，REU 保留为当前评分函数的单位。若希望评价它能否区分实验亲和力，需要在条件一致的结合测定数据上检验关系。精修前后分数与结构位移也应同时保存，使读者能够判断改善来自修复局部冲突，还是依赖更大的结构重排。

练习交付：把四份结果叠合，画出突变侧链与邻近残基；随后以三个新随机种子复算，整理每个候选的分数范围。能解释“为何允许这些运动、为何用这项分数比较”，就掌握了本课核心。
