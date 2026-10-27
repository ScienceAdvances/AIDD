# RDKit：把分子身份、三维构象和候选表准备好

> 本文对应网站教案的 Debian 版本。仓库练习请先读同目录 [README](README.md)：其中命令使用本仓库目录和可配置工具路径；正文的运行命令已按仓库路径调整；多工具任务统一通过可配置入口执行。


在 GNINA 课里，软件接收的是配体文件。真正动手时，问题往往更早出现：文件能打开，却少了电荷；同一个分子用了两种 SMILES；生成了三维坐标，却把原来的晶体姿势覆盖了。这一课用第七章的 **1IEP—伊马替尼**，完整做一遍读取、身份检查、三维构象生成和结果保存，再处理第六章的四条 p53 肽。

学完后应能独立交付一份带分子身份的 SDF、一张构象能量表和一张肽性质表。面试重点是说明每一步解决什么输入问题，以及性质、构象能量和结合能力之间的关系。

## 1. 版本与现有分析环境

本轮 Debian 使用 **RDKit 2026.09.1**（安装包版本 2026.9.1），复用已有分析环境：

```bash
cd .
python -c 'import rdkit; print(rdkit.__version__)'
```

RDKit 的身份检查和构象生成由 CPU 完成。输入、随机种子和依赖版本随输出保存，参考 API 见[官方文档](https://www.rdkit.org/docs/GettingStartedInPython.html)。

## 2. 输入是什么，为什么选这份文件

输入为 [reference.sdf](../../workflows/gnina/example_1iep/reference.sdf)，来源见第七章的[输入准备](../07_gnina/tutorial.md)。它保存伊马替尼的原子、键、电荷和晶体坐标，受体另存在 `receptor.pdb`。本课把参考 SDF 当作分子身份和已知姿势的对照，生成的构象写到新文件。

**SMILES** 把分子的连接关系写成字符串，手性和电荷也可写入其中。**SDF** 则可以同时保存这些化学信息、坐标和自定义属性。只拿 PDB 中的一组配体坐标做化学处理时，还需要补足可靠的键级、芳香性和电荷信息；本课已有 SDF，因此直接从它开始。

先亲自读取并检查：

```python
from rdkit import Chem

supplier = Chem.SDMolSupplier(
    'workflows/gnina/example_1iep/reference.sdf', removeHs=False
)
reference = next(iter(supplier))
if reference is None:
    raise ValueError('Cannot parse the reference SDF')
molecule = Chem.RemoveHs(reference)
print(Chem.MolToSmiles(molecule, isomericSmiles=True))
print(molecule.GetNumHeavyAtoms(), Chem.GetFormalCharge(molecule))
```

这段代码在教材根目录运行。`Mol` 是 RDKit 的分子对象；读取失败会得到 `None`，因此先检查再调用属性。`RemoveHs` 暂时移除显式氢，便于按重原子检查身份；氮上的形式电荷会保留。

本地读到 **37 个重原子、总形式电荷 +1**，规范异构 SMILES 为：

```text
Cc1ccc(NC(=O)c2ccc(CN3CC[NH+](C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1
```

其中 `[NH+]` 直接写出了质子化氮。后续对接若换成另一种质子化状态，受体中的静电作用和氢键关系也会变化，应在任务表中单独记录这一状态。

化学身份和三维姿势需要分别保存。两份 SDF 可以有相同原子、键与手性，却采用完全不同的坐标；同一组看起来相似的坐标，也可能因键级或质子化不同而代表不同化学状态。因此，本课保留原参考对象，再为新构象建立独立对象。这样进入对接后，能够明确分辨“起始形状”与“用来评价的实验姿势”。

## 3. 从连接关系生成一组新的三维构象

分子有许多单键可旋转，同一连接关系可对应不同的空间形状。**构象采样**就是生成这些形状。本课用 ETKDGv3：它把距离几何采样与已知的扭转角知识结合，先得到一组合理的起始坐标。[官方三维构象说明](https://www.rdkit.org/docs/GettingStartedInPython.html#working-with-3d-molecules)提供基础用法。

```python
from rdkit.Chem import AllChem

generated = Chem.AddHs(molecule)
generated.RemoveAllConformers()
parameters = AllChem.ETKDGv3()
parameters.randomSeed = 20261020
parameters.numThreads = 1
parameters.pruneRmsThresh = 0.5
conformer_ids = list(AllChem.EmbedMultipleConfs(
    generated, numConfs=20, params=parameters
))
if not conformer_ids:
    raise RuntimeError('No conformer was generated')
if not AllChem.MMFFHasAllMoleculeParams(generated):
    raise RuntimeError('MMFF parameters are incomplete')
energies = AllChem.MMFFOptimizeMoleculeConfs(
    generated, numThreads=1, maxIters=2000
)
```

先补显式氢，再清除构象，意味着新坐标从化学结构重新生成；原始 `reference` 对象仍保存晶体姿势。随机种子固定采样起点，单线程减少并行调度差异。`pruneRmsThresh=0.5` 让嵌入程序剔除过于相近的构象，所以请求 20 个并不保证保存 20 个。

**MMFF** 是分子力场，用键长、键角、扭转角和非键相互作用等能量项优化分子内部几何。这里每个返回项是 `(状态码, 能量)`，状态码 `0` 表示达到优化收敛条件。检查参数覆盖可以提前发现不支持的化学基团。

计算发生在独立配体上，因此这张能量表主要帮助排查高应变构象。进入蛋白口袋以后，配体会同时受到蛋白和溶剂环境的影响，候选姿势还要经过 GNINA 等后续流程评价。

最低 MMFF 能量的构象，只是在本次独立配体优化与参数条件下更受偏好。进入口袋后，原先不占优势的构象可能通过更好的蛋白接触得到补偿，结合时还涉及溶剂与构象熵。因此，准备构象库时需要在几何合理的前提下保留一定多样性，让后续对接有机会探索不同起点。

“18 个都收敛”则描述数值优化达到了停止条件，并不意味着找到了全局最低能量，更不意味着 18 个都能结合。构象数量、收敛状态和后续姿势质量分别回答覆盖、计算状态和任务效果，验收时要把三者连接起来。

## 4. 一次运行完整案例

完整实现保存在 [run_rdkit_case.py](../../scripts/run_rdkit_case.py)。它还为每个构象保存 ID、收敛状态和 MMFF 能量，并生成后面的肽性质表。

```bash
python scripts/run_rdkit_case.py \
  --input workflows/gnina/example_1iep/reference.sdf \
  --output runs/09_my_rdkit
```

输出目录必须是新目录。完成后先用 `run_summary.json` 核对版本、输入和电荷，再逐行检查 `conformer_energies.csv` 的构象 ID 与收敛状态码。随后读回 `imatinib_conformers.sdf`，确认各条结构保留同一化学身份；最后核对 `peptide_properties.csv` 是否包含 WT、D5E、S4E、L9K。这个顺序从输入身份走到输出数量，能区分化学问题、采样问题与文件保存问题。

```python
import csv
from rdkit import Chem

path = 'runs/09_my_rdkit'
poses = [m for m in Chem.SDMolSupplier(
    f'{path}/imatinib_conformers.sdf', removeHs=False
) if m is not None]
with open(f'{path}/conformer_energies.csv') as stream:
    rows = list(csv.DictReader(stream))
assert len(poses) == len(rows) > 0
assert all(int(row['mmff_status']) == 0 for row in rows)
print(len(poses))
```

## 5. 本地实际得到什么

2026-10-23 在 Debian 使用上述版本完成的记录见 [运行摘要](../../examples/debian13_20261023/rdkit/run_summary.json)和[构象能量表](../../examples/debian13_20261023/rdkit/conformer_energies.csv)。

本轮请求 20 个构象，实际保存 18 个，18 个均达到 MMFF 收敛条件。去重和嵌入成功数共同决定输出规模，因此验收应按真实保存数量进行。各构象共享 37 个重原子和 +1 的形式电荷，与本课参考配体的化学身份一致。

当前分子表示的分子量为 494.623 Da，RDKit 计算 logP 为 3.173，TPSA 为 87.48 Å²，氢键供体和受体数分别为 3 和 6。它们是按当前化学结构计算的描述符，用来说明组成和理化特征。规范 SMILES 往返后的 Morgan 相似度为 1.0，说明本次往返保留了这套指纹所表达的特征。


**logP** 描述中性分子在两相中的分配倾向；对这份带电表示算出的 RDKit 描述符，应连同输入状态报告。实验特定 pH 下的 logD 还包含电离平衡。**TPSA** 根据原子类型估计极性表面积，常用于比较分子的极性特征；跨膜和口服行为还受到电荷、构象和转运等因素影响。

**Morgan 指纹**把原子周围的局部连接环境编码成特征；本课使用半径 2、2048 位。**Tanimoto 相似度**比较两个指纹共有的特征比例。本课的 1.0 是一次身份往返检查；实际去重还应比较规范异构 SMILES、手性及盐型处理规则，以免位指纹碰撞或参数省略造成混淆。

指纹相似度可以用集合的交并比理解：共有特征占全部出现特征的比例越大，Tanimoto 越接近 1。但位指纹把许多化学环境压缩进有限位置，存在碰撞；指纹设置还可能忽略某些手性信息。因此，相似度 1.0 适合作为往返检查的一项证据，严格身份核对应同时查看规范异构 SMILES、手性、电荷与标准化策略。

## 6. 把同样的方法接到序列优化

第六章提名了四条肽。现在用 `MolFromFASTA` 建立它们的线性、游离端基、标准 L 氨基酸化学结构：

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

peptide = Chem.MolFromFASTA('ETFSDLWKLLPEN')
if peptide is None:
    raise ValueError('Invalid peptide sequence')
print(rdMolDescriptors.CalcMolFormula(peptide))
print(Descriptors.MolWt(peptide))
print(Chem.GetFormalCharge(peptide))
```

WT `ETFSDLWKLLPEN` 的分子量为 1591.783 Da。D5E 得到 `ETFSELWKLLPEN`，分子量为 1605.810 Da，Asp 改为 Glu 后侧链增加一个亚甲基；S4E 得到 `ETFEDLWKLLPEN`，分子量为 1633.820 Da，Ser 被较长的酸性侧链取代。

L9K 的序列为 `ETFSDLWKKLPEN`，分子量为 1606.798 Da，变化是把 Leu 换成可电离的 Lys。先把这些修改落实到明确的化学结构，才能继续讨论电荷、端基和制备条件。


本地构建后的四个 `Mol` 形式电荷都是 0，因为该 API 生成的结构采用其默认电荷表示。第六章基于侧链 pKa 估算的 pH 相关净电荷回答另一项问题：溶液条件下各可电离基团大致处于什么状态。解释两张表时要说明建模方法和 pH。

如果项目要求 N 端乙酰化、C 端酰胺化、环肽、D 氨基酸或非天然残基，应在明确的化学结构上进行编辑，并重新检查手性、端基和键连接。序列字符串在这些任务中需要配套的化学修饰记录。

形式电荷是当前结构图上写下的整数电荷，pH 相关净电荷则是对溶液中质子化平衡的估计。肽的默认结构图形式电荷为零，与某个 pH 下预计带负电并不矛盾；两者采用的描述层次不同。后续若要模拟或对接，必须选择具体质子化状态，不能把一个平均净电荷数字直接分配给每个原子。

## 7. 常见卡点与处理

读取得到 `None` 时，先查价态、键级和文件完整性，定位具体报错分子。不能嵌入或 MMFF 参数不全时，检查金属、非常规元素及化学状态，确定当前构象和参数方法是否适用。

构象数少于请求值时，应查看去重阈值与嵌入成功数，报告实际数量；同一分子去重失败时，先统一盐、互变异构、质子化和手性策略。跨机器略有差异时，把 RDKit 版本、种子与线程数一起记录，再比较构象覆盖和能量分布，避免只盯着某一个编号。


## 8. 面试时把这件事讲清楚

回顾本课，输入准备应从结构来源、键级、手性和质子化状态开始，保留晶体配体作为参考，再从连接关系生成独立构象。参数覆盖、优化状态和文件读回一起构成输入验收。本轮 18 个收敛构象提供了后续重对接的合理起点，是否恢复已知姿势则由 GNINA 结果评价。

这条操作线也适用于多肽：先把序列落实为明确化学结构，再讨论分子量、端基和修饰。结构保持与结合问题继续交给 [PyRosetta](../10_pyrosetta/tutorial.md)、复合物预测及实验。把每项工具承担的任务衔接好，才能从“文件能读取”走向“结果能解释”。

练习交付：把 `numConfs` 改成 50，使用新输出目录运行，比较保存数量和能量分布；再解释为什么去重阈值会影响构象库的规模。优先掌握“化学身份→合理构象→可追踪文件”这条操作线。
