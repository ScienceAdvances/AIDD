# Boltz：做复合物预测，正确使用 Boltz-2 亲和力输出

> 本文对应网站教案的 Debian 版本。仓库练习请先读同目录 [README](README.md)：其中命令使用本仓库目录和可配置工具路径；正文的运行命令已按仓库路径调整；多工具任务统一通过可配置入口执行。


这一课包含两项任务：把 MDM2 和 p53 两条序列交给 Boltz-2 做复合物预测；再用官方蛋白—小分子示例学习亲和力输入和输出单位。两项任务共享同一个程序，输出的科学含义分别解释。

学完后，应能制作 YAML、准备 MSA、运行多样本预测、检查链间置信度，并把亲和力值正确换算成 IC50 的模型预测尺度。面试中尤其要说明：结构置信度、结合分类概率和活性数值各自回答什么问题。

## 1. Debian 的 Boltz 环境与模型

本轮使用现有 `/opt/micromamba/envs/boltz`，Python 3.11、Boltz 包 2.2.1、Torch 2.6.0+cu124、RDKit 2026.9.1，选择 Boltz-2 模型。

```bash
cd .
python scripts/run_boltz_compatible.py predict --help
```

2.2.1 的亲和力解析调用 `AllChem.Descriptors.MolWt`，本环境 RDKit 的公开模块为 `rdkit.Chem.Descriptors`。兼容入口将该公开模块挂到预期命名空间，再进入原版 CLI，模型与分子量实现保持原样。

## 2. 第一份输入：两条蛋白链

教材保存了 [human_mdm2_p53_no_msa.yaml](../../workflows/modern_tools/boltz/human_mdm2_p53_no_msa.yaml)。内容与 AlphaFold 3 章使用同一组序列：

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: ETLVRPKPLLLKLLKSVGAQKDTYTMKEVLFYLGQYIMTKRLYDEKQQHIVYCSNDLLGDLFGVPSFSVKEHRKIYTMIYRNLVV
      msa: empty
  - protein:
      id: B
      sequence: ETFSDLWKLLPEN
      msa: empty
```

`id` 为预测中的链标识；序列分别为 85 和 13 个残基。`msa: empty` 明确采用无 MSA 模式，便于离线走通流程。第十三章在这一信息设置下没有可靠恢复参考；由于本轮缺少同模型有无 MSA 的成对对照，尚不能把差异单独归因于 MSA。

要使用 MSA 服务，另存一份 YAML，删除两行 `msa: empty`，并在命令中加 `--use_msa_server`。这会将待搜索序列提交给所配置的 MSA 服务。项目若已有本地 A3M，则按官方格式提供路径并保存比对来源；同一组候选尽量使用一致的信息准备策略。

## 3. 资源与运行前检查

官方分子库及两份模型已下载到 `downloads/boltz/`。结构与亲和力权重分别约 2.29GB、2.06GB，按发布的 SHA256 验证，见[权重清单](../../examples/debian13_20261023/cache/boltz/verified_weights.json)。

先核对 YAML 两条链的长度与空 MSA 设置，再检查输出中的 processed 输入、成功预测数量和实际权重。亲和力案例还需检查配体身份及 `properties`。解析成功只表示输入可展开，结构与亲和力推理仍要分别完成。

## 4. 运行复合物预测

首次运行会下载分子库和权重，2.2.1 的资源准备还涉及亲和力模型，预留数 GB 缓存空间。`--cache` 固定缓存位置，方便复用和记录来源。

```bash
python scripts/run_boltz_compatible.py predict \
  workflows/modern_tools/boltz/human_mdm2_p53_no_msa.yaml \
  --model boltz2 \
  --out_dir runs/14_my_structure \
  --cache downloads/boltz \
  --accelerator gpu --devices 1 --max_parallel_samples 1 --no_kernels \
  --seed 20261020 --num_workers 0 \
  --recycling_steps 3 --sampling_steps 200 --diffusion_samples 5 \
  --write_full_pae --write_full_pde
```

小体系 CPU 练习可以把 `--accelerator gpu` 改成 `cpu`。第一次先使用 `--diffusion_samples 1` 确认模型能完整运行，再恢复五样本。显存紧张时可减小 `--max_parallel_samples`，保持总采样数，并记录改变后的设置。

`recycling_steps` 控制表示迭代，`sampling_steps` 控制每次扩散采样的步数，`diffusion_samples` 控制结构样本数量。三者分别影响计算过程，比较候选时将其固定，才能判断分数差异主要来自候选还是计算规模。

过程验收依次看输入预处理、MSA 状态、权重加载、成功预测条目数和输出文件。遇到失败先检查日志与 `processed/` 中的输入记录，再定位模型阶段。

## 5. 读取五份结构与置信度

本轮空 MSA 的 MDM2/p53 五样本全部完成。按置信度最高选择的 model_0：iPTM=0.956，MDM2 CA RMSD=0.549 Å，受体对齐后 p53 CA RMSD=2.052 Å。五样本 p53 位移范围为 1.79–2.43 Å。

[逐样本表](../../examples/debian13_20261023/boltz_structure_comparison.csv)与[分析摘要](../../examples/debian13_20261023/boltz_analysis_summary.json)保存全部结果。叠合方法与 AF3 相同：先核对链身份，再只按 MDM2 CA 对齐，随后观察 p53。

这一案例中 Boltz 恢复的摆放更接近参考；不能由一个小体系推导普遍的模型优劣。高置信度和接近参考的结构也不等同于实验活性。

本轮五个样本的 iPTM 约为 0.935–0.956，模型对链间安排的置信度较高；按受体对齐后，p53 RMSD 为 1.79–2.43 Å，则是事后对参考恢复程度的检查。选择置信度最高的 model_0，是预先可执行的选择规则；事后选择 RMSD 最小的样本，则使用了参考答案。两者可以同时报告，但不能把后者写成模型在未知参考时的实际首选表现。

它与 AF3 的结果形成一个具体对照：相同序列和空 MSA 下，两套模型在这一体系的参考恢复不同。模型架构、权重与采样实现同时变化，因此这个对照不能确定是哪一项造成差异，也不能代表其他靶点的总体准确率。

## 6. 第二份输入：官方蛋白—小分子亲和力案例

[official_affinity.yaml](../../workflows/modern_tools/boltz/official_affinity.yaml)原样来自[官方示例](https://github.com/jwohlwend/boltz/blob/v2.2.1/examples/affinity.yaml)。A 链为官方给定的 384 残基蛋白，B 链由一个带手性的氨基酸小分子 SMILES 描述。关键片段是：

```yaml
  - ligand:
      id: B
      smiles: 'N[C@@H](Cc1ccc(O)cc1)C(=O)O'
properties:
  - affinity:
      binder: B
```

上面片段应放在官方完整 YAML 的对应位置。`properties.affinity.binder` 指向要评价的小分子链。当前亲和力任务面向蛋白—小分子；p53 肽和 VHH 使用蛋白链建模时，按结构置信度及相应实验终点评价。

```bash
python scripts/run_boltz_compatible.py predict \
  workflows/modern_tools/boltz/official_affinity.yaml \
  --model boltz2 \
  --out_dir runs/14_my_affinity \
  --cache downloads/boltz \
  --accelerator gpu --devices 1 --max_parallel_samples 1 --no_kernels --num_workers 0 \
  --seed 20261020 --use_msa_server \
  --diffusion_samples 5 --write_full_pae
```

官方文件没有提供预计算 MSA，因此这里启用服务。若改用离线 MSA，先在 YAML 中准备好相应比对路径，再运行。完成后查找 `affinity_*.json`，并保留对应复合物结构，检查配体化学身份与预测口袋。

`binder` 在输入中指明要评价哪个小分子实体，是任务选择字段，并不表示这个分子已经被实验证实会结合。结构输入、亲和力属性和输出文件需要逐项对应；尤其不能把这套蛋白—小分子输出的终点迁移到以 protein 链表示的 p53 肽或 VHH。

## 7. 亲和力单位是本课必须记住的内容

根据本版[官方输出定义](https://github.com/jwohlwend/boltz/blob/v2.2.1/docs/prediction.md#output)，`affinity_pred_value` 的尺度为：

\[
y = \log_{10}\bigl(IC_{50}\,[\mu M]\bigr)
\]

数值越低，模型预测的活性越强。由它换算得到 `IC50_µM=10**y`、`IC50_nM=1000*10**y`，以及无量纲的 `pIC50=6-y`。

可以先做三个单位算例：y=−3 时，IC50 为 10⁻³ µM，即 1 nM，pIC50 为 9；y=0 时为 1 µM，pIC50 为 6；y=2 时为 100 µM，pIC50 为 4。注意这里 y 越低越强，而 pIC50 越高越强，两种尺度的方向正好相反。


这些数值用于解释单位，独立模型推理结果应从自己的 JSON 读取。**IC50** 是特定实验条件下达到半数抑制的浓度；**Kd** 描述平衡结合。报告时保留输出终点名称和单位，实验比较再考虑测定条件。

另一个输出 `affinity_probability_binary` 位于 0–1，面向结合与非结合分类。官方将分类和数值亲和力用于不同阶段：分类帮助初筛，数值头更适合比较已经具有活性的分子优化。实际任务应先建立阳性、阴性和独立验证集，检查模型在目标体系中的区分能力。

```python
import json
from pathlib import Path

for path in Path('runs/14_my_affinity').rglob('affinity_*.json'):
    result = json.loads(path.read_text())
    y = result['affinity_pred_value']
    print(path.name, {'predicted_ic50_nM': 1000 * 10**y,
                      'predicted_pIC50': 6-y,
                      'binder_probability': result['affinity_probability_binary']})
```

分类头与数值头并行提供信息，却不是同一个量的两种单位。分类概率为 0.3534，不是 IC50=0.3534 µM，也不是抑制率 35.34%。数值头输出 y=2.2371197，才按本版定义换算为约 172.63 µM。实际使用时应先检验两种输出在相关阳性、阴性和活性系列里的表现，再决定它们如何参与提名。

IC50 的数值还依赖具体实验的底物浓度、机制和测定条件；Kd 则描述平衡结合。两者在满足特定假设时可能建立关系，当前通用模型输出没有提供足够条件完成这种转换。因此，报告时保留“模型预测 IC50 尺度”，比把它改名为“结合 Kd”更准确。

## 8. 本轮亲和力结果及模型间比较

官方蛋白—小分子案例使用 MSA server，实际完成结构与亲和力预测。`affinity_pred_value=2.2371197`，按本版输出定义换算为预测 IC50 **172.63 µM**，结合分类概率 **0.3534**。这是独立于 MDM2/p53 的输入，不能把该值当作 p53 肽的亲和力。

[原始 JSON](../../examples/debian13_20261023/boltz_affinity/boltz_results_official_affinity/predictions/official_affinity/affinity_official_affinity.json)保存原始值，[换算摘要](../../examples/debian13_20261023/boltz_analysis_summary.json)保留单位与来源。本轮结构使用空 MSA，亲和力案例使用服务器搜索；比较时要说明这一区别。

与第十五章衔接时，再检查 BoltzGen 的复合物 RMSD 和实际筛选标志。结构置信度、相对输入的几何保持和模型活性预测各自回答不同问题。

## 9. 排错与面试练习

NumPy 或 Torch 依赖冲突时，先核对本轮 Python 3.11 环境与已锁定依赖；若确需另建环境，再按既定边界确认。亲和力解析报 `AllChem.Descriptors` 时使用教材兼容入口，并记录该适配。MSA 缺失或服务失败时，检查 YAML 的 `msa` 字段和服务连接。

没有亲和力文件时，应核对 `properties.affinity.binder`、链类型及亲和力阶段日志，不能只凭结构文件出现就认为两项推理都完成。显存不足时先检查 token 规模并减少并行样本，保留总样本数；不同指标排名冲突时，回到结构置信度、分类概率和数值活性各自的终点。


实际项目中，可以先用已知阳性、阴性和活性系列检验 Boltz 输出在目标体系里的意义，再以一致的信息准备和采样协议评价新候选。结构置信度用于检查复合物假设，分类概率用于辅助区分候选，数值头保留其预测 IC50 尺度；三者不能互相替代。

比较 AF3、Boltz 和 GNINA 时，先统一分子体系与任务。复合物预测关注结构恢复和链间安排，已知口袋对接关注姿势搜索与排序，活性预测则需要独立实验数据评价。把任务说清楚，工具之间的比较才有可解释性。

练习交付：完成两份 YAML 的解析，再在 GPU 环境分别运行结构与亲和力案例，交付原始 JSON、结构、版本及单位换算表。把模型原始值和换算后的值同时保留。
