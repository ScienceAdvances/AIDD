# GNINA 小分子对接实战：输入制作、安装、运行与重对接评价

> 本文对应网站教案的 Debian 版本。仓库练习请先读同目录 [README](README.md)：其中命令使用本仓库目录和可配置工具路径；正文的运行命令已按仓库路径调整；多工具任务统一通过可配置入口执行。


前面的 HADDOCK 案例处理蛋白与肽。现在换成小分子：伊马替尼与 c-Abl。我们使用 AutoDock Vina 官方教程的 1IEP 输入，把它组织成 GNINA 工作流。

这里先做一个有答案的问题。实验已经告诉我们共晶配体位于哪里，我们另给待对接配体生成起始构象，再让软件尝试把它放回口袋。这个练习叫**重对接**。它帮助检查输入、搜索与排序能否恢复已知姿态，是建立对接流程时有用的一步。

本章已在 Debian GPU 完成六项真实对接，保留输入、日志、姿势和 RMSD。下面从输入与小规模检查开始，再解释完整对照结果。

## 1. GNINA 在整件事中负责什么

小分子对接要同时寻找配体在口袋里的位置、朝向和内部扭转角。本例的受体保持固定，配体仍可调整内部构象；“受体刚性”不表示配体是一块完全不能弯曲的物体。HADDOCK 的刚体阶段则主要先改变两组分的相对平移和旋转，随后再做柔性精修。二者都处理空间搜索，但分子类型、自由度和评价协议不同。

本例通过共晶配体知道口袋在哪，因而属于带已知位置先验的重对接练习。它检验在指定受体与化学状态下能否找回已知姿势；若目标是从化合物库发现新活性，还需加入活性与非活性分子、不同受体构象及独立实验评价。


小分子对接需要改变配体的平移、旋转和可旋转键，搜索可能的结合姿势，再给姿势排序。GNINA 源于 smina，后者源于 AutoDock Vina，并加入卷积神经网络，即 CNN，用于结构相关的评分等步骤。[官方实现](https://github.com/gnina/gnina)

评价这套方法时，需要选取相关体系的验证任务，检查输入质量、姿势恢复与排序表现。本章面向小分子–蛋白任务；VHH–抗原则需要采用适合蛋白复合物的表示和搜索方法。

本章先让传统经验评分参与对接，再用 CNN 对输出姿势重评分。你可以把搜索理解为“提出可能的摆法”，把重评分理解为“用另一个模型重新评价这些摆法”。重评分的候选范围由搜索阶段提供，因此搜索是否覆盖合理姿态会直接影响后续排序。

## 2. 三个文件，三个不同用途

进入教材的 `workflows/gnina` 目录，可以看到准备好的 1IEP 输入：

`receptor.pdb` 提供受体坐标，`reference.sdf` 保留共晶配体的位置，用于定义口袋和事后评价；真正进入搜索的 `ligand_input.sdf` 则采用独立生成的起始构象。先区分这三份文件，才能理解程序是在恢复已知姿势，还是仅仅沿用参考坐标。


参考配体保留与受体一致的共晶坐标，用于定位口袋和评价结果；待对接配体使用独立起始构象，让搜索从新的构象出发恢复结合姿势。

这份待对接配体已核对化学身份，有 37 个重原子，准备记录中的形式电荷为 +1。这个记录明确了当前输入的化学状态。换到新化合物时，还要检查键级、立体化学、质子化与互变异构等问题。

输入来源是 [Vina 官方基础案例](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html)。我们保留这套官方输入，另按 GNINA 接口编写运行与分析流程。

## 3. 搜索盒子为什么要围着口袋设置

对接需要先界定有限的搜索区域。搜索盒子规定允许尝试的位置范围；太小可能挡住正确姿态，太大则把有限计算花在许多不相关区域。

这里已经知道共晶配体的位置，所以用 `reference.sdf` 定义口袋，再在周围加缓冲。`--autobox_ligand` 根据所给文件的坐标范围建立盒子；本例的口袋位置由共晶配体提供。

如果把整份受体传给它，得到的会是覆盖受体范围的搜索区域。这种设置对应全蛋白范围的搜索。已知口袋的重对接和未知口袋搜索，是不同难度的问题。

## 4. Debian GPU 引擎与分析环境

本轮使用服务器已有 `/opt/bioinfra/bin/gnina`（1.3.3）与 RTX 3060。Python/RDKit 使用 `.venv-analysis`。从教材目录检查：

```bash
cd .
# 如需额外动态库路径，在个人工具配置的 extra_env 中指定。
gnina --version
python workflows/gnina/gnina_workflow.py --help
```

本轮初次启动缺少 cuDNN 动态库，按上面方式为当前进程指定现有环境的库目录后完成 GPU 计算。没有通过新建引擎环境解决这一问题。

## 5. 输入怎样制作：取官方结构，再生成独立配体构象

第一次先使用随教材保存的 `example_1iep/`。它有三份坐标和一份输入清单：

`receptor.pdb` 采用 Vina 官方 1IEP 教学受体，已经准备氢，其坐标必须与共晶参考处于同一框架。`reference.sdf` 保存官方共晶配体的实验位置，定位和评价时保留原坐标；`ligand_input.sdf` 依据同一化学身份重新生成构象，并记录键级、手性与形式电荷。`input_manifest.json` 随后把下载来源、指纹、RDKit 版本和生成参数联系起来，说明每份输入怎样得到。


想亲自从来源重新准备时，在 GNINA 目录执行：

```bash
python gnina_workflow.py fetch-example --work my_inputs
python gnina_workflow.py check \
  --sdf my_inputs/reference.sdf my_inputs/ligand_input.sdf
```

`fetch-example` 依次下载官方受体和配体，读取参考配体的化学身份，保留参考坐标，再为待对接配体重新生成起始构象：去除旧构象，补氢，使用 **ETKDGv3** 嵌入，再以 **MMFF** 做构象优化。随机种子为 `20261011`，优化收敛状态写进清单。

ETKDG 利用距离几何和经验结构知识产生合理三维起点；MMFF 用小分子力场调整局部几何。它们准备的是供搜索使用的起始构象。配体在蛋白口袋中的具体位置由后续对接探索。

这份案例输入应有 **37 个重原子，形式电荷 +1**，参考与待对接配体的化学身份一致。脚本沿用官方化学状态；为新化合物准备输入时，再按实验 pH、手性、互变异构与共价状态审查候选形式。

受体准备也需要具体判断。本例沿用官方已准备受体；换新 PDB 时，先选生物学相关链和口袋，再检查缺失原子、替代构象、质子化、金属、辅因子与结构水。每项处理都记录到自己的输入清单。

## 6. 把输入变成任务表，先验收再运行

以下主线继续使用 `example_1iep/`，便于对照当前教材的输入指纹。先检查化学结构并运行四项辅助测试：

```bash
python gnina_workflow.py check \
  --sdf example_1iep/reference.sdf example_1iep/ligand_input.sdf
python -m unittest -v test_workflow.py
```

测试覆盖同坐标系 RMSD、缺失分数和输出文件完整性等行为。然后打开 `jobs_1iep.tsv`，它有五列：任务名、受体、待对接配体、参考配体、随机种子。自己写一份 TSV 时，用 Python 明确保存制表符：

```bash
python - <<'PY'
import csv
from pathlib import Path

with Path('my_jobs.tsv').open('w', newline='') as stream:
    writer = csv.writer(stream, delimiter='\t')
    writer.writerow(['job_id', 'receptor', 'ligand', 'reference', 'seed'])
    for seed in (11, 23, 47):
        writer.writerow([f'abl_imatinib_seed{seed}',
                         'example_1iep/receptor.pdb',
                         'example_1iep/ligand_input.sdf',
                         'example_1iep/reference.sdf', seed])
PY
python gnina_workflow.py run \
  --jobs my_jobs.tsv --output my_dry_run \
  --gnina gnina --device 0 --cpu 4 --exhaustiveness 32 --num-modes 20 --dry-run
```

相对路径按 **TSV 所在目录**解释。若使用第 5 节新准备的输入，将三列文件路径改为 `my_inputs/...`。三个任务使用相同分子与参数，仅改变随机种子，用于观察搜索的重复性。

预演会解析每份输入、检查任务名与文件路径，并在 `run_manifest.json` 保存实际命令和输入指纹。确认这三项任务都生成后，再开始引擎计算。

## 7. 真正运行：先小规模检查，再增加搜索投入

在 Linux 的 GNINA 目录，先运行一个采样较小的任务：

```bash
gnina \
  --receptor example_1iep/receptor.pdb \
  --ligand example_1iep/ligand_input.sdf \
  --autobox_ligand example_1iep/reference.sdf --autobox_add 4 \
  --scoring vina --cnn_scoring rescore --pose_sort_order CNNscore \
  --cpu 4 --device 0 --seed 11 --exhaustiveness 2 --num_modes 5 \
  --out my_smoke.sdf > my_smoke.log 2>&1
```

`--autobox_ligand` 从参考配体坐标定义搜索盒，`--autobox_add 4` 在边界外加 4 Å 缓冲；实际尺寸还受 GNINA 的扩展设置影响，运行后查看日志。输入清单根据参考边界算出的中心约为 `(15.107, 53.977, 17.143)` Å。

`--scoring vina` 选择经验评分，`--cnn_scoring rescore` 对搜索输出进行 CNN 重评分，`--pose_sort_order CNNscore` 据此排序。最后一组参数指定四个 CPU 辅助线程、GPU 0、随机种子、搜索投入和最多五个姿势。

先查看 `my_smoke.log` 中的启动信息、盒子、搜索与评分，再确认 SDF 可读取：

```bash
python gnina_workflow.py check --sdf my_smoke.sdf
```

小规模任务成功后，使用任务表执行正式的三个种子：

```bash
python gnina_workflow.py run \
  --jobs my_jobs.tsv --output my_cnn_run \
  --gnina gnina --device 0 --cpu 4 --exhaustiveness 32 --num-modes 20 --cnn-scoring rescore
```

这次输出保存在 `my_cnn_run/`。辅助程序逐个任务执行，保存版本、完整参数、耗时、退出状态、输入指纹和有效输出数。每个 SDF 可以包含多个姿势；日志与清单用同一个任务名对应。

## 8. 运行后怎样形成一张可以判断的结果表

先提取每个姿势的真实分数，再计算各任务与共晶参考之间的 RMSD：

```bash
python gnina_workflow.py summarize --run my_cnn_run --output my_scores.csv
python - <<'PY'
import subprocess
import sys

for seed in (11, 23, 47):
    subprocess.run([
        sys.executable, 'gnina_workflow.py', 'rmsd',
        '--reference', 'example_1iep/reference.sdf',
        '--poses', f'my_cnn_run/abl_imatinib_seed{seed}.sdf',
        '--output', f'my_rmsd_seed{seed}.csv',
    ], check=True)
PY
```

姿势比较保持共同受体坐标框架。程序按化学身份建立对应原子映射，并在有需要时考虑对称等价；复杂化学仍需检查映射是否合理。

结果至少整理成下列字段：种子、实际姿势数、排序第一的 RMSD、所有输出中的最佳 RMSD、第一名的经验分数、CNNscore、CNNaffinity 和耗时。Top 1 描述最终选择的表现，最佳输出 RMSD 描述搜索集合是否包含接近参考的解。

本轮已经完成两个评分方案、三个种子共六次真实引擎运行，每次输出 20 个姿势。结果见下面的实测记录。

### 8.1 Debian GPU 对照结果

三种子为 11、23、47；两种方案为 `cnn_scoring=rescore` 与 `none`，每次 exhaustiveness=32、num_modes=20。三种子的首位重原子 RMSD 分别约 0.786、0.793、0.790 Å。该组输入下，CNN 重排序与经验能量排序的第一名一致。

[六任务对照表](../../examples/debian13_20261023/gnina_comparison.csv)保留第一名、最佳输出与真实分数，各姿势原始 SDF 在 `gnina_rescore_gpu/` 和 `gnina_none_gpu/`。本例找回了接近参考的姿势，但没有证据说明 CNN 在此例优于能量排序。

复现本轮完整对照可运行：

```bash
# 先在 config/tools.json 配置本章所用环境与模型资源
python practice.py run 07 --output runs/07_new --execute
```

## 9. 输出里的三类分数，分别在评价什么

读结果时，先区分三类分数的定义与用途。

传统经验分数通常以 kcal/mol 表示，在可比条件下较负往往更受该函数偏好，其物理近似会影响它与实际结合自由能的关系。`CNNscore` 侧重模型对姿势的评价，通常优先检查较高分的姿势。`CNNaffinity` 是亲和力相关的预测值，按 pK 尺度表达，需要用相关实验数据评价预测误差。

例如在实验 pKd 的定义下，pKd=6 对应 1 μM，pKd=9 对应 1 nM。报告时应注明数值来自模型预测还是实验测定，并给出相应的误差信息。定义与方法见 [GNINA 方法论文](https://link.springer.com/article/10.1186/s13321-021-00522-2)。

如果经验分数很负，CNNscore 却很低，可以回到配体身份、电荷、口袋、碰撞和参考姿态，检查分歧来自什么。能把冲突讲清楚，比只报告一个最低分更有价值。

例如，经验分数 −12.9 比 −10 更受当前函数偏好；CNNscore 接近 0.92 表示模型对该姿势给出较高评价。这两个数不能相互相减，也不能把 0.92 直接读成“实验结合成功率 92%”。CNNaffinity 即使采用 pK 尺度，也仍是训练数据与模型条件下的预测，需要在相关体系上检验校准和误差。

如果某条姿势的几项分数同时较好，可以提高检查它的优先级，但还应查看是否存在不合理接触、化学状态错误或口袋先验偏差。评分的用途是帮助排列待验证假设，具体结论要由对应的几何或实验终点支持。

## 10. 在共同受体坐标系中评价姿势

重对接有共晶答案，因此可以计算配体重原子 RMSD。这里要比较的是配体在受体上的位置与形状，所以受体坐标框架必须一致。

想象一个配体形状完全正确，却整体偏离口袋 3 Å。如果先把它单独平移回参考位置，再算 RMSD，就会得到很小的误差。你检查了内部形状，却把真正需要评价的相对位置误差消除了。

本章辅助程序不对配体单独拟合，并检查化学身份和原子映射；复杂对称性仍要人工审查。对应的三种子提取命令见第 8 节。

同时关注排序第一的 RMSD 和所有输出中的最小 RMSD。前者反映搜索与排序合起来的效果，后者告诉你输出集合里有没有接近参考的解。这与 HADDOCK 中分开检查集合最佳 DockQ 和评分第一名 DockQ 的思路相通。

教学上常用 2 Å 一类阈值判断是否接近参考，但它用于衡量姿势恢复。新化合物的结合与细胞药效分别通过相应实验评价。

同一个配体的对称原子可能有多种合理对应，因此 RMSD 还取决于原子映射。映射错误会把化学等价姿态判成很差，单独拟合配体又会把真正的摆放误差消除。正确做法是先保证化学身份和映射合理，再在共同受体坐标框架里比较。低 RMSD 只说明恢复了参考姿势；参考结构本身、蛋白构象和实验条件仍构成这项评价的边界。

## 11. 设计一次 GNINA 内部对照

本轮已经完成三个种子、两种评分模式共六次计算，每次输出 20 个姿势。复查时分别记录 Top 1 与最佳输出 RMSD、姿势数和耗时，再检查结构。下面是关闭 CNN 的内部对照命令，输入、盒子和采样设置保持一致：

```bash
python gnina_workflow.py run \
  --jobs my_jobs.tsv --output my_energy_run \
  --gnina gnina --device 0 --cpu 4 --exhaustiveness 32 --num-modes 20 \
  --cnn-scoring none
```

这项对照比较 GNINA 内部两种评分模式；与原版 Vina 比较时，需要另外准备和记录 Vina 工作流。如果要隔离“排序规则”的作用，还应冻结同一组姿势再分别评分；两次运行输出的第 N 条，不一定是同一个结构。

本轮[对照记录](../../examples/debian13_20261023/gnina_comparison.csv)已包含真实分数与 RMSD。先定义共同受体坐标框架及比较标准，再解释各个种子的结果。

本轮 `rescore` 和 `none` 两种模式在三个种子下选到的首位姿势相同，对应 RMSD 约为 0.786、0.793、0.790 Å。这说明在这个已知口袋案例里，两种模式都恢复了较接近参考的首位姿势；当前结果没有显示 CNN 带来了首位选择优势。

需要进一步区分“改了整个搜索流程”与“只改排序规则”。两次独立运行可能输出不同姿势集合；若要隔离排序的贡献，应冻结同一组结构，再分别计算分数并比较其排序。保存的 20 个姿势也是经过筛选的输出，不代表搜索曾访问过的所有构象。

## 操作中遇到的问题

程序启动失败时先核对 Linux、CPU 架构、动态依赖和发布包；输入下载或解析失败时检查文件完整性、SDF 化学结构及构象。若 TSV 中路径找不到，应按 TSV 所在目录解释相对路径，再检查制表符和文件名。

输出某种分数为空时，回到 SDF 属性名、CNN 模式和日志，保留缺失状态。RMSD 异常时则先核对化学身份、原子映射、共同受体坐标框架和口袋位置，确保比较的确实是同一分子的同一组原子。


## 读完后，试着自己解释

如果有人问“GNINA 比 Vina 好在哪里”，先说明它增加了 CNN 相关能力，再说明你的判断需要在哪个体系、哪种评价任务和哪个版本下验证。把结论对应到具体验证任务，才能说明该能力在什么条件下有帮助。

再用一句话解释共晶参考配体与待对接配体各自的作用。能说清两个文件的用途，就理解了这个案例的起点。最后读[综合课](../05_integration/tutorial.md)，把不同工具的结果连成合理的研究决策。
