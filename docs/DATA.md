# 输入、示例与可选数据包

`data/` 保存原始小型输入，`workflows/` 保存本课使用的约束和任务配置，`examples/` 保存 Debian 参考结果，`runs/` 供读者产生新结果。共享数据采用一份副本，各章通过 README 指向它，避免复制后版本不同。

## 随仓库提供的材料

多肽方向提供 1YCR 原始结构、小鼠 p53 序列、HADDOCK 输入结构与 AIR、固定版 ProteinMPNN 科学代码，以及六条 Debian 设计的解析、概率和序列输出。序列优化保留 197 条候选、八组评分审计 NPZ 和全部提名标记。

CAR-T 提供论文公开的原始工作簿与派生表，能够从单元格重新走到 204 个观测、49 个分组及三张图。小分子方向提供 1IEP 已准备受体、共晶参考、独立起始配体及三个种子的任务表；两种评分模式各 20 个姿势的真实输出也已保留。

MD 提供建模所需 PDB/MDP、逐帧结构和热力学表，以及一个恢复完整并居中的蛋白轨迹 `protein_centered.xtc` 和对应参考 PDB。这个小轨迹覆盖 201 帧，适合练习蛋白结构分析；它不包含全部水和离子，也不能替代含水原始轨迹来做溶剂分析。

进阶案例提供 RFdiffusion2 汇总指标、RFantibody 三阶段 QV 和最终结构、AF3 五样本原始输出、Boltz 的结构/亲和力原始结果，以及 BoltzGen 最终排名与过滤字段。没有新增结合、细胞或酶活性实验。

## 按需获取的资源

```bash
python practice.py fetch list
python practice.py fetch proteinmpnn_weights
python practice.py fetch gromacs_reference
```

下载器先核对保存的大小和 SHA256，再安装到目标目录；已存在而内容不同的文件不会被静默覆盖。`gromacs_reference` 下载的是官方原始参考轨迹，并非作者 Debian 新跑的轨迹。官方 main 分支归档若更新，旧哈希校验会失败，应回到来源确认，而非直接移除校验。

其他模型按各自官方方式获取。AF3 参数、PyRosetta 安装包和完整模型缓存没有放进 Git；运行相应工具时可能需要读者确认许可、预先准备资源，或允许工具按官方方式下载。

## 可选 Release 数据

维护者可将单独交付的 `AIDD-practice-data-debian13.tar.gz` 上传到 GitHub Releases。它包含作者 Debian 的完整 MD 运行文件、HADDOCK smoke/full 原始运行目录，以及 RFdiffusion2 主链、映射和序列/Chai 产物；不包含模型缓存或 Python 环境。该包独立于基础仓库，大小与 SHA256 见 `provenance/release-assets.json`。

下载后在仓库根目录解包：

```bash
tar -xzf AIDD-practice-data-debian13.tar.gz
```

解包到 `optional_data/debian13_20261023/`。完整 MD 位于其 `md/` 子目录，可传给 `scripts/analyze_md.py --local-run`；完整 HADDOCK 在 `haddock/runs/`；RFdiffusion2 在 `rfdiffusion2/`。Release 的实际下载地址由维护者发布后填写，不使用虚构下载链接。

所有分发文件的当前哈希保存在 SHA256SUMS；`provenance/source-files.json` 另记移植前的来源和哈希。历史执行记录中的路径已换成占位根目录，它们用于解释来源，不作为可直接执行的命令。
