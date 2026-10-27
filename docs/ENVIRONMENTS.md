# 环境与工具配置

基础分析建议使用 Linux 上的 Python 3.11 或 3.12。`env/requirements-analysis.txt` 给出适合读者安装的依赖范围，`env/requirements-debian-analysis.lock.txt` 记录作者已使用的版本。前者方便开始，后者用于尽量贴近教材环境；版本差异要随结果保留。

离线复算需要 NumPy、pandas、Matplotlib、Biopython、openpyxl 和 RDKit。它不需要 PyTorch 或任何设计模型。完整 MD 轨迹重分析再安装 `env/requirements-md.txt`，并提供 GROMACS。这里的安装命令由读者主动执行，仓库运行入口不会创建或修改环境。

## 每种工具使用自己的环境

```bash
cp config/tools.example.json config/tools.json
```

`analysis_python` 用于共享脚本和 ProteinMPNN 的 GPU 任务。若要重新运行 02/06，应把它指向同时具有基础分析依赖与 CUDA PyTorch 的解释器；只做 `practice.py analyze` 时，使用当前基础环境即可。ProteinMPNN 固定源码已经随附，权重另行下载。

`haddock3`、`gmx`、`gnina` 和 `boltzgen` 是可执行文件，可写已在 PATH 中的命令名或完整路径。`pyrosetta_python` 和 `boltz_python` 则指向各自环境的 Python。PyRosetta 由读者按[官方获取方式](https://www.pyrosetta.org/downloads)安装，本仓库不分发其安装包。

AF3 需要填写 `af3_python`、`af3_repo` 和 `af3_params`；先按 [AlphaFold 3 官方仓库](https://github.com/google-deepmind/alphafold3)获取代码与许可范围内的参数，确认 CCD 数据已构建。本文教材基线为 `3.0.5.dev59+g3033a1ce4` / 提交 `3033a1ce47a9d107be7db9e644208b63e2d84a92`。

RFdiffusion2 需要填写 `rfd_python`、`rfd_repo`、`ligand_python`、`chai_python`、`chai_predict`。其中 `chai_predict` 指向该固定仓库的 `lib/chai/predict.py`，并在对应环境安装仓库配套依赖；原生适配器替换容器启动命令，不替代模型实现。基线提交为 `d365cbf4db3958814a9f8e4f6f94fa309dfebc2b`，Chai 复核关闭可选 ESM 嵌入。建议仓库和输出路径仅使用 ASCII，避免 Hydra 对路径解析产生歧义。

RFantibody 需要 `rfantibody_repo`、`rfantibody_python`、`rfantibody_bin`，后者是含 `rfdiffusion`、`qvextract`、`qvscorefile` 的 bin 目录。固定提交为 `3e5a5bb1cb98befde9fc049ae3ff58cabb491f67`；执行器使用项目可写目录运行后续官方脚本，规避便利 CLI 切换到只读安装目录的问题。

`extra_env` 可保存特定环境要求，例如 cuDNN 库搜索路径或模型缓存位置。不要把访问令牌写进准备上传的示例配置；个人 `config/tools.json` 已被忽略。

## 教材版本基线

HADDOCK3 为 2026.10.1，GROMACS 为 2026.4 CUDA 构建，GNINA 为 1.3.3，RDKit 为 2026.09.1，PyRosetta 为 2026.3+releasequarterly.5e498f1409。ProteinMPNN 源码提交为 `8907e6671bfbfc92303b5f79c4b5e6ce47cdef57`，权重为 v_48_020。

Boltz 包为 2.2.1、模型为 Boltz-2；BoltzGen 为 0.3.2，实际执行 design→inverse_folding→folding→analysis→filtering 五阶段。其他版本的命令或模型路线可能不同，首次执行前查看命令预览，并将实际版本写入自己的实验记录。ESMFold2 未执行，不在本仓库的运行入口中。
