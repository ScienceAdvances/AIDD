# ESMFold2 兼容性评估（仅评估，尚未安装或运行）

评估日期：2026-10-23。按用户要求，ESMFold2 环境及驱动改动必须另行确认。本轮没有安装其环境、下载其权重、修改驱动或软件源。

## 本机与官方要求

|项目|本机|教材固定版本要求/证据|判断|
|---|---|---|---|
|驱动|550.163.01，nvidia-smi 显示 CUDA 12.4|官方运行时使用 CUDA 13；NVIDIA 要求驱动至少 580|当前驱动不满足|
|显存|RTX 3060，12GB|官方 H200 测量最小 66-token 例峰值 13.20 GiB；常驻 12.862 GiB|原生完整模型在本机极可能显存不足|
|主存|约 15 GiB，另有 31 GiB swap|官方未提供本配置的 CPU offload 验证|不能把交换空间当显存；卸载路线尚无可复现保证|
|磁盘|尚未下载 ESMFold2 专用资源|官方预计约 33GB，另需安装临时空间|需预留空间|
|现有 BoltzGen|0.3.2，已完成五阶段|教材较新固定提交有额外 ESMFold2 worker|现有结果不包含 ESMFold2|

显存数字来自 H200 上的官方测量，不是 RTX 3060 上的实际测试，也不构成所有硬件的严格下界。它足以说明：仅升级驱动并不能合理承诺此机可跑完整模型。[官方性能测量](https://raw.githubusercontent.com/HannesStark/boltzgen/a07fea1153111ef453aaea7a15697ee05de9291f/docs/esmfold2-performance.md)

运行时固定 Python 3.12、ESM 3.4.1.post1、完整 full2021 checkpoint 和 ESMC-6B；模型资源与现有 BoltzGen 环境隔离。[官方安装说明](https://raw.githubusercontent.com/HannesStark/boltzgen/a07fea1153111ef453aaea7a15697ee05de9291f/docs/esmfold2.md)、[固定依赖](https://raw.githubusercontent.com/HannesStark/boltzgen/a07fea1153111ef453aaea7a15697ee05de9291f/src/boltzgen/resources/runtime/esmfold2.txt)

CUDA 12 的次版本兼容机制不能跨到 CUDA 13。容器仍依赖宿主机驱动；`cuda-compat` 的官方支持范围是数据中心卡、指定 NGC Server Ready RTX 和 Jetson，不能据此宣称 GeForce RTX 3060 获得支持。[NVIDIA 驱动要求](https://docs.nvidia.com/deploy/cuda-compatibility/minor-version-compatibility.html)、[前向兼容范围](https://docs.nvidia.com/deploy/cuda-compatibility/forward-compatibility.html)

只读检查当前配置的软件源缓存：nvidia-driver 候选仍为 550.163.01，已配置 backports 也是 550 系列；没有发现可直接安装的 580 候选。未执行 apt update 或修改源。

## 可选方案和待确认边界

1. **推荐：在更大显存、驱动 >=580 的 GPU 上补跑 ESMFold2。** 按官方测量推断建议至少 24GB，32/48GB 留有更多余量；具体大小仍需按目标裁剪长度验证，不保证 24GB 能跑任意输入。可复用本轮候选，补做源序列映射、ESMFold2 scoring、analysis 和 filtering。外部机器位置、数据传输与任何费用需用户指定/批准；本轮没有上传或付费。
2. **坚持本机：先批准独立的实验性移植环境，再评估 CUDA 12 / CPU 分段卸载。** 这不是官方已验证路线；目前没有足够证据承诺完整 checkpoint 能在 12GB GPU / 15GB RAM 下成功。需要逐项验证精度、资源占用与输出，失败也应保留记录。不得以 Fast checkpoint 或其他模型代替教材 full2021 并宣称同一流程。
3. **只升级驱动：当前不推荐单独做。** 它只能解决 CUDA 门槛，显存限制仍在。若用户仍希望升级，必须先给出驱动来源、内核模块兼容、回退与重启安排，再另行批准具体系统改动。

## 补跑时的输入与验收

- 使用教材固定的 BoltzGen 提交 a07fea1153111ef453aaea7a15697ee05de9291f；不覆盖已完成的 0.3.2 结果。
- 旧版本产物缺少完整 `source_context`。必须从原始输入重建每条链的完整实验构建体序列、保留位置及最终设计序列；不能用裁剪片段冒充完整来源序列。必要时在新目录重新生成。
- 固定 full2021 / ESMC-6B 权重、20 loops、200 steps、5 samples 和 seed 0；保存输入哈希、方向性 ipSAE、五样本评分、结构和 PAE。
- 新主线按 ESMFold2 最小方向性 ipSAE 排名；本轮 Boltz iPTM 不能填充该字段。现有两条候选的 Boltz 结构筛选均未通过，补跑 ESMFold2 也不会抹掉这个事实。

输入映射、评分及续跑规则依据[教材固定版本文档](https://raw.githubusercontent.com/HannesStark/boltzgen/a07fea1153111ef453aaea7a15697ee05de9291f/docs/esmfold2.md)。

## 本轮用户决定

用户已选择“先交付已完成结果，ESMFold2 保留待办”。以上方案为后续参考，不构成安装、升级、上传数据或付费的授权。本次交付不等待 ESMFold2。
