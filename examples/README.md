# Debian 示例结果

本目录是练习的只读参考输入。`proteinmpnn/` 保存六条真实设计及概率；`haddock_stages/` 保存各阶段 capri_ss.tsv；`debian13_20261023/` 保存其他案例的结构、统计和过滤记录。使用 `python practice.py analyze 章号 --output runs/新名称` 复算，新文件不写回这里。

JSON 或 CSV 中的 `${REPOSITORY_ROOT}` 和 `${ORIGINAL_PROJECT_ROOT}` 表示经过移植整理的历史路径，仅用于来源说明，不能作为 shell 环境变量直接回放。当前运行命令由 practice.py 和 config/tools.json 生成。原始字节哈希见 provenance/source-files.json，当前字节哈希见 SHA256SUMS。
