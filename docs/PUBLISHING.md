# 把练习仓库上传到 GitHub

上传本目录的内容即可。完整运行目录、模型权重、个人配置、下载缓存和新实验输出已经用 .gitignore 排除；单独的数据压缩包适合放 GitHub Releases。

先在 GitHub 建立一个空仓库，仓库名可使用 AIDD-practice。随后在本目录执行下列命令；将 YOUR_ACCOUNT 换成自己的账号：

```bash
python3 scripts/check_distribution.py
python3 practice.py verify
git init -b main
git add .
git status --short
git commit -m "Add chapter-based AIDD practice code and Debian examples"
git remote add origin https://github.com/YOUR_ACCOUNT/AIDD-practice.git
git push -u origin main
```

发布前在 `git status` 中确认没有个人 config/tools.json、runs/、模型权重或下载缓存。原创内容的版权说明保存在 COPYRIGHT.md，第三方许可证各自保留；如希望授予标准开源许可，由维护者另行选择并添加 LICENSE，本次整理没有代替作者授予新的开放许可。

可选数据包通过 Releases 的附件上传。基础仓库说明引用固定附件名与哈希，不预设账号；发布后可把真实 Release 地址补充到 docs/DATA.md。GitHub 普通 Git 会阻止超过 100 MiB 的单文件，浏览器上传单文件上限为 25 MiB，依据[官方大文件说明](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)。基础仓库按小文件分发；可选数据包不加入普通 Git 历史。

本地可上传目录与 ZIP 都已准备好。本次整理没有创建远程仓库、推送代码或发布 Release。
