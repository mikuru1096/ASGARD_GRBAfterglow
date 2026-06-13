# 私有网页文档发布

ASGARD 的网页文档由 `mkdocs.yml` 从 `doc/` 目录构建，并通过 GitHub Actions 发布到 GitHub Pages。目标远端是 `asgard-private`，站点可见性应在 GitHub Pages 设置中设为 `Private`，使访问权限跟随私有仓库的 collaborators / teams。

## 发布目标

- GitHub repository: `mikuru1096/ASGARD_private`
- Git remote: `asgard-private`
- Source: GitHub Actions
- Build command: `mkdocs build --strict`
- Artifact directory: `site`

## GitHub 设置

在 `mikuru1096/ASGARD_private` 中执行一次性设置：

1. 打开 `Settings -> Pages`。
2. 在 `Build and deployment` 中把 `Source` 设为 `GitHub Actions`。
3. 在 Pages visibility 中选择 `Private`。
4. 在 `Settings -> Collaborators and teams` 中只保留允许阅读文档的合作者或团队。

这些设置是 GitHub 仓库状态，不在源码文件中保存。若 Pages 设置页没有 `Private` 选项，说明当前仓库或组织没有可用的私有 Pages 访问控制能力，不能用 GitHub Pages 满足“仅合作者可见”。

## 本地构建检查

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run --with "mkdocs<2" --with "mkdocs-material>=9.5" mkdocs build --strict'
```

本地构建产物位于 `site/`，属于生成文件，不应提交。

## 发布流程

文档改动完成后，将相关提交推送到私有远端：

```bash
rtk git push asgard-private main
```

推送后，`Publish private documentation` workflow 会构建并部署 Pages。部署完成后，只有 `ASGARD_private` 仓库有权限的用户能打开站点。

## 维护原则

- 文档源码继续放在 `doc/`，不要把 `site/` 提交进仓库。
- 新文档需要加入 `mkdocs.yml` 的 `nav`，否则严格构建会失败或网页入口不完整。
- 修改代码结构、调用链或 public API 时，同步更新 `doc/code_overview.md`、`doc/source_tree.md`、`doc/call_chain.md`、`doc/public_api.md` 中受影响页面；这些页面已经在 `mkdocs.yml` 的网页导航中注册。
- 删除文档或运行时中转层时，同时搜索 `doc/`、`README.md` 和 `mkdocs.yml`，确保网页文档没有旧入口、旧模块名或旧调用链。
- Private visibility 只控制网页访问，不替代仓库权限审计；敏感内容仍应只放在私有仓库。
