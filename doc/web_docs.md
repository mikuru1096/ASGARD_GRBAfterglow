# 网页文档发布

ASGARD 的网页文档由 `mkdocs.yml` 从 `doc/` 目录构建，并通过 GitHub Actions 发布到 GitHub Pages。目标远端是 `asgard-private`，站点可见性应在 GitHub Pages 设置中设为 `Private`，使访问权限跟随私有仓库的 collaborators / teams。

当前另有一条 HEtools 托管路径：`https://hetools.cn/asgard-doc/`。该路径由 HEtools 的 frp 隧道路由到 `100.108.14.93` 上的本地静态文档服务，不覆盖 HEtools 首页的 Streamlit 应用。

公开 README 与网页文档需要互相指向：

- README URL: `https://github.com/mikuru1096/ASGARD_private#readme`
- Web docs URL: `https://hetools.cn/asgard-doc/`

## 发布目标

### GitHub Pages

- GitHub repository: `mikuru1096/ASGARD_private`
- Git remote: `asgard-private`
- Source: GitHub Actions
- Build command: `mkdocs build --strict`
- Artifact directory: `site`

### HEtools

- Public URL: `https://hetools.cn/asgard-doc/`
- Tailscale host: `100.108.14.93`
- SSH user: `wangyun`
- HEtools app root: `/home/wangyun/Desktop/HEtools_web_beta_v03`
- Static webroot: `/home/wangyun/Desktop/asgard_doc_webroot/asgard-doc`
- frpc config: `/home/wangyun/frp_0.58.0_linux_amd64/frpc.ini`
- runner: `/home/wangyun/Desktop/run_HEtools.py`
- Streamlit port: `127.0.0.1:8501`
- static docs port: `127.0.0.1:8502`

`frpc.ini` contains two HTTP proxies on the same domains. The Streamlit proxy keeps `locations = /`; the ASGARD documentation proxy uses `locations = /asgard-doc`. frp uses the longer URL prefix for `/asgard-doc/...`, while the HEtools app continues to serve `/`.

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

### GitHub Pages 发布

文档改动完成后，将相关提交推送到私有远端：

```bash
rtk git push asgard-private main
```

推送后，`Publish private documentation` workflow 会构建并部署 Pages。部署完成后，只有 `ASGARD_private` 仓库有权限的用户能打开站点。

### HEtools 发布

HEtools 发布不使用仓库内的 `site/` 目录。先在 `/tmp` 构建，再把产物打包成带 `asgard-doc/` 前缀的静态目录：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && rm -rf /tmp/asgard_mkdocs_site /tmp/asgard_doc_webroot_payload /tmp/asgard_doc_webroot_payload.tgz && uv run --with "mkdocs<2" --with "mkdocs-material>=9.5" mkdocs build --strict --site-dir /tmp/asgard_mkdocs_site && mkdir -p /tmp/asgard_doc_webroot_payload/asgard-doc && cp -a /tmp/asgard_mkdocs_site/. /tmp/asgard_doc_webroot_payload/asgard-doc/ && tar -C /tmp/asgard_doc_webroot_payload -czf /tmp/asgard_doc_webroot_payload.tgz asgard-doc'
```

上传并解包到远端后，必须满足：

- `/home/wangyun/Desktop/asgard_doc_webroot/asgard-doc/index.html` 存在。
- `/home/wangyun/frp_0.58.0_linux_amd64/frpc verify -c /home/wangyun/frp_0.58.0_linux_amd64/frpc.ini` 通过。
- `run_HEtools.py` 同时启动 `frpc`、`run.sh` 和 `python3 -m http.server 8502 --bind 127.0.0.1 --directory /home/wangyun/Desktop/asgard_doc_webroot`。
- `ss -ltnp` 能看到 `127.0.0.1:8501` 和 `127.0.0.1:8502`。

公网验收：

```bash
rtk bash -lc 'source ~/.wsl_env && curl -L --max-time 30 -sS https://hetools.cn/asgard-doc/ | grep -m1 "ASGARD Documentation" && curl -I --max-time 30 https://hetools.cn/asgard-doc/assets/stylesheets/main.484c7ddc.min.css'
```

HEtools 根路径也要同步检查，确保 Streamlit 首页没有被文档路由接管：

```bash
rtk bash -lc 'source ~/.wsl_env && curl -L --max-time 30 -sS https://hetools.cn/ | grep -m1 "Streamlit"'
```

## 维护原则

- 文档源码继续放在 `doc/`，不要把 `site/` 提交进仓库。
- 新文档需要加入 `mkdocs.yml` 的 `nav`，否则严格构建会失败或网页入口不完整。
- 修改代码结构、调用链或 public API 时，同步更新 `doc/code_overview.md`、`doc/source_tree.md`、`doc/call_chain.md`、`doc/public_api.md` 中受影响页面；这些页面已经在 `mkdocs.yml` 的网页导航中注册。
- 删除文档或运行时中转层时，同时搜索 `doc/`、`README.md` 和 `mkdocs.yml`，确保网页文档没有旧入口、旧模块名或旧调用链。
- 修改 README 或网页文档发布路径时，同步更新 `README.md`、`doc/index.md` 和本页的互链 URL。
- Private visibility 只控制网页访问，不替代仓库权限审计；敏感内容仍应只放在私有仓库。
- `https://hetools.cn/asgard-doc/` 跟随 HEtools 公开访问面；若文档必须只对合作者可见，应继续使用 GitHub Pages private visibility，或在 HEtools/frp/nginx 层加入明确认证。
