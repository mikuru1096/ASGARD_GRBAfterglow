# 网页文档发布

MkDocs 源在 `doc/`，导航和主题在 `mkdocs.yml`。网页只发布当前跟踪文档，不把
本地构建物、benchmark 输出或外部仓库复制进站点。

## 本地检查

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tools/check_text_encoding.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run --with "mkdocs<2" --with "mkdocs-material>=9.5" mkdocs build --strict --site-dir /tmp/asgard_mkdocs_site'
```

使用 `/tmp` 输出，避免仓库内旧 `site/` 影响清理。strict build 必须无缺页、坏链
或导航孤儿。

## GitHub Pages

Pages 的正式来源是仓库配置的 workflow/branch。发布前确认：

1. `mkdocs.yml` 只引用存在页面；
2. public API、backend limits 和示例与当前代码一致；
3. 本地 strict build 通过；
4. 工作树不含生成的 `site/`；
5. 目标 branch 和 remote 正确。

只推送源码和配置，由既有发布流程生成站点。

## HEtools

HEtools 部署是独立运行环境，不是本文档源码的第二权威副本。更新时从已提交的同一
revision 构建，保留既有 `/gpt-academic` 与 ASGARD 路由边界，不在文档中记录主机
凭据、临时目录或会过期的手工复制命令。

若页面由 HEtools 反向代理，先验证静态资源 base URL、路由前缀和刷新后直达链接。
部署故障应在部署配置解决，不向 Markdown 加环境特判。

## 维护规则

- 页面只描述当前公开能力；计划进入 `TODO.md`，确认缺陷进入 `BUG.md`；
- 源码树和调用链只在 `code_overview.md` 维护；
- 算法公式只在算法总纲及专题页维护；
- 删除页面时同步导航和所有相对链接；
- 不提交生成站点、访问日志、截图缓存或部署密钥。
