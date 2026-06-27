# 网页文档发布

ASGARD 的网页文档由 `mkdocs.yml` 从 `doc/` 目录构建，并通过 GitHub Actions 发布到 GitHub Pages。目标远端是 `asgard-private`，站点可见性应在 GitHub Pages 设置中设为 `Private`，使访问权限跟随私有仓库的 collaborators / teams。

当前另有一条 HEtools 托管路径：`https://hetools.cn/asgard-doc/`。该路径由 HEtools 的 frp 隧道路由到 `100.108.14.93` 上的本地静态文档服务，不覆盖 HEtools 首页的 Streamlit 应用。

公开 README 与网页文档需要互相指向：

- README URL: `https://github.com/mikuru1096/ASGARD_private#readme`
- Web docs URL: `https://hetools.cn/asgard-doc/`

当前网页导航已包含 `magnetized_rs_dg1d_tutorial.md`，作为磁化反向激波、`upstream_sigma` 物理闭合、DG1D 高阶输运、当前收敛阶和接口默认设置的专题教程；也包含 `prompt_internal_shock_tutorial.md`，作为 prompt internal-shock snapshot 的物理推导、代码入口、formal plotting 和边界说明。新增或改名该类专题页时，必须同步 `mkdocs.yml` 并跑 strict build。

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
- GPT_Academic port: `127.0.0.1:8502`
- static docs port: `127.0.0.1:8503`
- Latest docs backup: `/home/wangyun/Desktop/asgard_doc_webroot_backups/asgard-doc_20260619_152936_docs_publish`

`frpc.ini` contains multiple HTTP proxies on the same domains. The Streamlit proxy keeps `locations = /`; the ASGARD documentation proxy uses `locations = /asgard-doc`; the GPT_Academic proxy uses `locations = /gpt-academic`. frp uses the longer URL prefix first, while the HEtools app continues to serve `/`.

### GPT_Academic 挂载

GPT_Academic 通过 HEtools 域名公开入口挂载，但必须先登录才能进入 UI：

- Public URL: `https://hetools.cn/gpt-academic/`
- Source root on local Windows: `E:\chatpaper\gpt_academic`
- Remote app root: `/home/wangyun/Desktop/gpt_academic_web`
- Remote Python env: `/home/wangyun/Desktop/gpt_academic_web/.venv`
- Private config: `/home/wangyun/Desktop/gpt_academic_web/config_private.py`
- Login credential record: `/home/wangyun/Desktop/gpt_academic_login.txt`
- frpc route: `locations = /gpt-academic`, `local_port = 8502`
- Latest deployment backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/gpt_academic_deploy_20260626_111000`
- Latest account/workspace backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/gpt_account_workspace_20260626_155745`
- Latest home entry backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/home_gpt_academic_link_20260626_161458`
- Latest sidebar style backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/home_gpt_academic_sidebar_style_20260626_162129`
- Latest sidebar icon backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/home_gpt_academic_icon_20260626_162429`
- Latest online approval backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/gpt_account_admin_20260626_163137`
- Latest login-entry backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/gpt_login_entry_20260626_164519`
- Latest account form/password-rule backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/gpt_login_entry_20260626_170151`
- Latest secure download backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/gpt_file_tokens_20260626_174450`
- Latest admin credential config backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/gpt_set_jia_password_20260626_175644`
- Latest admin page builtin-account UI backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/gpt_admin_builtin_ui_20260626_175941`
- Latest old-style arXiv ID backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/gpt_old_arxiv_id_20260626_212711`
- Latest old LaTeX documentstyle backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/gpt_documentstyle_latex_20260626_215143`
- Latest legacy file-link redirect backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/gpt_file_redirect_20260626_215927`
- Latest arXiv cache browser backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/gpt_arxiv_cache_browser_20260626_221750`
- Latest arXiv cache MathJax backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/gpt_arxiv_cache_mathjax_20260626_222738`
- Latest LaTeX font-switch fix backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/gpt_latex_font_switch_fix_20260627_094843`
- Latest legacy AAS `aaspp4` compatibility backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/gpt_aaspp4_compat_20260627_104057`

部署包必须排除 `.venv/`、`__pycache__/`、`gpt_log/`、`history/`、`private_upload/`、临时目录和明文私密配置。线上配置放入 `config_private.py`，权限保持 `600`，包含 `WEB_PORT = 8502`、`CUSTOM_PATH = "/gpt-academic"`、`AUTO_OPEN_BROWSER = False`、`AUTHENTICATION`、模型名和 API key。不要把 API key 或登录密码写入源码、文档或聊天记录。

GPT_Academic 只允许绑定远端本机 `127.0.0.1:8502`，不能公开监听 `0.0.0.0`。当前部署版已把 `shared_utils/fastapi_server.py` 的 uvicorn 绑定地址固定为 `127.0.0.1`；重新同步上游代码时必须重新检查这一点。ASGARD 文档服务因此迁移到 `127.0.0.1:8503`，公网 `/asgard-doc/` 路径不变。

GPT_Academic 的远端 LaTeX 编译能力使用 rootless TeX Live 2026，不安装系统级 TeX Live：

- TinyTeX root: `/home/wangyun/.TinyTeX`
- TinyTeX bin: `/home/wangyun/.TinyTeX/bin/x86_64-linux`
- Available commands: `pdflatex`、`xelatex`、`bibtex`、`latexmk`、`tlmgr`、`kpsewhich`
- Installed scheme: `scheme-full`
- Installed package count after full install: `5110`
- Full install log: `/home/wangyun/Desktop/texlive_scheme_full_install.log`
- Runner patch backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/latex_runtime_20260626_135955`
- arXiv 0711.1980 compile fix backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/gpt_academic_latex_fix_20260626_141443`

GPT_Academic 源码会直接调用 `pdflatex`，按需调用 `xelatex`，并在参考文献场景调用 `bibtex`，所以只部署单文件 Tectonic 不够。`run_HEtools.py` 的 `gpt_academic_env()` 必须把 TinyTeX bin 加到 GPT_Academic 进程的 `PATH`；修改或重建 runner 后要检查正在运行的 GPT_Academic 进程环境是否仍包含该路径。`2604.23041` 暴露了 MNRAS/arXiv 模板对 `orcidlink`、`pdflscape`、`ctable`、`grffile`、`ae`、`lineno` 和 `babel-english` 等包的依赖；当前远端已扩展到 `scheme-full`，不再按单篇论文逐个补包。

`0711.1980` 的翻译编译问题暴露了两个 GPT_Academic LaTeX 处理边界：

- `crazy_functions/latex_fns/latex_actions.py` 的编译器选择必须把 `ctex` / `ctexart` / `ctexbook` / `ctexrep` 视为 `xelatex` 触发条件，否则 `pdflatex` 会在 CTeX Fandol 字体集处失败。
- `crazy_functions/latex_fns/latex_actions.py` 的 LaTeX 切分必须保护 `\be...\ee` / `\bea...\eea` 方程宏；`crazy_functions/latex_fns/latex_toolbox.py::fix_content` 必须把翻译片段中的裸 `&` 转义为 `\&`，并检查 `\be/\ee` 宏数量不被 LLM 输出改变。
- `crazy_functions/latex_fns/latex_toolbox.py::fix_content` 还必须修复旧式字体切换命令紧贴正文的问题，例如把 `{\it修正方程分析}` 转成 `{\it 修正方程分析}`；否则 XeLaTeX 会把 `\it修正...` 当作未定义控制序列。该问题由 `2604.23041` 暴露，远端已用 `merge_result.pkl` 重新合成 `merge_translate_zh.tex` 并生成 `translation/translate_zh.pdf`。
- `astro-ph/0112087` 暴露了另一类旧 arXiv 源码边界：`\documentstyle[12pt,aaspp4,graphicx]{article}` 里的 `aaspp4` 不在 TeX Live `scheme-full` 中。`crazy_functions/latex_fns/latex_toolbox.py` 的旧 `documentstyle` 转换不能生成 `\usepackage{aaspp4}`，而应插入最小 AAS 兼容层：`amsmath/amssymb`、`\affil`、`\keywords`、`\acknowledgments` 等。该文件还使用旧式 `${\mathbf \hat \phi}$` 数学写法，需要转成 `${\hat{\boldsymbol{\phi}}}$`；LLM 译文中出现的 `\cite{...}}）` 这类 citation/ref 后多余右括号，也由 `fix_content` 在标点前消解。远端已重新生成 `astro-ph/0112087/workfolder/merge.pdf`、`merge_translate_zh.pdf` 和 `translation/translate_zh.pdf`。

修复后，用 `gpt_log/arxiv_cache/0711.1980/workfolder/merge_result.pkl` 重新合成的 `codex_remerge_after_latex_fix.tex` 可由 `xelatex` 编译出 PDF。

GPT_Academic 的 arXiv LaTeX 翻译入口支持新旧两类 arXiv ID。`crazy_functions/Latex_Function.py` 会把 `0711.1980`、`0711.1980v2`、`astro-ph/9808007`、`astro-ph/9808007v1`、`https://arxiv.org/abs/...` 和 `https://arxiv.org/pdf/...` 规范化为 `/abs/` URL。旧格式 ID 的缓存目录仍保留 `astro-ph/9808007` 层级；下载源码 tar 文件名使用 `astro-ph_9808007.tar`，避免把 ID 中的 `/` 误当成本地文件名目录。远端 smoke 已验证 `astro-ph/9808007` 可下载源码并解包出 `.tex` 文件。

旧 arXiv 源码可能仍使用 LaTeX2.09 的 `\documentstyle`，例如 `astro-ph/9808007` 的 `SA_astro-ph2.tex`。`crazy_functions/latex_fns/latex_toolbox.py` 的主文件识别同时接受 `\documentclass` 和 `\documentstyle`；合并源码时会把 `\documentstyle[12pt,graphicx]{article}` 转成 `\documentclass[12pt]{article}` 加 `\usepackage{graphicx}`，再沿用现有 `ctex` 注入和编译链。远端 smoke 已验证 `astro-ph/9808007` 可定位主文件并生成包含 `documentclass`、`graphicx` 和 `ctex` 的 merged TeX。

GPT_Academic 提供已缓存 arXiv 翻译浏览页：`https://hetools.cn/gpt-academic/arxiv-cache`。该页只允许已登录用户访问，扫描共享 `gpt_log/arxiv_cache`，仅展示已有中文译文 PDF 的论文；PDF 优先使用 `translation/translate_zh.pdf`，否则使用 `workfolder/merge_translate_zh.pdf`。标题和摘要只从本地 `workfolder/merge.tex` 或 `extract/**/*.tex` 提取，不主动联网补全；TeX 数学片段保持原样并由 MathJax 渲染，避免把 `$\Omega_a$` 清理成 `$_a$` 这类残缺形式；下载链接统一使用 `/gpt-academic/dl/<token>`，不得输出 `/file=/home/...` 明文路径。

翻译分段调用 LLM API 的默认并发上限由 `DEFAULT_WORKER_NUM.txt` 控制，当前本地和远端均为 `100`。`config.py::CONCURRENT_COUNT` 是 Gradio 队列并发，当前也为 `100`，不需要为了翻译分段并发单独调整；具体模型仍受 `request_llms/bridge_all.py::model_info[model]["can_multi_thread"]` 约束，不支持并发的模型会按原逻辑降为单线程。

GPT_Academic 验收命令：

```bash
rtk bash -lc 'source ~/.wsl_env && curl --max-time 30 -sS https://hetools.cn/gpt-academic/ | grep -m1 "data-hetools-account-request-script"'
rtk bash -lc 'source ~/.wsl_env && curl --max-time 30 -sS https://hetools.cn/gpt-academic/account-request | grep -m1 "name=\"password\""'
rtk bash -lc 'source ~/.wsl_env && curl -L --max-time 30 -sS -o /tmp/gpt_bad_login.html -w "%{http_code}\n" -X POST -d "username=wrong&password=wrong" https://hetools.cn/gpt-academic/login'
rtk bash -lc 'source ~/.wsl_env && curl --max-time 30 -sS -o /tmp/gpt_admin.html -w "%{http_code}\n" https://hetools.cn/gpt-academic/admin/accounts'
rtk bash -lc 'source ~/.wsl_env && curl -L --max-time 30 -sS -o /tmp/gpt_plain_file.body -w "%{http_code}\n" "https://hetools.cn/gpt-academic/file=/home/wangyun/Desktop/gpt_academic_web/gpt_log/jia/downloadzone/2026-06-26-15-22-40-merge_translate_zh.pdf"'
rtk bash -lc 'source ~/.wsl_env && curl --max-time 30 -sS -o /tmp/gpt_arxiv_cache.html -w "%{http_code}\n" https://hetools.cn/gpt-academic/arxiv-cache'
```

LaTeX runtime 验收命令：

```bash
ssh wangyun@100.108.14.93 'source ~/.profile && pdflatex --version | head -n 1 && xelatex --version | head -n 1 && bibtex --version | head -n 1 && latexmk --version | head -n 1'
ssh wangyun@100.108.14.93 'pid=$(pgrep -f "/home/wangyun/Desktop/gpt_academic_web/.venv/bin/python main.py" | head -n 1); tr "\0" "\n" < /proc/$pid/environ | grep "^PATH=.*TinyTeX"'
```

匿名访问和错误账号都不得进入 GPT_Academic UI，也不得触发模型调用。`config.py`、`config_private.py`、`docker-compose.yml` 等路径在公网应返回登录页或不可访问响应，不能返回源文件内容。

GPT_Academic 账号申请入口挂在同一子路径下：

- Account request URL: `https://hetools.cn/gpt-academic/account-request`
- Admin approval URL: `https://hetools.cn/gpt-academic/admin/accounts`
- Pending requests: `/home/wangyun/Desktop/gpt_academic_web/gpt_log/admin/account_requests.jsonl`
- Approved accounts: `/home/wangyun/Desktop/gpt_academic_web/gpt_log/admin/approved_accounts.json`
- Account helper: `/home/wangyun/Desktop/gpt_academic_web/tools/manage_accounts.py`

账号申请页是公开表单，只写入 pending request，不创建账号，也不进入 Gradio UI 或触发模型调用。申请表字段为账号、密码、邮箱和用途说明；密码要求至少 8 位，并且必须同时包含字母和数字。在线审批页只允许已登录的 `jia` 访问；管理员批准后会使用申请人填写的密码写入 approved accounts、创建用户工作目录，并刷新当前进程的 `gradio_app.auth`，新账号不需要重启即可登录。在线审批页也支持拒绝申请和撤销已批准账号，审批页不回显用户申请密码。`jia` 属于 `config_private.py::AUTHENTICATION` 内置管理员账号，不属于普通 approved accounts，审批页会单独显示内置管理员账号且不提供撤销按钮。

命令行仍可在远端查看和批准账号：

```bash
ssh wangyun@100.108.14.93 'cd /home/wangyun/Desktop/gpt_academic_web && .venv/bin/python tools/manage_accounts.py requests'
ssh wangyun@100.108.14.93 'cd /home/wangyun/Desktop/gpt_academic_web && .venv/bin/python tools/manage_accounts.py approve USERNAME --password STRONG_PASSWORD'
```

命令行批准不会触发进程内 auth 刷新，仍需要重启 GPT_Academic 进程或重启 `/home/wangyun/Desktop/run_HEtools.py`。`account_requests.jsonl` 和 `approved_accounts.json` 权限必须保持 `600`；不要把申请密码或生成密码写入仓库文档。

GPT_Academic 工作区按登录用户隔离：`gpt_log/<user>`、`private_upload/<user>` 和下载目录只允许当前登录用户访问。为了保留文献翻译复用能力，`gpt_log/arxiv_cache` 仍是共享缓存目录；`autogen` 也保留为系统自动生成目录。普通用户不再跨账号访问 `default_user` 工作区。

GPT_Academic 敏感下载链接不再暴露 `/file=/home/...` 明文服务器路径。`shared_utils/fastapi_server.py` 会把 HTML、JSON 和纯文本响应里的敏感 `/gpt-academic/file=...` 链接改写为 `/gpt-academic/dl/<token>`；token 使用 `gpt_log/admin/file_token_secret` 做 HMAC 签名，密钥权限必须保持 `600`。如果旧页面或历史消息里仍残留敏感 `/file=...` 链接，登录且授权后会自动 302 跳转到 `/dl/<token>`，未登录或跨用户访问仍返回 401/403，不直接暴露文件。

HEtools 首页 `https://hetools.cn/` 的正文 welcome 区和 sidebar 都暴露 `GPT_Academic` 入口。sidebar 使用 `st.sidebar.page_link` 并放在模块列表最后，与其它 HEtools module 入口保持同一观感。首页只提供跳转链接，登录保护仍由 `/gpt-academic/` 的 Gradio authentication 执行；登录页内直接注入整行 `申请账号` 按钮，指向 `/gpt-academic/account-request`。`GET /gpt-academic/` 直接入口会清理旧登录态并显示登录页；`POST /gpt-academic/login` 仍由 Gradio 原生登录路由处理，HEtools 中间件只在成功响应上追加一次性短期标记，允许登录成功后的首次入口跳转进入 UI，避免登录循环。

### HEtools 托管的 ASGARD 代码

HEtools 的 afterglow 页面从 `base_func/afterglow_base_func.py` 调用 ASGARD。线上运行代码目录是：

- ASGARD code root: `/home/wangyun/Desktop/HEtools_web_beta_v03/base_func/ASGARD_GRBAfterglow-main`
- HEtools glue file: `/home/wangyun/Desktop/HEtools_web_beta_v03/base_func/afterglow_base_func.py`
- Current deployed ASGARD HEAD: `4b6678a`
- Latest deployment backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/asgard_code_update_20260615_190356`

更新线上 ASGARD 代码时，只替换 `ASGARD_GRBAfterglow-main` 和 `afterglow_base_func.py` 中 ASGARD 胶水入口。不要改 HEtools 页面主体、`menu.py`、`pages/Afterglow_modeling.py` 的 UI 逻辑，也不要改 `/asgard-doc/` 文档路由。新的胶水入口必须继续提供 `constants` 和 `ASGARD_fs_fluxdensity(frequency, **kwargs)`，返回 `(Tobs, Fluxes)`，其中 `Fluxes` 形状保持为 `(Num_frequency, Num_Tobs)`。

线上替换前必须在远端 staging 目录完成：

- `-Wline-truncation` 检查覆盖 `Constants`、`Dynamics_forward`、`electron_forward_fullhide_1d`、`electron_radiation`、`SED_interpolation`、`radiation_gamma_gamma_absorption`、`hadronic_forward_1d` 的源闭包。
- 用 `/home/wangyun/anaconda3/envs/mylab/bin/python build_extensions.py --force` 编译上述模块。
- 调用 `ASGARD_fs_fluxdensity` 做最小 live-shape 检查，要求 `Tobs.shape == (180,)`、`Fluxes.shape == (1, 180)` 且全部 finite。
- 通过 staging 检查后再停止服务、备份旧目录和旧胶水、替换正式目录、重启 `/home/wangyun/Desktop/run_HEtools.py`。

替换后必须同时验收：

- `https://hetools.cn/` 返回 Streamlit 首页。
- `https://hetools.cn/asgard-doc/` 返回 ASGARD 文档。
- `127.0.0.1:8501`、`127.0.0.1:8502`、`127.0.0.1:8503` 正常监听。
- live `afterglow_base_func.py` 可导入，`ASGARD_fs_fluxdensity` 最小 live-shape 检查通过。

### HEtools Afterglow_modeling 页面

`https://hetools.cn/Afterglow_modeling` 由 `pages/Afterglow_modeling.py` 提供 UI，由 `base_func/afterglow_base_func.py` 提供 ASGARD / VegasAfterglow / jetsimpy 计算入口。当前页面功能基线：

- ASGARD 和 VegasAfterglow 支持可配置时间网格、额外单频 Hz、ASGARD `Formal` / `Quick` 精度预设、lightcurve 和 SED 输出。
- 外部介质支持 `ISM` 和 `Wind`。默认保持 `ISM`；`Wind` 使用 ASGARD `WindMedium(a_star, density_floor_cm3, density_cap_cm3)` 和 VegasAfterglow `Wind(A_star, floor, cap)`，只开放当前 backend 明确支持的 \(k=2\) 星风。页面参数为 `log A_star`、`log wind floor`、可选 `log wind cap`。
- SSC 作为显式开关提供给 ASGARD 和 VegasAfterglow，默认关闭；开启后 ASGARD 返回 forward synchrotron + SSC，总 flux 形状不变，VegasAfterglow 通过 `Radiation(..., ssc=True)` 计算总辐射。jetsimpy 保持 legacy comparison，不接 SSC。
- SED 图默认用 `nuFnu` 显示能量输出，可切换回 `Fnu`；开启 SSC 时 ASGARD 和 VegasAfterglow SED 同时画 total、sync、SSC 三类曲线，CSV 同时导出 `Fnu` 和 `nuFnu` 分量。SED y 轴以当前已绘制曲线的最大正有限峰值为上界，固定显示 6 dex 动态范围。
- 模型按 ASGARD -> VegasAfterglow -> jetsimpy 串行运行；不并行抢占 CPU。
- `st.cache_data(max_entries=64)` 缓存 raw arrays，不缓存 Matplotlib figure；cache key 包含模型名、物理参数、介质参数、时间网格、频率、ASGARD 精度预设、`include_ssc` 和 `HETOOLS_DEPLOYED_ASGARD_HEAD.txt`。
- 页面输出 `Lightcurve`、`SED`、`Downloads`、`Diagnostics` tabs；CSV/PNG 下载由页面即时从 raw arrays 和 figure 生成。
- jetsimpy 保持 legacy lightcurve comparison，不扩展 SED。

ASGARD SSC 线上依赖 `src/Radiation/radiation_ssc_spectrum*.so`。若重新部署 ASGARD 代码且该扩展缺失，需要在远端 ASGARD root 运行 `-Wline-truncation` 源闭包检查后执行：

```bash
rtk bash -lc 'source ~/.wsl_env && ssh wangyun@100.108.14.93 "cd /home/wangyun/Desktop/HEtools_web_beta_v03/base_func/ASGARD_GRBAfterglow-main && PATH=/home/wangyun/anaconda3/envs/mylab/bin:\$PATH /home/wangyun/anaconda3/envs/mylab/bin/python build_extensions.py --force --module radiation_ssc_spectrum"'
```

最近一次页面升级备份：

- Main upgrade backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/afterglow_upgrade_20260615_195413`
- Compatibility fix backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/afterglow_upgrade_fix_20260615_195625`
- Page state fix backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/afterglow_page_state_fix_20260615_195830`
- Width API fix backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/afterglow_width_fix_20260615_195940`
- SSC switch backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/afterglow_ssc_20260615_200922`
- SED SSC component backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/afterglow_sed_ssc_visibility_20260615_201754`
- SED nuFnu mode backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/afterglow_sed_nufnu_20260615_202242`
- Wind medium backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/afterglow_wind_medium_20260615_202424`
- SED adaptive y-axis backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/afterglow_sed_adaptive_y_20260615_203017`
- Vegas SED components and 6 dex y-axis backup: `/home/wangyun/Desktop/HEtools_web_beta_v03_backups/afterglow_vegas_sed_y6_20260615_203850`

页面改动后必须验收：

- `py_compile menu.py pages/Afterglow_modeling.py base_func/afterglow_base_func.py` 通过。
- `ASGARD_fs_fluxdensity` 标量和数组频率输入都保持 `(180,)`、`(Nfreq, 180)` 形状且 finite。
- ASGARD / VegasAfterglow lightcurve 和 SED 最小检查必须通过。
- `include_ssc=True` 时，ASGARD / VegasAfterglow lightcurve 和 SED 必须输出 finite 且形状不变。
- `medium_type="Wind"` 时，ASGARD / VegasAfterglow lightcurve 和 SED 必须通过；legacy 9 参数 tuple 继续默认 `ISM`。
- SED `nuFnu` 模式下，ASGARD SSC 分量在 payload/CSV/图像中可见，且 total 等于 sync + SSC。
- SED 自适应 y 轴改动后，`streamlit.testing.v1.AppTest` 勾选 `Wind + SSC + nuFnu` 运行无 exception。
- VegasAfterglow SED 诊断显示单点 `sed_time_s` 调用与完整时间网格截面差异约 `10^-3`；页面保留原 Vegas total 求解路径，但在 `include_ssc=True` 时导出并绘制 `sync` / `SSC` 分量，且验证 total 等于 sync + SSC。
- 重复调用同一 `cached_lightcurve_payload` 的第二次请求明显命中缓存。
- `streamlit.testing.v1.AppTest` 初始加载和默认 `Run` 均无 exception，并出现 `Lightcurve`、`SED`、`Downloads`、`Diagnostics` tabs。
- 公网 `https://hetools.cn/`、`https://hetools.cn/Afterglow_modeling`、`https://hetools.cn/asgard-doc/` 均返回 200。

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
- `run_HEtools.py` 同时启动 `frpc`、`run.sh`、GPT_Academic 和 `python3 -m http.server 8503 --bind 127.0.0.1 --directory /home/wangyun/Desktop/asgard_doc_webroot`。
- `ss -ltnp` 能看到 `127.0.0.1:8501`、`127.0.0.1:8502` 和 `127.0.0.1:8503`。

公网验收：

```bash
rtk bash -lc 'source ~/.wsl_env && curl -L --max-time 30 -sS https://hetools.cn/asgard-doc/ | grep -m1 "ASGARD 文档" && curl -I --max-time 30 https://hetools.cn/asgard-doc/assets/stylesheets/main.484c7ddc.min.css'
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
