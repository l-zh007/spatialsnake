# Spatialsnake：GitHub 发布与 Bioconda recipe 提交流程

本文用于把当前 Spatialsnake 源码发布到 GitHub，并把
`bioconda/recipes/spatialsnake/meta.yaml` 提交到
`bioconda/bioconda-recipes`。命令默认从项目根目录执行；先设置：

```bash
export PROJECT_ROOT="$(git rev-parse --show-toplevel)"
cd "$PROJECT_ROOT"
```

官方资料：

- [Bioconda 初始设置](https://bioconda.github.io/contributor/setup.html)
- [Bioconda 提交流程](https://bioconda.github.io/contributor/workflow.html)
- [Bioconda recipe 规范](https://bioconda.github.io/contributor/guidelines.html)
- [Bioconda 本地测试](https://bioconda.github.io/contributor/building-locally.html)

## 0. 当前状态与发布阻塞项

当前 recipe 是一个草案，不能直接提交。发布前必须处理以下事项：

- 最终 GitHub 仓库已确认为 `l-zh007/spatialsnake`。发布前仍需用
  `git remote -v` 确认本地 remote、README、`pyproject.toml` 和 recipe
  均指向该仓库。
- `spatialsnake/command_line.py`、`setup.py` 和 recipe 当前版本都是
  `0.0.1`。当前源码已明显晚于已发布的 `0.0.1`，应发布新版本，而不是覆盖旧版本。
- recipe 中的 `0.0.1` URL 和 SHA256 对应旧 sdist；该归档不包含当前新增的
  workflow/helper 文件。必须先发布新的不可变源码归档，再填写新 URL 和 SHA256。
- `requirements.txt` 把 `spatialdata`、`spatialdata-io`、
  `spatialdata-plot` 等声明为必需依赖，但 recipe 的 `requirements.run`
  没有完整声明这些依赖。Bioconda 包不能用“安装后再运行 pip/CRAN/GitHub
  下载”代替必需运行依赖；缺失依赖需要先进入 conda-forge/Bioconda，或把相应功能
  明确拆成可选功能。
- recipe 当前有大量精确版本 pin。只保留上游确实要求或为兼容性验证过的约束；
  其余尽量使用范围，避免与 Bioconda 全局 pin 冲突。

## 1. 先填写发布信息

执行命令前填写下表，后续所有占位符都以此为准。

| 字段 | 填写值 |
| --- | --- |
| GitHub 用户名 | `<GITHUB_USER>` |
| GitHub 仓库名 | `spatialsnake`（确认后填写） |
| 新版本号 | `<VERSION>`，例如 `0.0.2` |
| Git tag | `v<VERSION>` |
| 源码归档来源 | PyPI sdist / GitHub tag tarball（二选一） |
| 源码 URL | `<SOURCE_URL>` |
| SHA256 | `<SOURCE_SHA256>` |
| Bioconda fork URL | `https://github.com/<GITHUB_USER>/bioconda-recipes.git` |
| recipe 维护者 | `<BIOCONDA_MAINTAINER>` |

Shell 变量示例：

```bash
export GITHUB_USER="<GITHUB_USER>"
export GITHUB_REPO="spatialsnake"
export RELEASE_VERSION="<VERSION>"
export BIOCONDA_MAINTAINER="<BIOCONDA_MAINTAINER>"
```

## 2. GitHub 上传范围

应上传：

- `spatialsnake/`：Python、R、Snakemake 源码和 workflow 环境文件；
- `tests/`：项目测试；
- `resources/`：示例配置、logo 所需资源和 `resources/data/`；
- `spatialsnake-logo.png`、README、CHANGELOG、LICENSE；
- `pyproject.toml`、`setup.py`、`MANIFEST.in` 和依赖文件；
- `.github/`、`.snakemake-workflow-catalog.yml`；
- `bioconda/`：本项目保存的 recipe 草案和本文档。

不上传：

- `build/`、`dist/`、`*.egg-info/`；
- `__pycache__/`、`*.pyc`、pytest/Snakemake/Conda 缓存；
- 本机虚拟环境、日志、结果目录、IDE 配置和凭据文件；
- `docker/`：本项目按当前约定保存在仓库外；
- token、私钥、`.env`、个人绝对路径或未脱敏的临床/受试者信息。

当前 `resources/data/` 约 437 MB，最大单文件约 10.1 MB，没有触及 GitHub
100 MB 的单文件硬限制。但数据会显著增大 clone 体积；若以后出现大文件或频繁更新，
改用 Git LFS 或研究数据仓库，并在 README 中记录许可证、来源和校验值。

检查待提交内容：

```bash
git status --short
git ls-files --others --exclude-standard
find . -path ./.git -prune -o -type f -size +95M -print
```

敏感信息检查（输出必须为空；合法的 GitHub Actions secret 引用除外）：

```bash
rg -l --hidden -g '!.git/**' \
  -e 'BEGIN (RSA |EC |OPENSSH )?PRIVATE KEY' \
  -e 'AKIA[0-9A-Z]{16}' \
  -e 'github_pat_[A-Za-z0-9_]{20,}' \
  -e 'ghp_[A-Za-z0-9_]{20,}' \
  -e 'sk-[A-Za-z0-9]{20,}' \
  -e 'xoxb-[A-Za-z0-9-]{10,}' .

rg -n --hidden -g '!.git/**' '/users/[^/[:space:]]+|/home/[^/[:space:]]+' .
```

如果凭据曾经进入 Git 历史，仅从当前文件删除并不安全：先立即撤销/轮换凭据，再用
`git-filter-repo` 清理历史；不要在不理解影响时强制推送共享分支。

## 3. 统一项目名称和版本

先确认远端：

```bash
git remote -v
git remote get-url origin
```

若最终仓库是 `https://github.com/<GITHUB_USER>/spatialsnake`：

```bash
git remote set-url origin "git@github.com:${GITHUB_USER}/${GITHUB_REPO}.git"
```

同步以下文件中的版本和仓库 URL：

- `spatialsnake/command_line.py` 的 `__version__`；
- `setup.py` 的 `version`（建议后续删除重复版本源，只保留动态版本）；
- `bioconda/recipes/spatialsnake/meta.yaml` 的 version；
- `pyproject.toml`、README 和 recipe 的 Homepage/Issues/dev_url。

核对版本只能出现预期值：

```bash
rg -n '__version__|version[[:space:]]*=' \
  spatialsnake/command_line.py setup.py pyproject.toml \
  bioconda/recipes/spatialsnake/meta.yaml
```

## 4. 在干净环境验证源码

先运行不依赖外部数据的测试：

```bash
python -m compileall -q spatialsnake tests
python -m pytest -q
```

构建新的 sdist 和 wheel：

```bash
python -m pip install --upgrade build twine
python -m build
python -m twine check dist/*
```

确认 sdist 包含 CLI 所需文件：

```bash
tar -tzf "dist/spatialsnake-${RELEASE_VERSION}.tar.gz" | \
  rg 'spatialsnake/(command_line.py|config.yaml|workflow/Snakefile)'

tar -tzf "dist/spatialsnake-${RELEASE_VERSION}.tar.gz" | \
  rg 'workflow/envs/install_packages.yaml|install_packages_conda.R|install_packages_github_conda.R'
```

`resources/data/` 位于仓库中，但当前不在 Python 包目录内，也未由 `MANIFEST.in`
纳入 sdist。这适合“GitHub 示例数据”定位；若程序安装后必须访问它，应先明确许可证和
包体积策略，再将其放入包数据或外部数据仓库，不能假设 Bioconda 包会自动包含它。

## 5. 提交并发布 GitHub 版本

提交前人工检查 diff，尤其不要把移动到仓库外的文件重新加入：

```bash
git status --short
git diff --check
git diff --stat
git diff -- . ':!resources/data'
```

添加明确允许上传的内容：

```bash
git add .gitignore .github README.md CHANGELOG.md LICENSE \
  MANIFEST.in pyproject.toml setup.py requirements.txt requirements-extended.txt \
  spatialsnake tests resources spatialsnake-logo.png bioconda
git status --short
```

确认 staged 文件中没有生成物：

```bash
git diff --cached --name-only | \
  rg '(^|/)(build|dist|__pycache__|\.pytest_cache|\.snakemake|\.conda_pkgs|docker)(/|$)|\.py[co]$' \
  && echo '发现不应提交的文件，请先移出暂存区' \
  || echo 'staged 文件范围正常'
```

提交和推送：

```bash
git commit -m "Prepare Spatialsnake ${RELEASE_VERSION}"
git push -u origin main
git tag -a "v${RELEASE_VERSION}" -m "Spatialsnake ${RELEASE_VERSION}"
git push origin "v${RELEASE_VERSION}"
```

不要复用或移动已公开的 tag。若还要发布到 PyPI，先上传 TestPyPI 验证，再上传 PyPI；
Bioconda 只需要一个公开、稳定、带 SHA256 的源码归档。

## 6. 填写 `meta.yaml`

本项目只向 Bioconda PR 提交：

```text
recipes/spatialsnake/meta.yaml
```

不要把本仓库的 `bioconda/README.md`、本文档、测试数据、Docker 文件或第三方包 recipe
复制进 `bioconda-recipes`。

### 6.1 版本、URL 和 SHA256

PyPI sdist 方案：

```yaml
{% set name = "spatialsnake" %}
{% set version = "<VERSION>" %}

source:
  url: https://pypi.io/packages/source/s/spatialsnake/spatialsnake-{{ version }}.tar.gz
  sha256: <SOURCE_SHA256>
```

GitHub tag tarball 方案：

```yaml
{% set version = "<VERSION>" %}

source:
  url: https://github.com/<GITHUB_USER>/spatialsnake/archive/refs/tags/v{{ version }}.tar.gz
  sha256: <SOURCE_SHA256>
```

下载并计算哈希：

```bash
curl -L "<SOURCE_URL>" -o "/tmp/spatialsnake-${RELEASE_VERSION}.tar.gz"
sha256sum "/tmp/spatialsnake-${RELEASE_VERSION}.tar.gz"
tar -tzf "/tmp/spatialsnake-${RELEASE_VERSION}.tar.gz" | head
```

把结果原样填入 recipe。新版本 `build.number` 必须从 `0` 开始；只修改同一上游版本的
recipe 时才递增 build number。

### 6.2 构建和依赖

纯 Python 包可保留：

```yaml
build:
  noarch: python
  number: 0
  script: {{ PYTHON }} -m pip install . -vv --no-deps --no-build-isolation
```

逐项处理 `requirements.txt` 和实际 import：

| 上游依赖 | conda 包名 | host | run | 版本约束依据 | 已验证 |
| --- | --- | --- | --- | --- | --- |
| `<PYTHON_PACKAGE>` | `<CONDA_PACKAGE>` | 是/否 | 是/否 | `<UPSTREAM_OR_TEST>` | `[ ]` |
| `<R_OR_BIOC_PACKAGE>` | `<R-/BIOCONDUCTOR-NAME>` | 是/否 | 是/否 | `<UPSTREAM_OR_TEST>` | `[ ]` |

原则：

- 构建工具放 `host`，运行时直接或间接需要的包放 `run`；
- recipe 不能指定依赖来自某个特定 channel；依赖必须已在 defaults、conda-forge、
  Bioconda，或在同一 PR 中提交合规 recipe；
- 不要在 build、post-link 或 import 时联网执行 pip、CRAN、GitHub 安装；
- `spatialsnake install-packages` 可以是用户显式选择的扩展功能，但不能用于补齐
  `requirements.txt` 已声明的核心依赖；
- `pip`/`setuptools` 不应仅因为安装脚本使用过就自动放入 `run`。只有程序运行时
  确实调用它们才保留，并在 PR 中解释；
- 直接 CLI/import 测试涉及的包必须是 `run` 依赖；
- 精确 pin 必须有兼容性证据，优先写最小必要范围。

检查 conda 中是否存在依赖：

```bash
conda search -c conda-forge -c bioconda '<CONDA_PACKAGE>'
```

当前 recipe 特别需要重新判定这些被排除的核心包：

```text
spatialdata, spatialdata-io, spatialdata-plot, squidpy,
napari-spatialdata, cell2location, scvi-tools, torch
```

如果核心功能缺少任何一个包，先为缺失依赖选择“提交依赖 recipe、上游改为可选依赖、
或暂缓 Bioconda”之一；不要提交一个只有 `spatialsnake --version` 能通过的残缺包。

### 6.3 测试和元数据

至少保留：

```yaml
test:
  imports:
    - spatialsnake
  commands:
    - spatialsnake --version
    - spatialsnake --help
```

再增加一个无需网络和外部大数据、能覆盖核心 workflow 初始化的轻量 smoke test。
`test.commands` 只能依赖 `run` 依赖；需要额外测试依赖时使用同目录 `run_test.sh`，
并注意 mulled test 不会获得 `test.requires`。

核对 `about`：

- `home`、`dev_url`、`doc_url` 均可访问且指向正确项目；
- `license: MIT` 与仓库 LICENSE 一致，`license_file: LICENSE` 存在于 sdist；
- `summary` 准确；
- `extra.recipe-maintainers` 是有效 GitHub 用户名。

## 7. 提交到 Bioconda

Fork 并克隆：

```bash
git clone "https://github.com/${GITHUB_USER}/bioconda-recipes.git"
cd bioconda-recipes
git remote add upstream https://github.com/bioconda/bioconda-recipes.git
git checkout master
git pull upstream master
git push origin master
git checkout -b "add-spatialsnake-${RELEASE_VERSION}"
```

只复制 recipe：

```bash
mkdir -p recipes/spatialsnake
cp "${PROJECT_ROOT}/bioconda/recipes/spatialsnake/meta.yaml" \
  recipes/spatialsnake/meta.yaml
git status --short
```

建立本地测试环境并运行官方检查：

```bash
mamba create -n bioconda-test -c conda-forge -c bioconda bioconda-utils
conda activate bioconda-test
bioconda-utils lint --git-range master
bioconda-utils build --docker --mulled-test --git-range master
```

没有 Docker 时可先运行：

```bash
bioconda-utils build --git-range master
```

但最终仍要以 Bioconda PR CI 的容器化结果为准。

确认 PR 只包含一个文件：

```bash
git diff --name-only master...HEAD
git diff -- recipes/spatialsnake/meta.yaml
```

提交并推送：

```bash
git add recipes/spatialsnake/meta.yaml
git commit -m "Add spatialsnake ${RELEASE_VERSION}"
git push -u origin "add-spatialsnake-${RELEASE_VERSION}"
```

向 `bioconda/bioconda-recipes:master` 开 PR。CI 失败时只修改 recipe 或确有必要的
同目录测试/补丁文件，重新 lint/build 后推送。非 Bioconda 成员在 CI 通过后评论：

```text
@BiocondaBot please add label
```

## 8. 合并后的验证

等待 channel 同步后，在全新环境验证：

```bash
mamba create -n spatialsnake-release -c conda-forge -c bioconda \
  "spatialsnake=${RELEASE_VERSION}"
conda activate spatialsnake-release
spatialsnake --version
spatialsnake --help
python -c "import spatialsnake; print(spatialsnake.__file__)"
```

再运行 recipe 中定义的核心 smoke test。若核心功能仍要求执行
`spatialsnake install-packages` 才能启动，应将其视为 recipe 依赖设计未完成，而不是
正常安装步骤。

## 9. 最终勾选表

- [ ] GitHub 仓库名、remote 和项目 URL 已统一。
- [ ] 版本号在源码、构建元数据和 recipe 中一致。
- [ ] GitHub tag/发布归档不可变，SHA256 已从下载文件重新计算。
- [ ] sdist 包含 CLI 需要的 workflow、YAML、R 和 Snakemake 文件。
- [ ] `resources/data/` 已确认可公开、许可证清楚且没有受试者隐私信息。
- [ ] staged diff 不含 Docker、构建产物、缓存、字节码和凭据。
- [ ] 每个核心依赖都在 `requirements.run` 中可由 conda 解析。
- [ ] 构建、安装和测试过程不联网安装 pip/CRAN/GitHub 包。
- [ ] `bioconda-utils lint` 通过。
- [ ] Docker build 与 mulled test 通过。
- [ ] Bioconda PR 只包含 `recipes/spatialsnake/meta.yaml`（及确有必要的同目录文件）。
- [ ] 合并后在全新环境完成安装和核心 smoke test。
