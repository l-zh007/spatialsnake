# test3 Conda 环境 Docker 打包与分发手册

本文档说明如何把 `test3` Conda 环境、当前 Spatialsnake 源码及其 Python/R
双语言依赖封装为 Docker 镜像，并完成本地验证、镜像仓库上传和终端用户运行。

## 1. 已确认的运行基线

本次检查日期为 2026-06-12，源环境与宿主机信息如下：

| 项目 | 已确认值 |
| --- | --- |
| Conda 环境 | `/users/zhenghlin/miniconda3/envs/test3` |
| 源环境大小 | 约 13 GB |
| 操作系统 | Ubuntu 22.04.5 LTS |
| CPU 架构 | `linux/amd64`，宿主机为 `x86_64` |
| glibc | 2.35 |
| Python | 3.12.11 |
| R | 4.4.0 |
| Snakemake | 9.8.1 |
| PyTorch | 2.8.0+cu128 |
| CUDA 宿主机测试 | `torch.cuda.is_available() == True`，两张 RTX 4090 D |
| Docker Engine | 28.4.0 |
| Docker Buildx | 0.28.0 |
| 基础镜像 | `ubuntu:22.04`，固定到 SHA-256 摘要 |
| 环境归档 | 6.2 GB，SHA-256 `0e6a48279f77e2fa5b67c5326155f939496b025f240fc08c2155ce0d59c8a1b6` |
| 已验证本地镜像 | `spatialsnake-test3:2026-06-12`，12.7 GB |
| 已验证镜像 ID | `sha256:ffd71988f99a8fac2dd06828c2520323efae4a0e480a36458c6983efd97cbdbf` |

镜像当前只面向 `linux/amd64`。`conda-pack` 要求打包端与目标端操作系统类型
一致，因此不能把该 Linux 环境直接用于 Windows 容器或 macOS 原生容器。

## 2. Spatialsnake 架构和依赖链

Spatialsnake 是一个以 Snakemake 为调度层的空间转录组分析流程。主要层次如下：

1. CLI 层：`spatialsnake/command_line.py`
2. 工作流层：`spatialsnake/workflow/Snakefile` 与 `workflow/rules/*.smk`
3. Python 分析层：SpatialData、Scanpy、Squidpy、scvi-tools、cell2location、
   CellPhoneDB、CellCharter、BANKSY、LIANA、pySCENIC、PyDESeq2 等
4. R 分析层：Seurat、spacexr/RCTD、CellChat、clusterProfiler、
   ComplexHeatmap、AnnotationDbi、org.Hs.eg.db、org.Mm.eg.db 等
5. 数据层：以 SpatialData/Zarr 为统一对象，兼容 Visium、Visium HD、
   Visium Segment、Xenium、MERFISH/MERSCOPE 和 Stereo-seq

核心流程包括：

- 数据导入与统一格式转换
- 质量控制、过滤、归一化和批次校正
- 聚类、重聚类、手动注释和再注释
- cell2location 与 RCTD 注释
- 差异表达、富集分析
- CellPhoneDB、LIANA、CellChat、CellCharter、BANKSY、pySCENIC
- SpatialData 的拆分、合并与格式转换

`test3` 同时包含 Conda、pip、CRAN/Bioconductor 和 GitHub R 包。仅使用
`conda export` 无法完整保存 GitHub R 包，因此本方案以 `conda-pack` 的环境
归档作为镜像输入，并额外保存四类审计材料：

- `environment-linux-64.yml`：Conda 与外部包的完整 YAML 表示
- `environment-from-history.yml`：用户显式安装的 Conda 需求
- `explicit-linux-64.txt`：同平台 Conda 包 URL 锁定
- `pip-freeze.txt` 与 `r-packages.tsv`：pip 和 R 的精确版本清单

## 3. 容器架构设计

镜像构建采用以下顺序：

```text
Ubuntu 22.04 固定摘要
        |
        +-- ca-certificates + tini
        |
        +-- /opt/conda/envs/test3
        |      从 conda-pack 归档展开
        |      执行 conda-unpack 修复前缀
        |
        +-- 当前 Spatialsnake 源码
        |      pip install --no-deps --no-build-isolation .
        |
        +-- Python/R/CLI 构建期冒烟测试
        |
        +-- 非 root 用户 spatialsnake
               工作目录 /work
```

关键设计决定：

- 固定 Ubuntu 基础镜像摘要，避免同一 tag 随时间漂移。
- 归档 `test3` 的实际内容，确保非 Conda R 包也进入镜像。
- 使用 BuildKit 只读构建挂载展开归档，压缩包不会作为额外的 6.2 GB 镜像层保留。
- 构建时重新安装当前仓库源码，避免环境中的旧源码覆盖当前修改。
- 默认使用非 root 用户，容器入口为 `spatialsnake`。
- 原始数据和结果不写进镜像，使用 bind mount 挂载到 `/work`。
- 不在镜像中保存 Docker Hub token、GitHub token、数据库口令或私有数据。
- 构建上下文由 `Dockerfile.dockerignore` 严格限制，不会发送 `.git`、
  `.conda_pkgs`、日志和构建缓存。

## 4. 目录内容

```text
docker/test3/
├── Dockerfile
├── Dockerfile.dockerignore
├── README.zh-CN.md
├── export-test3.sh
├── conda_pack_compat.py
├── build-image.sh
├── verify-image.sh
├── push-image.sh
├── smoke_test.py
├── artifacts/
│   ├── .gitignore
│   └── test3-linux-64.tar.gz       # 生成物，不进入 Git
└── locks/
    ├── environment-linux-64.yml
    ├── environment-from-history.yml
    ├── explicit-linux-64.txt
    ├── pip-freeze.txt
    ├── r-packages.tsv
    ├── system-info.txt
    └── SHA256SUMS
```

## 5. 构建端前置检查

进入项目根目录：

```bash
cd /users/zhenghlin/beifen2/spatialdata
```

确认环境、平台和 Docker：

```bash
conda env list
conda run -n test3 python --version
conda run -n test3 Rscript -e 'cat(R.version.string, "\n")'
conda run -n test3 snakemake --version
docker version
docker buildx version
uname -m
getconf GNU_LIBC_VERSION
df -h .
docker system df
```

建议至少预留：

- 环境归档空间：8 到 15 GB
- Docker 构建临时空间：20 到 40 GB
- 最终镜像空间：本次实测 12.7 GB，仍建议至少预留 15 GB

不要在空间不足时运行 `docker system prune -a`。该命令会删除其他项目未使用的
镜像和构建缓存，应由管理员或资源所有者单独评估。

## 6. 源环境预处理和审计

### 6.1 验证 Python 和 R

```bash
/users/zhenghlin/miniconda3/envs/test3/bin/python \
  docker/test3/smoke_test.py
```

可选的 pip 元数据检查：

```bash
/users/zhenghlin/miniconda3/envs/test3/bin/python -m pip check
```

当前 `test3` 的 `pip check` 会报告：

- `pybanksy 1.3.4` 声明依赖 `python-igraph`，环境安装的分发名是 `igraph`
- `pybanksy 1.3.4` 声明 `numpy<2.0`，环境实际为 `numpy==2.2.6`

本次实际导入 `banksy` 成功，但这不等于所有 BANKSY 数据路径均已完成生物学
结果验证。不要仅为了消除元数据警告就把 NumPy 降级，因为 SpatialData、
Scanpy、Numba 和其他已固定组件可能受到影响。应在独立环境中完成兼容性回归后
再调整版本。

### 6.2 安装 conda-pack

Conda-pack 官方建议安装到 root/base 环境，而不是修改待打包环境：

```bash
conda install -n base -c conda-forge conda-pack
conda-pack --version
```

若企业环境禁止修改 base，可创建独立工具环境：

```bash
conda create -n conda-pack-tool -c conda-forge \
  python=3.12 conda-pack setuptools=80 -y
conda activate conda-pack-tool
```

### 6.3 导出清单并生成归档

默认按环境名查找 `test3`：

```bash
./docker/test3/export-test3.sh
```

若有多个 Conda 安装，建议明确指定：

```bash
CONDA_EXE=/users/zhenghlin/miniconda3/bin/conda \
CONDA_ENV_PREFIX=/users/zhenghlin/miniconda3/envs/test3 \
COMPRESS_LEVEL=1 \
./docker/test3/export-test3.sh
```

参数说明：

- `CONDA_EXE`：负责读取环境记录的 Conda 可执行文件。
- `CONDA_ENV_PREFIX`：待打包环境绝对路径，优先级高于环境名。
- `CONDA_ENV_NAME`：未给出 prefix 时使用，默认 `test3`。
- `COMPRESS_LEVEL`：gzip 压缩级别，默认 `1`。提高级别会减小归档，但耗时增加。
- `CONDA_PACK_PYTHON`：特殊情况下指定安装了 `conda-pack` 的 Python。

脚本完成后检查：

```bash
ls -lh docker/test3/artifacts/test3-linux-64.tar.gz
cd docker/test3
sha256sum --check --ignore-missing locks/SHA256SUMS
cd ../..
```

`test3` 存在一个不完整的
`conda-meta/spatialsnake-0.0.1-py_0.json` 记录，并且 pip 覆盖了
`scipy`/`numba` 的少量 `direct_url.json` 元数据。兼容包装器仅在读取阶段把
畸形记录视为普通文件，并启用 conda-pack 的 `ignore_missing_files`；它不会
改写源环境。Python/R 冒烟测试仍是归档可接受的前置条件。

## 7. Dockerfile 说明

基础镜像固定为：

```dockerfile
ARG BASE_IMAGE=ubuntu:22.04@sha256:4f838adc7181d9039ac795a7d0aba05a9bd9ecd480d294483169c5def983b64d
```

如需更新基础镜像，先拉取并记录新摘要：

```bash
docker pull ubuntu:22.04
docker image inspect ubuntu:22.04 --format '{{json .RepoDigests}}'
```

更新摘要后必须重新执行完整构建和验证。不要只把 tag 改成 `latest`。

环境安装位置固定为：

```text
/opt/conda/envs/test3
```

Dockerfile 展开归档后执行：

```bash
/opt/conda/envs/test3/bin/conda-unpack
```

该步骤负责修复从源 prefix 到容器 prefix 的文本和二进制前缀。执行后不要再移动
该目录。归档通过 Dockerfile `RUN --mount=type=bind` 只在构建步骤中读取，不会
写入最终镜像层；这是镜像从 19.3 GB 降至 12.7 GB 的关键。

## 8. 本地镜像构建

默认构建：

```bash
./docker/test3/build-image.sh
```

默认镜像名和 tag：

```text
spatialsnake-test3:2026-06-12
```

生产发布建议显式命名：

```bash
IMAGE_NAME=spatialsnake-test3 \
IMAGE_TAG=0.0.1-test3-20260612 \
APP_UID="$(id -u)" \
APP_GID="$(id -g)" \
./docker/test3/build-image.sh
```

参数说明：

- `IMAGE_NAME`：本地镜像仓库名。
- `IMAGE_TAG`：不可变版本标识，建议包含应用版本、环境名和日期。
- `APP_UID`/`APP_GID`：镜像默认非 root 用户的 UID/GID。

等价的手动命令：

```bash
docker build \
  --pull \
  --file docker/test3/Dockerfile \
  --tag spatialsnake-test3:0.0.1-test3-20260612 \
  --build-arg IMAGE_VERSION=0.0.1-test3-20260612 \
  --build-arg BUILD_DATE="$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
  --build-arg VCS_REF="$(git rev-parse --short=12 HEAD)" \
  --build-arg APP_UID="$(id -u)" \
  --build-arg APP_GID="$(id -g)" \
  .
```

构建完成后：

```bash
docker image ls spatialsnake-test3
docker image inspect spatialsnake-test3:0.0.1-test3-20260612
docker image history spatialsnake-test3:0.0.1-test3-20260612
```

## 9. 镜像验证

执行完整导入和 CLI 冒烟测试：

```bash
IMAGE_NAME=spatialsnake-test3 \
IMAGE_TAG=0.0.1-test3-20260612 \
./docker/test3/verify-image.sh
```

该脚本还会以宿主机 UID/GID 挂载临时工作目录，验证终端用户能够在 `/work`
写出结果文件。当前已完成的 CPU 验证包括：

- 构建期及运行期 Python/R 全依赖导入
- `spatialsnake --version` 与默认 `--help`
- 非 root 用户运行
- bind mount 写入和宿主文件权限
- 镜像层检查，确认环境归档未重复留在镜像中

单独检查版本：

```bash
docker run --rm \
  spatialsnake-test3:0.0.1-test3-20260612 \
  --version
```

检查 RCTD 所需 R 包：

```bash
docker run --rm \
  --entrypoint Rscript \
  spatialsnake-test3:0.0.1-test3-20260612 \
  -e 'library(Seurat); library(spacexr); library(schard); sessionInfo()'
```

检查 CellChat：

```bash
docker run --rm \
  --entrypoint Rscript \
  spatialsnake-test3:0.0.1-test3-20260612 \
  -e 'library(CellChat); cat(as.character(packageVersion("CellChat")), "\n")'
```

检查镜像内清单：

```bash
docker run --rm \
  --entrypoint bash \
  spatialsnake-test3:0.0.1-test3-20260612 \
  -lc 'ls -lh /opt/spatialsnake-manifest && cat /opt/spatialsnake-manifest/system-info.txt'
```

真实数据验证应至少选择一个小型样本，依次测试：

1. `integrate`
2. `preprocess`
3. `clustering`
4. 一条 R 路径，例如 `annotation --anno_algorithm RCTD`
5. 一条扩展 Python 路径，例如 `advance_analysis --runpipe liana`

冒烟测试验证软件可导入，不替代真实数据、统计结果和图形结果验收。

## 10. GPU 容器前置条件

源环境的 PyTorch 已在宿主机直接验证可识别两张 GPU，但当前 Docker daemon
执行 `docker run --gpus all ...` 返回：

```text
could not select device driver "" with capabilities: [[gpu]]
```

这表示宿主机 Docker 尚未配置 NVIDIA Container Toolkit，不是镜像内缺少
PyTorch CUDA 包。该环节需要管理员权限，当前任务未修改 Docker daemon。

管理员应按 NVIDIA 官方文档安装 toolkit，然后执行：

```bash
sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker
```

先验证 Docker GPU：

```bash
docker run --rm --runtime=nvidia --gpus all ubuntu:22.04 nvidia-smi
```

再验证本镜像：

```bash
VERIFY_GPU=1 \
IMAGE_NAME=spatialsnake-test3 \
IMAGE_TAG=0.0.1-test3-20260612 \
./docker/test3/verify-image.sh
```

只分配第 0 张 GPU：

```bash
docker run --rm \
  --gpus '"device=0"' \
  --entrypoint python \
  spatialsnake-test3:0.0.1-test3-20260612 \
  -c 'import torch; print(torch.cuda.get_device_name(0))'
```

## 11. 运行 Spatialsnake

### 11.1 准备宿主机工作目录

```bash
mkdir -p "$HOME/spatialsnake-run/data" "$HOME/spatialsnake-run/results"
cd "$HOME/spatialsnake-run"
touch sample.txt
```

推荐结构：

```text
spatialsnake-run/
├── data/
├── results/
├── sample.txt
└── preprocess.yaml
```

### 11.2 CPU 示例

挂载整个工作目录：

```bash
docker run --rm \
  --mount type=bind,src="$PWD",dst=/work \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --cpus 16 \
  --memory 64g \
  spatialsnake-test3:0.0.1-test3-20260612 \
  single_analysis sample.txt visium \
  --option=preprocess \
  --configfile=preprocess.yaml \
  --jobs=16
```

参数说明：

- `--rm`：任务结束后删除容器，不删除镜像和宿主机结果。
- `--mount`：把当前目录映射到容器 `/work`。
- `--user`：让输出文件归当前宿主机用户所有。
- `--env HOME=/tmp`：自定义 UID 运行时提供可写 HOME。
- `--cpus`/`--memory`：限制容器资源，避免占满宿主机。
- `--jobs`：Snakemake 最大并行作业数，应不大于可用 CPU。

如果希望原始数据只读，可分别挂载：

```bash
docker run --rm \
  --mount type=bind,src="$PWD/data",dst=/work/data,readonly \
  --mount type=bind,src="$PWD/results",dst=/work/results \
  --mount type=bind,src="$PWD/sample.txt",dst=/work/sample.txt,readonly \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  spatialsnake-test3:0.0.1-test3-20260612 \
  single_analysis sample.txt visium --option=integrate --jobs=8
```

Docker bind mount 默认可写。对原始数据、参考数据库和配置文件使用 `readonly`
可以降低误修改风险。

### 11.3 GPU 示例

```bash
docker run --rm \
  --gpus all \
  --shm-size=16g \
  --mount type=bind,src="$PWD",dst=/work \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  spatialsnake-test3:0.0.1-test3-20260612 \
  single_analysis sample.txt visium \
  --option=annotation \
  --anno_algorithm=cell2Location \
  --device=cuda \
  --jobs=8
```

`--shm-size` 对 PyTorch 多进程数据加载有帮助，应根据宿主机内存调整。

### 11.4 进入调试 Shell

```bash
docker run --rm -it \
  --entrypoint bash \
  --mount type=bind,src="$PWD",dst=/work \
  spatialsnake-test3:0.0.1-test3-20260612
```

容器内无需执行 `conda activate`，`PATH` 已指向
`/opt/conda/envs/test3/bin`。

## 12. 上传到 Docker Hub

### 12.1 创建仓库

在 Docker Hub 中创建仓库，例如：

```text
YOUR_DOCKERHUB_USER/spatialsnake-test3
```

如镜像包含不能公开分发的依赖、数据或许可证受限组件，应创建私有仓库并先完成
许可证审查。本镜像不应包含实验数据，但环境软件本身仍需逐项遵守上游许可证。

### 12.2 安全登录

交互式登录：

```bash
docker login
```

自动化场景使用访问令牌和标准输入，不要把 token 写在命令行参数或 Dockerfile：

```bash
printf '%s' "$DOCKERHUB_TOKEN" \
  | docker login --username YOUR_DOCKERHUB_USER --password-stdin
```

### 12.3 标记和推送

使用脚本：

```bash
REGISTRY=docker.io \
REGISTRY_NAMESPACE=YOUR_DOCKERHUB_USER \
REGISTRY_REPOSITORY=spatialsnake-test3 \
IMAGE_NAME=spatialsnake-test3 \
IMAGE_TAG=0.0.1-test3-20260612 \
./docker/test3/push-image.sh
```

等价手动命令：

```bash
docker image tag \
  spatialsnake-test3:0.0.1-test3-20260612 \
  YOUR_DOCKERHUB_USER/spatialsnake-test3:0.0.1-test3-20260612

docker image push \
  YOUR_DOCKERHUB_USER/spatialsnake-test3:0.0.1-test3-20260612
```

验证远端摘要：

```bash
docker buildx imagetools inspect \
  YOUR_DOCKERHUB_USER/spatialsnake-test3:0.0.1-test3-20260612
```

只有在版本 tag 验证完成后，才建议追加 `latest`：

```bash
docker image tag \
  YOUR_DOCKERHUB_USER/spatialsnake-test3:0.0.1-test3-20260612 \
  YOUR_DOCKERHUB_USER/spatialsnake-test3:latest

docker image push YOUR_DOCKERHUB_USER/spatialsnake-test3:latest
```

生产和论文复现应使用不可变版本 tag 或镜像 digest，不应只依赖 `latest`。

## 13. 上传到私有仓库

假设仓库地址为 `registry.example.org`：

```bash
docker login registry.example.org

REGISTRY=registry.example.org \
REGISTRY_NAMESPACE=bioinformatics \
REGISTRY_REPOSITORY=spatialsnake-test3 \
IMAGE_NAME=spatialsnake-test3 \
IMAGE_TAG=0.0.1-test3-20260612 \
./docker/test3/push-image.sh
```

目标镜像为：

```text
registry.example.org/bioinformatics/spatialsnake-test3:0.0.1-test3-20260612
```

私有仓库若使用自签名证书，应把受信任 CA 正确配置到 Docker daemon。不要通过
关闭 TLS 校验来绕过证书问题。

## 14. 终端用户拉取与使用

拉取指定版本：

```bash
docker pull \
  YOUR_DOCKERHUB_USER/spatialsnake-test3:0.0.1-test3-20260612
```

验证：

```bash
docker run --rm \
  YOUR_DOCKERHUB_USER/spatialsnake-test3:0.0.1-test3-20260612 \
  --version
```

按 digest 拉取可获得最强的内容固定：

```bash
docker pull \
  YOUR_DOCKERHUB_USER/spatialsnake-test3@sha256:<REMOTE_DIGEST>
```

运行：

```bash
cd "$HOME/spatialsnake-run"

docker run --rm \
  --mount type=bind,src="$PWD",dst=/work \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  YOUR_DOCKERHUB_USER/spatialsnake-test3:0.0.1-test3-20260612 \
  single_analysis sample.txt visium \
  --option=clustering \
  --jobs=16
```

用户不需要安装 Conda、Python 或 R，只需要兼容的 Docker Engine。GPU 用户还
需要宿主机 NVIDIA 驱动和 NVIDIA Container Toolkit。

## 15. 离线分发备选方案

如果目标机器不能访问镜像仓库：

```bash
docker image save \
  spatialsnake-test3:0.0.1-test3-20260612 \
  | gzip -1 > spatialsnake-test3_0.0.1-test3-20260612.tar.gz

sha256sum \
  spatialsnake-test3_0.0.1-test3-20260612.tar.gz \
  > spatialsnake-test3_0.0.1-test3-20260612.tar.gz.sha256
```

目标机器：

```bash
sha256sum --check \
  spatialsnake-test3_0.0.1-test3-20260612.tar.gz.sha256

gzip -dc \
  spatialsnake-test3_0.0.1-test3-20260612.tar.gz \
  | docker image load
```

镜像 tar 通常很大，传输前应确认文件系统、介质和网关的单文件大小限制。

## 16. 发布和维护策略

每次环境或源码更新都应：

1. 在 `test3` 中完成 Python/R 和真实数据回归。
2. 重新运行 `export-test3.sh`。
3. 记录新的 `locks/SHA256SUMS`。
4. 使用新版本 tag 构建，不覆盖旧版本 tag。
5. 运行 CPU 冒烟测试。
6. 在已配置 GPU runtime 的主机运行 GPU 测试。
7. 推送版本 tag。
8. 记录远端镜像 digest。
9. 最后再决定是否更新 `latest`。

建议版本格式：

```text
<spatialsnake版本>-<环境名>-<YYYYMMDD>
```

例如：

```text
0.0.1-test3-20260612
```

## 17. 常见问题

### 17.1 `test3-linux-64.tar.gz: not found`

先运行：

```bash
./docker/test3/export-test3.sh
```

并确认归档位于：

```text
docker/test3/artifacts/test3-linux-64.tar.gz
```

### 17.2 `requested access to the resource is denied`

检查：

```bash
docker login
docker image ls
```

目标 tag 的 namespace 必须是有推送权限的 Docker Hub 用户或组织。

### 17.3 容器输出文件权限不正确

构建时设置：

```bash
APP_UID="$(id -u)" APP_GID="$(id -g)" ./docker/test3/build-image.sh
```

或运行时设置：

```bash
--user "$(id -u):$(id -g)" --env HOME=/tmp
```

镜像默认用户 UID/GID 为构建时的 `APP_UID`/`APP_GID`。若终端用户与该数值不同，
直接挂载权限较严格的宿主目录可能出现 `Permission denied`；运行时传入
`--user` 是最通用的处理方式。

### 17.4 GPU 不可用

依次检查：

```bash
nvidia-smi
docker run --rm --gpus all ubuntu:22.04 nvidia-smi
docker run --rm --gpus all \
  --entrypoint python \
  <IMAGE> \
  -c 'import torch; print(torch.cuda.is_available())'
```

第二步失败时，应修复宿主机 NVIDIA Container Toolkit，而不是在镜像中重装
PyTorch。

### 17.5 镜像过大

该镜像优先保证 `test3` 的完整功能，包含 CUDA Python 包、R、Bioconductor、
编译工具和 Jupyter 组件，因此体积较大。减小体积应另建精简环境，并按功能拆分
CPU、GPU、RCTD/CellChat、pySCENIC 等镜像。不要直接从当前归档删除未知文件，
否则可能破坏二进制依赖。

### 17.6 ARM64 主机

当前归档包含 `linux-64` 二进制，不能通过 Docker Buildx 自动转换成 ARM64
可用环境。应在 `linux-aarch64` 上重新解析和安装依赖，并单独完成 R/Python
兼容性验证。

## 18. 官方资料

- Docker 构建、tag 与发布：
  https://docs.docker.com/get-started/docker-concepts/building-images/build-tag-and-publish-an-image/
- Dockerfile 构建最佳实践：
  https://docs.docker.com/build/building/best-practices/
- `docker image push`：
  https://docs.docker.com/reference/cli/docker/image/push/
- `docker login` 与 `--password-stdin`：
  https://docs.docker.com/reference/cli/docker/login/
- Docker bind mount：
  https://docs.docker.com/engine/storage/bind-mounts/
- Docker GPU：
  https://docs.docker.com/engine/containers/gpu/
- Conda 环境共享与导出：
  https://docs.conda.io/projects/conda/en/stable/user-guide/tasks/manage-environments.html
- `conda export`：
  https://docs.conda.io/projects/conda/en/stable/commands/export.html
- conda-pack：
  https://conda.github.io/conda-pack/
- NVIDIA Container Toolkit：
  https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html
- NVIDIA Docker GPU 样例：
  https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/sample-workload.html
