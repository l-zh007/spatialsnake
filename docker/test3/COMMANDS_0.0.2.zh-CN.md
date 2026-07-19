# Spatialsnake 0.0.2 test3 Docker 构建、推送、拉取命令

以下命令假设项目根目录为：

```bash
cd /users/zhenghlin/beifen2/spatialdata
```

## 1. 可选：重新打包 test3 环境

```bash
CONDA_EXE=/users/zhenghlin/miniconda3/bin/conda \
CONDA_ENV_PREFIX=/users/zhenghlin/miniconda3/envs/test3 \
CONDA_PACK_EXE=/users/zhenghlin/miniconda3/envs/snakemake_env/bin/conda-pack \
COMPRESS_LEVEL=1 \
./docker/test3/export-test3.sh
```

## 2. 构建本地镜像

```bash
cd /users/zhenghlin/beifen2/spatialdata

IMAGE_NAME=spatialsnake-test3 \
IMAGE_TAG=0.0.2-test3-20260719 \
APP_UID="$(id -u)" \
APP_GID="$(id -g)" \
./docker/test3/build-image.sh
```

等价手动构建命令：

```bash
cd /users/zhenghlin/beifen2/spatialdata

DOCKER_BUILDKIT=1 docker build \
  --pull \
  --file docker/test3/Dockerfile \
  --tag spatialsnake-test3:0.0.2-test3-20260719 \
  --build-arg IMAGE_VERSION=0.0.2-test3-20260719 \
  --build-arg BUILD_DATE="$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
  --build-arg VCS_REF="$(git rev-parse --short=12 HEAD)" \
  --build-arg APP_UID="$(id -u)" \
  --build-arg APP_GID="$(id -g)" \
  .
```

## 3. 登录 Docker Hub

交互式登录：

```bash
docker login -u YOUR_DOCKERHUB_USER
```

Token 登录：

```bash
printf '%s' "$DOCKERHUB_TOKEN" \
  | docker login --username YOUR_DOCKERHUB_USER --password-stdin
```

## 4. 推送到 Docker Hub

使用脚本：

```bash
cd /users/zhenghlin/beifen2/spatialdata

REGISTRY=docker.io \
REGISTRY_NAMESPACE=YOUR_DOCKERHUB_USER \
REGISTRY_REPOSITORY=spatialsnake-test3 \
IMAGE_NAME=spatialsnake-test3 \
IMAGE_TAG=0.0.2-test3-20260719 \
./docker/test3/push-image.sh
```

等价手动推送命令：

```bash
docker image tag \
  spatialsnake-test3:0.0.2-test3-20260719 \
  YOUR_DOCKERHUB_USER/spatialsnake-test3:0.0.2-test3-20260719

docker image push \
  YOUR_DOCKERHUB_USER/spatialsnake-test3:0.0.2-test3-20260719
```

可选推送 latest：

```bash
docker image tag \
  YOUR_DOCKERHUB_USER/spatialsnake-test3:0.0.2-test3-20260719 \
  YOUR_DOCKERHUB_USER/spatialsnake-test3:latest

docker image push YOUR_DOCKERHUB_USER/spatialsnake-test3:latest
```

## 5. 目标机器拉取

```bash
docker pull \
  YOUR_DOCKERHUB_USER/spatialsnake-test3:0.0.2-test3-20260719
```

私有仓库先登录：

```bash
docker login -u YOUR_DOCKERHUB_USER

docker pull \
  YOUR_DOCKERHUB_USER/spatialsnake-test3:0.0.2-test3-20260719
```

按 digest 固定拉取：

```bash
docker pull \
  YOUR_DOCKERHUB_USER/spatialsnake-test3@sha256:<REMOTE_DIGEST>
```

## 6. 运行示例

```bash
cd "$HOME/spatialsnake-run"

docker run --rm \
  --mount type=bind,src="$PWD",dst=/work \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  YOUR_DOCKERHUB_USER/spatialsnake-test3:0.0.2-test3-20260719 \
  single_analysis sample.txt visium \
  --option=clustering \
  --jobs=16
```
