# 使用 GitHub Actions 构建并推送 test3 镜像

适用场景：服务器没有 sudo 权限，rootless Docker 也无法稳定连接 Docker Hub，但可以把代码和大文件分片上传到 GitHub。

## 方案概要

```text
提交 docker/test3 构建文件和 workflow
        |
        +-- 把 test3-linux-64.tar.gz 分片上传到 GitHub Release
        |
        +-- GitHub Actions 下载分片并合并
        |
        +-- GitHub Actions 构建 Docker 镜像
        |
        +-- GitHub Actions 登录 Docker Hub 并推送镜像
```

注意：`test3-linux-64.tar.gz` 约 6.2 GB，不能直接提交到 Git，也不能作为单个 GitHub Release asset 上传。需要分片。

## 1. 网页端配置 Docker Hub Secrets

在 GitHub 仓库页面进入：

```text
Settings -> Secrets and variables -> Actions -> New repository secret
```

添加两个 secret：

```text
DOCKERHUB_USERNAME = jacking007
DOCKERHUB_TOKEN    = Docker Hub Personal Access Token
```

Docker Hub PAT 在 Docker Hub 网页创建：

```text
Account settings -> Personal access tokens
```

## 2. 本地/服务器生成分片

```bash
cd /users/zhenghlin/beifen2/spatialdata

./docker/test3/split-test3-artifact.sh
```

默认每片 `1900M`，低于 GitHub Release 单文件 2 GiB 限制。

## 3. 创建 GitHub Release 并上传分片

如果服务器有 `gh`：

```bash
export HTTPS_PROXY="${https_proxy:-${HTTPS_PROXY:-}}"
export HTTP_PROXY="${http_proxy:-${HTTP_PROXY:-}}"

gh auth login

gh release create test3-env-0.0.2-20260719 \
  --repo l-zh007/spatialsnake \
  --title "test3 packed environment 0.0.2 20260719" \
  --notes "Packed conda environment chunks for GitHub Actions Docker build."

gh release upload test3-env-0.0.2-20260719 \
  docker/test3/artifacts/test3-linux-64.tar.gz.part-* \
  docker/test3/artifacts/test3-linux-64.tar.gz.parts.sha256 \
  --repo l-zh007/spatialsnake \
  --clobber
```

如果没有 `gh`，使用 GitHub 网页：

```text
Releases -> Draft a new release
Tag: test3-env-0.0.2-20260719
Title: test3 packed environment 0.0.2 20260719
Attach binaries:
  docker/test3/artifacts/test3-linux-64.tar.gz.part-000
  docker/test3/artifacts/test3-linux-64.tar.gz.part-001
  docker/test3/artifacts/test3-linux-64.tar.gz.part-002
  docker/test3/artifacts/test3-linux-64.tar.gz.part-003
  docker/test3/artifacts/test3-linux-64.tar.gz.parts.sha256
Publish release
```

## 4. 提交 workflow 和 docker/test3 构建文件

```bash
cd /users/zhenghlin/beifen2/spatialdata

git add .gitignore .github/workflows/dockerhub-test3.yml docker/test3
git status --short
git commit -m "ci: build and publish test3 docker image"
git push origin main
```

不要提交 `docker/test3/artifacts/test3-linux-64.tar.gz` 或分片；它们由 `.gitignore` 忽略。

## 5. 网页端触发 GitHub Actions

进入 GitHub 仓库页面：

```text
Actions -> Build and Push test3 Docker Image -> Run workflow
```

填写：

```text
release_tag = test3-env-0.0.2-20260719
image_tag   = 0.0.2-test3-20260719
push_latest = false
```

点击 `Run workflow`。

## 6. 拉取镜像

Actions 成功后：

```bash
docker pull jacking007/spatialsnake-test3:0.0.2-test3-20260719
```

如果 Docker Hub 仓库是 private，先登录：

```bash
docker login -u jacking007
docker pull jacking007/spatialsnake-test3:0.0.2-test3-20260719
```
