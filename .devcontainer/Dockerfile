# .devcontainer/Dockerfile

# Base C++ image
FROM mcr.microsoft.com/devcontainers/cpp:1-debian-12

# Install FFmpeg (for video) and TBB (for C++ parallel algorithms)
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends ffmpeg libtbb-dev

# Set the default working directory for when the container runs
WORKDIR /workspace