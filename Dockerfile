# syntax=docker/dockerfile:1

FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# ---- Pinned tool versions ----
ARG GCC_VERSION=14
ARG CMAKE_VERSION=3.29.6
ARG NINJA_VERSION=1.12.1
ARG CONAN_VERSION=2.7.0
ARG PYTHON_VENV=/opt/venv

# Keep Conan's cache outside /workspace so a bind mount over
# /workspace does not hide it.
ENV CONAN_HOME=/opt/conan

# Make the virtual environment the default Python environment.
ENV PATH="${PYTHON_VENV}/bin:${PATH}"

# ---- System dependencies ----
#
# These are installed by apt rather than allowing Conan's xorg
# recipe to install them. This makes the Docker image self-contained
# and avoids Conan system-package-manager errors.
#
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    wget \
    unzip \
    git \
    pkg-config \
    build-essential \
    python-is-python3 \
    python3 \
    python3-pip \
    python3-venv \
    # SSH client needed for fetching packages from SSH-based remotes
    openssh-client \
    \
    # GCC
    gcc-${GCC_VERSION} \
    g++-${GCC_VERSION} \
    \
    # Vulkan / OpenGL
    libvulkan-dev \
    vulkan-tools \
    mesa-vulkan-drivers \
    libgl1-mesa-dev \
    libegl-dev \
    libgles-dev \
    libglvnd-dev \
    \
    # Wayland
    libwayland-dev \
    wayland-protocols \
    \
    # X11
    libx11-dev \
    libx11-xcb-dev \
    libxext-dev \
    libxfixes-dev \
    libxi-dev \
    libxrender-dev \
    libfontenc-dev \
    libice-dev \
    libsm-dev \
    libxaw7-dev \
    libxcomposite-dev \
    libxcursor-dev \
    libxdamage-dev \
    libxinerama-dev \
    libxkbfile-dev \
    libxmu-dev \
    libxmuu-dev \
    libxpm-dev \
    libxrandr-dev \
    libxres-dev \
    libxss-dev \
    libxt-dev \
    libxtst-dev \
    libxv-dev \
    libxxf86vm-dev \
    \
    # Qt dependencies: XCB
    libxcb1-dev \
    libxcb-util-dev \
    libxcb-util0-dev \
    libxcb-glx0-dev \
    libxcb-render0-dev \
    libxcb-render-util0-dev \
    libxcb-icccm4-dev \
    libxcb-image0-dev \
    libxcb-keysyms1-dev \
    libxcb-randr0-dev \
    libxcb-shape0-dev \
    libxcb-sync-dev \
    libxcb-xfixes0-dev \
    libxcb-xinerama0-dev \
    libxcb-dri2-0-dev \
    libxcb-dri3-dev \
    libxcb-present-dev \
    libxcb-composite0-dev \
    libxcb-cursor-dev \
    libxcb-ewmh-dev \
    libxcb-res0-dev \
    libxkbcommon-dev \
    libxkbcommon-x11-dev \
    libfontconfig1-dev \
    libfreetype-dev \
    \
    && rm -rf /var/lib/apt/lists/*

# ---- Set GCC as the default compiler ----
RUN update-alternatives --install \
        /usr/bin/gcc gcc /usr/bin/gcc-${GCC_VERSION} 100 \
    && update-alternatives --install \
        /usr/bin/g++ g++ /usr/bin/g++-${GCC_VERSION} 100 \
    && update-alternatives --install \
        /usr/bin/cc cc /usr/bin/gcc-${GCC_VERSION} 100 \
    && update-alternatives --install \
        /usr/bin/c++ c++ /usr/bin/g++-${GCC_VERSION} 100

# ---- Install pinned CMake ----
RUN ARCH="$(uname -m)" \
    && wget -q \
        "https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-${ARCH}.sh" \
        -O /tmp/cmake.sh \
    && sh /tmp/cmake.sh --skip-license --prefix=/usr/local \
    && rm /tmp/cmake.sh \
    && cmake --version

# ---- Install pinned Ninja ----
RUN wget -q \
        "https://github.com/ninja-build/ninja/releases/download/v${NINJA_VERSION}/ninja-linux.zip" \
        -O /tmp/ninja.zip \
    && unzip -o /tmp/ninja.zip -d /usr/local/bin \
    && chmod +x /usr/local/bin/ninja \
    && rm /tmp/ninja.zip \
    && ninja --version

# ---- Install pinned Conan ----
RUN python3 -m venv "${PYTHON_VENV}" \
    && "${PYTHON_VENV}/bin/pip" install --no-cache-dir --upgrade pip \
    && "${PYTHON_VENV}/bin/pip" install --no-cache-dir "conan==${CONAN_VERSION}" \
    && conan --version

# ---- Configure Conan ----
RUN mkdir -p "${CONAN_HOME}" \
    && conan profile detect --force \
    && CONAN_PROFILE="${CONAN_HOME}/profiles/default" \
    && sed -i 's/^compiler.version=.*/compiler.version=14/' "${CONAN_PROFILE}" \
    && sed -i 's/^compiler.libcxx=.*/compiler.libcxx=libstdc++11/' "${CONAN_PROFILE}" \
    && printf '\n[conf]\n' >> "${CONAN_PROFILE}" \
    && printf 'tools.cmake.cmaketoolchain:generator=Ninja\n' >> "${CONAN_PROFILE}" \
    && printf 'tools.system.package_manager:mode=check\n' >> "${CONAN_PROFILE}" \
    && conan remote update conancenter --url https://center2.conan.io \
    && conan remote list \
    && cat "${CONAN_PROFILE}"

WORKDIR /workspace

# ---- Pre-cache Conan dependencies ----
COPY conanfile.py .

RUN conan install . \
    --build=missing \
    -s build_type=Release

# ---- Entrypoint ----
COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh

RUN chmod +x /usr/local/bin/docker-entrypoint.sh

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

CMD ["/bin/bash"]
