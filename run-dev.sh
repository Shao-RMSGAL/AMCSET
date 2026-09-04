#!/usr/bin/env bash
#
# run-dev.sh — Build, start, and enter the AMCSET dev container with
#              GPU-accelerated rendering and GUI forwarding.
#
# Usage:
#   ./run-dev.sh            # build image if needed, then enter the container
#   ./run-dev.sh --rebuild  # force a rebuild of the image first
#
set -euo pipefail

IMAGE_NAME="amcset-dev"
CONTAINER_NAME="amcset-dev-container"
WORKDIR_IN_CONTAINER="/workspace"
CONAN_CACHE_VOLUME="amcset-conan-cache"   # persists deps across runs

# ---- Optional forced rebuild ----
if [[ "${1:-}" == "--rebuild" ]]; then
    echo "Forcing image rebuild..."
    docker build --no-cache -t "${IMAGE_NAME}" .
    shift
elif ! docker image inspect "${IMAGE_NAME}" >/dev/null 2>&1; then
    echo "Image '${IMAGE_NAME}' not found. Building..."
    docker build -t "${IMAGE_NAME}" .
fi

# ---- GPU flags ----
GPU_FLAGS=()
if [[ -d /dev/dri ]]; then
    echo "Found /dev/dri — enabling hardware-accelerated rendering."
    GPU_FLAGS+=(--device /dev/dri)
    if getent group render >/dev/null 2>&1; then
        GPU_FLAGS+=(--group-add "$(getent group render | cut -d: -f3)")
    fi
    if getent group video >/dev/null 2>&1; then
        GPU_FLAGS+=(--group-add "$(getent group video | cut -d: -f3)")
    fi
else
    echo "Warning: /dev/dri not found. Falling back to software rendering."
fi

# ---- Display / GUI passthrough ----
DISPLAY_FLAGS=()

# Note: with `set -u`, guard every host variable with :- to avoid aborting.
HOST_XDG_RUNTIME_DIR="${XDG_RUNTIME_DIR:-}"
HOST_WAYLAND_DISPLAY="${WAYLAND_DISPLAY:-}"
HOST_DISPLAY="${DISPLAY:-}"

if [[ -n "${HOST_WAYLAND_DISPLAY}" && -n "${HOST_XDG_RUNTIME_DIR}" \
      && -S "${HOST_XDG_RUNTIME_DIR}/${HOST_WAYLAND_DISPLAY}" ]]; then
    # ----- Wayland session -----
    echo "Detected Wayland display: ${HOST_WAYLAND_DISPLAY}"
    DISPLAY_FLAGS+=(
        -e XDG_RUNTIME_DIR=/tmp
        -e WAYLAND_DISPLAY="${HOST_WAYLAND_DISPLAY}"
        -v "${HOST_XDG_RUNTIME_DIR}/${HOST_WAYLAND_DISPLAY}:/tmp/${HOST_WAYLAND_DISPLAY}"
        -e QT_QPA_PLATFORM="wayland;xcb"
    )
    # X11 fallback for XWayland apps, if an X socket exists
    if [[ -n "${HOST_DISPLAY}" && -d /tmp/.X11-unix ]]; then
        DISPLAY_FLAGS+=(
            -e DISPLAY="${HOST_DISPLAY}"
            -v /tmp/.X11-unix:/tmp/.X11-unix
        )
    fi
elif [[ -n "${HOST_DISPLAY}" ]]; then
    # ----- X11 session -----
    echo "Detected X11 display: ${HOST_DISPLAY}"
    if command -v xhost >/dev/null 2>&1; then
        xhost +local:docker >/dev/null 2>&1 || \
            echo "Warning: 'xhost +local:docker' failed; GUI may not appear."
    else
        echo "Warning: 'xhost' not available; GUI may not appear."
    fi
    DISPLAY_FLAGS+=(
        -e DISPLAY="${HOST_DISPLAY}"
        -e QT_QPA_PLATFORM=xcb
        -v /tmp/.X11-unix:/tmp/.X11-unix
    )
else
    echo "Warning: No usable WAYLAND_DISPLAY or DISPLAY detected."
    echo "The GUI will not appear and conanfile display detection may fail."
fi

# ---- Persist the Conan cache across container runs ----
# The image sets CONAN_HOME=/opt/conan. Mounting a named volume there means
# dependencies downloaded/built inside the container are reused next time.
docker volume inspect "${CONAN_CACHE_VOLUME}" >/dev/null 2>&1 || \
    docker volume create "${CONAN_CACHE_VOLUME}" >/dev/null

# ---- Remove any stale container with the same name ----
echo "Removing any stale containers..."
docker rm -f "${CONTAINER_NAME}" >/dev/null 2>&1 || true
echo "Done removing."

# ---- Start and enter the container ----
echo "Starting docker container..."
docker run -d \
    --name "${CONTAINER_NAME}" \
    "${GPU_FLAGS[@]}" \
    "${DISPLAY_FLAGS[@]}" \
    -v "$(pwd):${WORKDIR_IN_CONTAINER}" \
    -v "${CONAN_CACHE_VOLUME}:/opt/conan" \
    -w "${WORKDIR_IN_CONTAINER}" \
    "${IMAGE_NAME}" \
    sleep infinity

echo "Container started."
echo "Attach to it from VS Code using:"
echo "    Dev Containers: Attach to Running Container"
