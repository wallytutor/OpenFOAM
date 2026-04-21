#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

handleOpenFOAM() {
    if [ ! -d /opt/openfoam13 ]; then
        local repo="dl.openfoam.org"
        local keyasc="/etc/apt/trusted.gpg.d/openfoam.asc"

        sudo sh -c "wget -O - https://$repo/gpg.key > $keyasc"
        sudo add-apt-repository http://$repo/ubuntu

        sudo apt-get update
        sudo apt-get install -y openfoam13
    fi

    source /opt/openfoam13/etc/bashrc
}

handleUvInstall() {
    if command -v uv >/dev/null 2>&1; then
        echo "Already installed: uv $(uv --version)"
    else
        curl -LsSf https://astral.sh/uv/install.sh | sh
        source $HOME/.local/bin/env
    fi
}

handleUvMode() {
    # In WSL, if not working in /home/<user>/, the the following before
    # running `uv` to avoid symlink issues:
    if grep -qi "microsoft" /proc/version; then
        if realpath "$PWD" | grep -q "^$(realpath "$HOME")"; then
            echo "Inside HOME"
        else
            echo "Running in WSL outside HOME: setting UV_LINK_MODE=copy"
            export UV_LINK_MODE=copy
        fi
    fi
}

handleEnvironment() {
    local targetEnv="${SCRIPT_DIR}/.venv"

    if [ -n "$VIRTUAL_ENV" ] && [ "$VIRTUAL_ENV" != "$targetEnv" ]; then
        echo "Inside a virtual environment: $(which python)"
        echo "Deactivating current virtual environment"
        deactivate
    fi

    if [ -d "$targetEnv" ]; then
        echo "Activating existing virtual environment"
        source "$targetEnv/bin/activate"
    else
        echo "Creating new virtual environment with uv"
        uv venv --python 3.12 "$targetEnv"
        source "$targetEnv/bin/activate"
        uv pip install -r requirements.txt
    fi
}

# Install `uv` if not already installed:
handleUvInstall

# Handle WSL mode if applicable:
handleUvMode

# Handle virtual environment:
handleEnvironment
