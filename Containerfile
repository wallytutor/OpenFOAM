# ############################################################################
#  Containerfile.foam
# ############################################################################

FROM ubuntu:24.04

# ----------------------------------------------------------------------------
# Basics and avoid locale issues
# ----------------------------------------------------------------------------

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y locales && \
    locale-gen en_US.UTF-8

ENV LANG=en_US.UTF-8     \
    LANGUAGE=en_US:en    \
    LC_CTYPE=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

RUN update-locale LANG=en_US.UTF-8

# ----------------------------------------------------------------------------
# APT installed tools
# ----------------------------------------------------------------------------

# If the image build breaks (stops installing some of the packages), then it
# may be that the database needs to be updated first; place the failing command
# just after a(nother) copy of the following line:
RUN apt-get update

# Base development toolkit (may be required by Python packages):
RUN apt-get install -y build-essential git gcc g++ gfortran make cmake
RUN apt-get install -y libboost-dev liblapack-dev libopenblas-dev

# Install extra apt packages:
# - software-properties-common for adding repositories
# - wget for getting keys
# - curl for downloading Rust
# - neovim for code editing
# - libssl-dev required by typst
# - pkg-config required by typst
RUN apt-get install -y \
    software-properties-common \
    wget \
    curl \
    neovim \
    libssl-dev \
    pkg-config

# Install OpenFOAM:
RUN sh -c "wget -O - https://dl.openfoam.org/gpg.key > /etc/apt/trusted.gpg.d/openfoam.asc"
RUN add-apt-repository http://dl.openfoam.org/ubuntu
RUN apt-get update && apt-get install -y openfoam13

# Permanently source OpenFOAM environment:
RUN echo "source /opt/openfoam13/etc/bashrc" >> /etc/bash.bashrc

# ----------------------------------------------------------------------------
# For Python/Rust
# ----------------------------------------------------------------------------

# Install uv from the official image
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

ENV CARGO_HOME=/opt/cargo
ENV RUSTUP_HOME=/opt/rustup

RUN curl https://sh.rustup.rs -sSf | sh -s -- -y --no-modify-path \
    --default-toolchain stable \
    --profile minimal \
    --default-host x86_64-unknown-linux-gnu

# Add Cargo to PATH for all container sessions:
ENV PATH="/opt/cargo/bin:$PATH"

RUN cargo install maturin
RUN cargo install --locked typst-cli

ARG QUARTO_URL=https://github.com/quarto-dev/quarto-cli/releases/download
ARG QUARTO_VERSION=1.10.18
ARG QUARTO_DEB=quarto-${QUARTO_VERSION}-linux-amd64.deb

RUN wget ${QUARTO_URL}/v${QUARTO_VERSION}/${QUARTO_DEB} \
    && dpkg -i ${QUARTO_DEB} \
    && rm ${QUARTO_DEB}

# ----------------------------------------------------------------------------
# FINAL STEPS
# ----------------------------------------------------------------------------

# Clean up and make image smaller:
RUN apt-get clean  && rm -rf /var/lib/apt/lists/* && rm -rf /var/tmp/build*

WORKDIR /home/ubuntu

############################################################################
# EOF
############################################################################