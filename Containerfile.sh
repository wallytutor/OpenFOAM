#!/usr/bin/env bash

# _name_ of Containerfile._name_:
name="foam"

# Name the project based on the current user:
project="${name}-$(whoami)"

# Build the container image:
podman build -t "${project}" -f "Containerfile" .

# Dump container to portable .tar file:
podman save -o "${project}.tar" "localhost/${project}"

# Convert container into apptainer:
apptainer build "${project}.sif" "docker-archive://${project}.tar"
