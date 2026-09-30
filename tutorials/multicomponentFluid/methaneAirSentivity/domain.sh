#!/usr/bin/env bash

# Run from this directory:
cd ${0%/*} || exit 1

# Source tutorial run functions:
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# To open with ParaView:
touch case.foam

# Generate the base mesh:
uv run python domain.py > log.mesh

# Convert base Gmsh mesh to OpenFOAM:
runApplication gmshToFoam domain.msh

# Perform preliminary renumbering:
runApplication renumberMesh

# Fix patch types after splitting:
file=constant/polyMesh/boundary
foamDictionary $file -entry entry0/back/type  -set wedge > /dev/null
foamDictionary $file -entry entry0/front/type -set wedge > /dev/null
foamDictionary $file -entry entry0/walls/type -set wall  > /dev/null

# Check mesh:
runApplication checkMesh

#------------------------------------------------------------------------------