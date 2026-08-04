#!/usr/bin/env bash

# Run from this directory:
cd ${0%/*} || exit 1

# To open with ParaView:
touch case.foam

# Convert base Gmsh mesh to OpenFOAM:
gmshToFoam wedge-coaxial-nozzle.msh

# Perform preliminary renumbering:
renumberMesh

# Check global mesh:
checkMesh -allTopology

# Fix patch types after splitting:
file=constant/polyMesh/boundary
foamDictionary $file -entry entry0/back/type  -set wedge
foamDictionary $file -entry entry0/front/type -set wedge
foamDictionary $file -entry entry0/walls/type -set wall

# Check mesh:
checkMesh
