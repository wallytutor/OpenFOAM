#!/usr/bin/env bash

# Run from this directory:
cd ${0%/*} || exit 1

# To open with ParaView:
touch case.foam

# Convert base Gmsh mesh to OpenFOAM:
gmshToFoam wedge-horizontal-cht.msh

# Perform preliminary renumbering:
renumberMesh

# Check global mesh:
checkMesh -allTopology

# Split regions:
splitMeshRegions -cellZones

# Clean original converted mesh:
rm -rf constant/polyMesh

# Fix patch types after splitting:
file=constant/solid/polyMesh/boundary
foamDictionary $file -entry entry0/solidFront/type -set wedge
foamDictionary $file -entry entry0/solidBack/type  -set wedge
foamDictionary $file -entry entry0/wallSolid/type  -set wall

file=constant/fluid/polyMesh/boundary
foamDictionary $file -entry entry0/fluidFront/type -set wedge
foamDictionary $file -entry entry0/fluidBack/type  -set wedge
foamDictionary $file -entry entry0/surfFluid/type  -set wall

# Check mesh:
checkMesh -region solid
checkMesh -region fluid
