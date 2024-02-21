#!/bin/env bash

# Run from strictly this directory.
cd ${0%/*} || exit 1

# Ensure clean start.
rm -rf case/

# Copy reference case here.
cp -avr ../reference case/

echo """
# Fix patches types.
foamDictionary constant/polyMesh/boundary \
    -entry entry0/frontAndBack/type -set empty

foamDictionary constant/polyMesh/boundary \
    -entry entry0/walls/type        -set wall
""" > case/patches.sh

# Change directory and run case.
cd case && chmod u+x Allrun.sh && ./Allrun.sh
