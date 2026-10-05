#!/usr/bin/env bash

cd ${0%/*} || exit 1

# Path to output directory:
OUTPUT="../model/constant"

# Temporary directory:
mkdir -p unused/

# Generate the global thermodynamics file:
if [ ! -f "thermo.dat" ]; then
    uv run yaml2ck                            \
        --phase-name    gri30                 \
        --thermo        thermo.dat            \
        --mechanism     unused/_reactions.inp \
        --transport     unused/_transport.dat \
        --sort-elements alphabetical          \
        --sort-species  alphabetical          \
        --overwrite                           \
        --validate                            \
        gri30.yaml

    # Patches for compatibility with chemkinToFoam:
    sed -i '/200.000   1000.000  6000.000/d' thermo.dat
    sed -i 's|G300.000|G200.000|g' thermo.dat
    sed -i 's|REACTIONS CAL/MOLE MOLE|REACTIONS CAL/MOLE MOLES|g' reactions.inp
fi

# Clean temporary directory:
rm -rf unused/

# Generate transport data:
if [ ! -f "transport" ]; then
    uv run python convert.py
fi

# Convert to OpenFOAM format:
chemkinToFoam         \
    reactions.inp     \
    thermo.dat        \
    transport         \
    $OUTPUT/reactions \
    $OUTPUT/thermo
