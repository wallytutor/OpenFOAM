#!/usr/bin/env bash

cd ${0%/*} || exit 1

# Temporary directory:
mkdir -p unused/

# Generate the global thermodynamics file:
uv run yaml2ck                           \
    --phase-name    gri30                \
    --thermo        thermo.dat           \
    --mechanism     unused/reactions.inp \
    --transport     unused/transport.dat \
    --sort-elements alphabetical         \
    --sort-species  alphabetical         \
    --overwrite                          \
    --validate                           \
    gri30.yaml


# Generate transport data:
if [ ! -f "transport" ]; then
    uv run python convert.py
fi

# Generate YAML file:
if [ ! -f "wb.yaml" ]; then
    uv run ck2yaml -v \
        --input reactions.inp \
        --thermo thermo.dat \
        --transport unused/transport.dat \
        --name gas \
        --output wb.yaml
fi

# XXX after ck2yaml (which requires the default ranges in thermo.dat)!
# Patches for compatibility with chemkinToFoam:
sed -i '/200.000   1000.000  6000.000/d' thermo.dat
sed -i 's|G300.000|G200.000|g' thermo.dat
sed -i 's|REACTIONS CAL/MOLE MOLE|REACTIONS CAL/MOLE MOLES|g' reactions.inp

# Convert to OpenFOAM format (HERE!!!):
chemkinToFoam     \
    reactions.inp \
    thermo.dat    \
    transport     \
    reactions     \
    thermo

# Clean temporary directory:
rm -rf unused/