#!/usr/bin/env bash

cd ${0%/*} || exit 1

if [ ! -f "BFER_methane.yaml" ]; then
    url=https://chemistry.cerfacs.fr/chemistry-repo/pub/mechanisms
    wget --no-check-certificate $url/BFER_methane/BFER_methane.yaml
fi

if [ ! -f "reactions.inp" ] || [ ! -f "thermo.dat" ]; then
    uv run yaml2ck                    \
        --phase-name    CH4_BFER_multi \
        --mechanism     reactions.inp  \
        --thermo        thermo.dat     \
        --transport     transport.dat  \
        --sort-elements alphabetical   \
        --sort-species  alphabetical   \
        --overwrite                    \
        --validate                     \
        BFER_methane.yaml

    # Patches for compatibility with chemkinToFoam:
    sed -i '/200.000   1000.000  5000.000/d' thermo.dat
    sed -i 's|G300.000|G200.000|g' thermo.dat
    sed -i 's|REACTIONS CAL/MOLE MOLE|REACTIONS CAL/MOLE MOLES|g' reactions.inp
fi

# Generate transport data:
uv run python convert.py

# For BFER use this:
chemkinToFoam             \
    reactions.inp         \
    thermo.dat            \
    transport             \
    ../constant/reactions \
    ../constant/thermo

# For Westbrook use this:
# chemkinToFoam             \
#     WB_2step.inp          \
#     WB_2step.dat          \
#     transport             \
#     ../constant/reactions \
#     ../constant/thermo
