#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Ensure remove case base files.
rm -rf 0/ constant/ system/ All*

# Also remove files that might be generated at runtime.
rm -rf [0-9]*
rm -rf constant/polyMesh/
rm -rf logging
rm -rf postProcessing
rm -rf processors*
rm -rf VTK

# Copy all files from default case
cp -avr ../case/* .

# Rename zero/ directory to 0/
mv -f zero/ 0/ && rm -rf zero/

# Copy specific cloudProperties to constant/
cp cloudProperties constant/
