# Search database

## Compiling OpenFOAM docs

```bash
sudo apt install \
    doxygen \
    graphviz \
    texlive \
    texlive-latex-recommended \
    texlive-latex-extra

cp -avr /opt/openfoam13/doc openFOAM
cd openFOAM

# GENERATE_LATEX = YES

# XXX: update the Doxyfile
# doxygen -u

./Allmake latex

# cd latex
# make
```