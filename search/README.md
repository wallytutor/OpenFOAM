# Search database

## Compiling OpenFOAM docs

```bash
sudo apt install \
    doxygen \
    texlive \
    texlive-latex-recommended \
    texlive-latex-extra

cp -avr /opt/openfoam13/doc openFOAM
cd openFOAM

# GENERATE_LATEX = YES

./Allmake latex

# cd latex
# make
```