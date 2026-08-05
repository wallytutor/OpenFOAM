# Search database

## OpenFOAM *CFD Book*

```bash
# Example recursive crawl converting HTML pages to a clean text/PDF bundle
wget --recursive --level=2 --no-clobber --page-requisites --convert-links --no-parent http://example.com
```

## OpenFOAM Documentation

```bash
# Install required tools:
sudo apt install \
    doxygen \
    graphviz \
    pandoc \
    texlive \
    texlive-latex-recommended \
    texlive-latex-extra

# Copy documentation directory:
cp -avr /opt/openfoam13/doc docs/ && cd docs/Doxygen

# Enable LaTeX generation:
sed -i 's/\(GENERATE_LATEX[[:space:]]*=[[:space:]]*\)NO/\1YES/' Doxyfile

# Update the Doxyfile to latest standard:
doxygen -u

# Generate the documentation:
./Allmake

# Build final PDF:
cd latex && make
```

## Gmsh Documentation
