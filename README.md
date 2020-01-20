# Basic MAJA AOT/SRE outputs analysis toolbox

MAJA outputs can be either from production C++ version or Python prototype.

## Python environment

Recommended way to have a python conda env suited to this code:
```
conda env create -f conda-scatt.yml
```
(but check the path in the end of the file)


Otherwise, the recommended settings are:
python=3.6
```
conda install -c conda-forge gdal matplotlib scipy
```

## Documentation compilation

To compile the sphinx documentation, configure sphinx with:
```
conda install sphinx sphinx_rtd_theme sphinxcontrib sphinxcontrib-applehelp sphinxcontrib-devhelp sphinxcontrib-htmlhelp sphinxcontrib-qthelp sphinxcontrib-serializinghtml sphinxcontrib-websupport
```

Then:
```
make html
```
