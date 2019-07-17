# Very basic scatter plot utility to cross compare MAJA AOT outputs

**Read the doc in: ```_build/html/index.html```**

```
Usage: scatt.py [-h] [-v] [-r n] runA runB
positional arguments:
  runA                  XML file describing run A
  runB                  XML file describing run A
optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Set verbosity to INFO level
  -r, --resampling      Integer n such that products are resampled to 1/n
Example:
    ./scatt.py M01_MAJA.xml M01_MAQT.xml -r 12
```
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
