[![appveyor](https://ci.appveyor.com/api/projects/status/github/DTOcean/dtocean-electrical?branch=master&svg=true)](https://ci.appveyor.com/project/DTOcean/dtocean-electrical)
[![codecov](https://codecov.io/gh/DTOcean/dtocean-electrical/branch/master/graph/badge.svg)](https://codecov.io/gh/DTOcean/dtocean-electrical)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/64b45b67ff304732be79435c3fb20751)](https://www.codacy.com/project/H0R5E/dtocean-electrical/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=DTOcean/dtocean-electrical&amp;utm_campaign=Badge_Grade_Dashboard&amp;branchId=11704041)
[![release](https://img.shields.io/github/release/DTOcean/dtocean-electrical.svg)](https://github.com/DTOcean/dtocean-electrical/releases/latest)

# DTOcean Electrical Sub-Systems Module

This package provides the Electrical Sub-Systems design module for the DTOcean 
tools. It can design the electrical network of an array of fixed or floating 
wave or tidal ocean energy converters (constrained by the environment and 
device design), and calculate the cost and electrical losses. It optimises the 
design for minimum cost per unit exported power. 

See [dtocean-app](https://github.com/DTOcean/dtocean-app) or [dtocean-core](
https://github.com/DTOcean/dtocean-app) to use this package within the DTOcean
ecosystem.

* For python 2.7 only.

## Installation

Installation and development of dtocean-electrical uses the [Anaconda 
Distribution](https://www.anaconda.com/distribution/) (Python 2.7)

### Conda Package

To install:

```
$ conda install -c dataonlygreater dtocean-electrical
```

### Source Code

Conda can be used to install dependencies into a dedicated environment from
the source code root directory:

```
$ conda create -n _dtocean_electro python=2.7 pip
```

Activate the environment, then copy the `.condrc` file to store installation  
channels:

```
$ conda activate _dtocean_electro
$ copy .condarc %CONDA_PREFIX%
```

Install [polite](https://github.com/DTOcean/polite) into the environment. For 
example, if installing it from source:

```
$ cd \\path\\to\\polite
$ conda install --file requirements-conda-dev.txt
$ pip install -e .
```

Finally, install dtocean-electrical and its dependencies using conda and pip:

```
$ cd \\path\\to\\dtocean-electrical
$ conda install --file requirements-conda-dev.txt
$ pip install -e .
```

To deactivate the conda environment:

```
$ conda deactivate
```

### Tests

A test suite is provided with the source code that uses [pytest](
https://docs.pytest.org).

If not already active, activate the conda environment set up in the [Source 
Code](#source-code) section:

```
$ conda activate _dtocean_electro
```

Install packages required for testing to the environment (one time only):

```
$ conda install -y pytest pytest-mock openpyxl xlrd xlwt
```

Run the tests:

``` 
$ py.test tests
```

### Uninstall

To uninstall the conda package:

```
$ conda remove dtocean-electrical
```

To uninstall the source code and its conda environment:

```
$ conda remove --name _dtocean_electro --all
```

## Usage

Example scripts are available in the "examples" folder of the source code.

```
$ cd examples
$ python electrical_run.py
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first to
discuss what you would like to change.

See [this blog post](
https://www.dataonlygreater.com/latest/professional/2017/03/09/dtocean-development-change-management/)
for information regarding development of the DTOcean ecosystem.

Please make sure to update tests as appropriate.

## Credits

This package was initially created as part of the [EU DTOcean project](
https://www.dtoceanplus.eu/About-DTOceanPlus/History) by:

 * Adam Collin at [the University of Edinburgh](https://www.ed.ac.uk/)
 * Mathew Topper at [TECNALIA](https://www.tecnalia.com)

It is now maintained by Mathew Topper at [Data Only Greater](
https://www.dataonlygreater.com/).

## License

[GPL-3.0](https://choosealicense.com/licenses/gpl-3.0/)
