[![appveyor](https://ci.appveyor.com/api/projects/status/github/DTOcean/dtocean-electrical?branch=master&svg=true)](https://ci.appveyor.com/project/DTOcean/dtocean-electrical)
[![codecov](https://codecov.io/gh/DTOcean/dtocean-electrical/branch/master/graph/badge.svg)](https://codecov.io/gh/DTOcean/dtocean-electrical)
[![Lintly](https://lintly.com/gh/DTOcean/dtocean-electrical/badge.svg)](https://lintly.com/gh/DTOcean/dtocean-electrical/)
[![release](https://img.shields.io/github/release/DTOcean/dtocean-electrical.svg)](https://github.com/DTOcean/dtocean-electrical/releases/latest)

# DTOcean Electrical Sub-Systems Module

## Installation

### Install Anaconda Python

The recommended Python distribution for testing the Electical Sub-Systems module
is [Anaconda Python](https://www.continuum.io/downloads). Download the 
**Windows Python 2.7** installer and execute it.

Note: It is recommended that you install Anaconda for **All Users** (not the
default local user) and that you allow Anaconda to become the system Python
version. Other options may work, but we can not provide support if Anaconda is
installed in a different manner.

### Add Public Anaconda Cloud channel

To download some of the dependencies the following channel must be
added to Ananconda:

```
conda config --add channels https://conda.anaconda.org/topper
```

### Set up an Anaconda environment

Using a windows command prompt enter the following commands:

```
conda create -n dtocean_elec python pip pytest ipython-notebook
```

then, to activate the environment:

```
activate dtocean_elec
```

or

```
C:\Anaconda\Scripts\activate.bat dtocean_elec
```

**You must activate the environment every time you use the module.**

### Install Package Dependencies

```
conda install descartes matplotlib networkx numpy openpyxl pandas pypower scipy shapely-win-py27 xlrd xlwt
```

The package "polite" can be downloaded from [SETIS](https://setis.ec.europa.eu/dt-ocean/).

```
cd path\to\polite
winmake.bat install
```

### Install Module Package

```
cd path\to\dtocean-electrical
winmake.bat install
```

Remember to replace "path\to\" with the real path to the folder containing the
module source code.

## Testing

### Unit Tests

Unit tests will run automatically as part of the package installation process. 
They can be re-run using the command:

```
activate dtocean_elec
cd path\to\package
winmake.bat test
```

Note, you only need to activate the Anaconda environment once per session.

### Jupyter Notebooks

Examples of using the mdoule are given in [Jupyter Notebooks](http://jupyter.org/)
which are found in the "notebooks" folder of the dtocean-electrical source code.
The notebooks should be started from the Anaconda environment as follows:

```
activate dtocean_elec
start jupyter notebook
```

Note, you only need to activate the Anaconda environment once per session.

**It is important that the test data found in the "sample_data" directory is
copied into the same directory where the notebooks are being executed from**.
You can customise this directory using the config file described
[here](http://jupyter-notebook.readthedocs.io/en/latest/config.html)
and setting the "notebook_dir" variable.

Once the test data has been placed alongside the notebook, the notebook can be
executed in the normal way.

## DTOcean Project

DTOcean - "Optimal Design Tools for Ocean Energy Arrays" is funded by the 
European Commission’s 7th Framework Programme. Grant agreement number: 608597
