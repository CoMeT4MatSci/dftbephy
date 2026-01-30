# Installation

## Prerequisites

- NumPy and SciPy (`pip install numpy scipy`)
- DFTB+ ( For detailed installation instructions, please visit [DFTB+ Recipes](https://dftbplus-recipes.readthedocs.io/en/latest/introduction.html).)
- [Phonopy](https://phonopy.github.io/phonopy/install.html) (`conda install -c conda-forge phonopy`)
- [HSD](https://github.com/dftbplus/hsd-python) for DFTBephy input (`hsd-python` is available via `conda`)
- Cython for faster routines (`pip install Cython`)
- Other packages might be necessary: mpi4py, openmpi, spglib, h5py (Available via `pip` and `conda`)


> **Recommended installation method:**<br>
While not strictly required, we recommend using a conda environment to keep dependencies isolated, avoid conflicts with system-wide packages, and make the setup reproducible. <br>
Using `conda` or `mamba` is often the smoothest way to install packages because it provides precompiled binaries and handles compiled dependencies consistently. In this tutorial we preferred `mamba` since it resolves and installs dependencies faster than `conda`. Use `pip` only for packages not available on `conda-forge`.

## How to install DFTBephy

- Pull / download latest version from [GitHub](https://github.com/CoMeT4MatSci/dftbephy).
- Run `python setup.py build_ext --inplace` in the terminal to build faster routines.
- Run `pip install -e .` in the terminal to install package but keep it editable in the current directory.

## Usage in a Conda environment

If conda is not already installed, you can download the appropriate installer from [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main) and then follow these steps:

```
cd ~/Downloads
bash Miniconda3-latest-Linux-x86_64.sh
source ~/miniconda3/bin/activate
```

After running the last command, the base conda environment should be active—you should see (base) at the beginning of your terminal prompt. Before installing packages, make sure the conda-forge channel is enabled (needed if you want to install DFTB+ from conda-forge):

```
conda config --add channels conda-forge
```

Setting `channel_priority` to `strict` makes conda resolve packages as consistently as possible from a single, highest-priority channel. This reduces the chance of mixing packages from `defaults` and `conda-forge`, which can otherwise lead to version or dependency incompatibilities.

```
conda config --set channel_priority strict
```

If you prefer not to change global conda settings (or if `channel_priority strict` is not set), you can explicitly select the channel per command, e.g. `mamba install -c conda-forge numpy`.


### Environment Setup for DFTBephy

The required Phonopy version (`phonopy<2.41`) is not available as a conda build for newer Python versions. To avoid this, create a new conda environment using Python 3.10, then activate it.

```
conda create --name dftbephy "python=3.10"
conda activate dftbephy
```

Optionally, you can install DFTB+ inside this conda environment. If DFTB+ is already available on the machine, you can skip this step.

```
mamba install 'dftbplus=*=nompi_*'
mamba install  dftbplus-tools dftbplus-python
```

When `dftbephy` environment is still active, install the prerequisites for DFTBephy:

```
mamba install numpy scipy cython spglib openmpi mpi4py h5py setuptools hsd-python "phonopy<2.41"
```

Unfortunately, not all versions of phonopy are compatible with the current version of DFTBephy. The latest phonopy version that is still compatible is **Phonopy Version 2.40.0**.


> **Note:** Phonopy version 2.40.0 works with deprecation warnings:<br>
In the latest phonopy versions, some of the classes have been updated (see the phonopy changelog [here](https://phonopy.github.io/phonopy/changelog.html#)).
Some of our older scripts still use the legacy `get_*` calls and may require minor edits in order to run without warnings or errors.
Updating these calls is a part of the ongoing work in progress in our code base.

After downloading DFTBephy code from [GitHub](https://github.com/CoMeT4MatSci/dftbephy), you can run the following commands while the `dftbephy` conda environment is active.
```
unzip dftbephy-master.zip
cd dftbephy-master
python setup.py build_ext --inplace
pip install -e .
```


