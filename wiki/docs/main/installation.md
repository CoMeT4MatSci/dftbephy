# Installation

## Prerequisites

- NumPy and SciPy (available via `pip` and `conda`)
- DFTB+ ( for detailed installation instructions, please visit [DFTB+ Recipes](https://dftbplus-recipes.readthedocs.io/en/latest/introduction.html).)
- [Phonopy](https://phonopy.github.io/phonopy/install.html) (`conda install -c conda-forge phonopy`)
- [HSD](https://github.com/dftbplus/hsd-python) for DFTBephy input (`hsd-python` is available via `conda`)
- Cython for faster routines end extensions (available via `pip` and `conda`)
- Other required runtime dependencies: mpi4py, openmpi, spglib, h5py (available via `pip` and `conda`)


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

Create a new conda environment using a recent Python version and then activate it:

```
conda create --name dftbephy
conda activate dftbephy
```

> **Optional:**
> In this tutorial we use `python=3.14`. You can pin the Python version for the environment as:
> ```bash
> conda create -n dftbephy "python=3.14"
> ```


Optionally, you can install DFTB+ inside this conda environment. If DFTB+ is already available on the machine, you can skip this step.

```
mamba install 'dftbplus=*=nompi_*'
mamba install  dftbplus-tools dftbplus-python
```

When `dftbephy` environment is still active, install the prerequisites for DFTBephy:

```
mamba install numpy scipy cython "spglib<2.7" openmpi mpi4py h5py setuptools hsd-python "phonopy<=2.48"
```

Unfortunately, not all versions of phonopy and spglib are compatible with the current version of DFTBephy. Please install **spglib<2.7** and **phonopy<=2.48.0**.


After downloading DFTBephy code from [GitHub](https://github.com/CoMeT4MatSci/dftbephy), you can run the following commands while the `dftbephy` conda environment is active.
```
unzip dftbephy-master.zip
cd dftbephy-master
python setup.py build_ext --inplace
pip install -e .
```


