# Installation

## Prerequisites

Before and during the DFTBephy workflow:
- DFTB+ ( for detailed installation instructions, please visit [DFTB+ Recipes](https://dftbplus-recipes.readthedocs.io/en/latest/introduction.html).)
- [Phonopy](https://phonopy.github.io/phonopy/install.html) (`conda install -c conda-forge "phonopy>=3.0.0" `)

Before install/build DFTBephy:
- **Build-system dependencies**: `numpy`, `Cython`, `setuptools`
- **Other required runtime dependencies**: `scipy`, `spglib<2.7` , `h5py`, `hsd`
- (Optional) Dependencies for MPI version: `mpi4py`, `openmpi`

Build-system dependencies are needed while building DFTBephy itself. Therefore, they must be installed prior to the installation step.

Runtime dependencies are needed when running calculations or scripts after installation. These dependencies are specified in the `pyproject.toml` file and do not need to be installed manually by the user, as they are handled by the build system.
_The main runtime dependencies are:_
`scipy` — scientific computing routines
`spglib<2.7` — symmetry analysis of crystal structures
`h5py` — reading and writing HDF5 files
`hsd` — parsing DFTBephy input files [HSD](https://github.com/dftbplus/hsd-python)
`cython` — also required at runtime for faster routines end extensions

All DFTBephy dependencies are available via both `pip` and `conda`. External tools, in particular DFTB+, may need to be installed separately depending on the workflow.

> **Recommended installation method:**<br>
While not strictly required, we recommend using a conda environment to keep dependencies isolated, avoid conflicts with system-wide packages, and make the setup reproducible. <br>
Using `conda` or `mamba` is often the smoothest way to install packages because it provides precompiled binaries and handles compiled dependencies consistently. In this tutorial we preferred `mamba` since it resolves and installs dependencies faster than `conda`. Use `pip` only for packages not available on `conda-forge`.

## How to install DFTBephy

- Pull / download latest version from [GitHub](https://github.com/CoMeT4MatSci/dftbephy).
- Install build-system dependencies (`pip install numpy cython setuptools`).
- Run `python setup.py build_ext --inplace` in the terminal to build faster routines using Cython.
- Run `pip install -e .` in the terminal to install package but keep it editable in the current directory (**Serial Version**).
- (Optional) Run `python -m pip install -e ".[openmpi]"` to install the **MPI version**.

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

When `dftbephy` environment is still active, install all the prerequisites for DFTBephy in a single step:

```
mamba install numpy cython setuptools "spglib<2.7" "phonopy>=3.0.0" scipy h5py hsd-python openmpi mpi4py
```
Unfortunately, not all versions of phonopy and spglib are compatible with the current version of DFTBephy. Please install **spglib<2.7** and **phonopy>=3.0.0**.


After downloading DFTBephy code from [GitHub](https://github.com/CoMeT4MatSci/dftbephy), you can run the following commands while the `dftbephy` conda environment is active.
```
unzip dftbephy-master.zip
cd dftbephy-master
python setup.py build_ext --inplace
pip install -e . --no-deps
```
Since all dependencies have already been installed via conda/mamba, here, we installed DFTBephy without re-installing dependencies (`--no-deps`).

