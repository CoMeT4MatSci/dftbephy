# dftBephy
Calculating electron-phonon couplings (EPCs) with DFTB.

# Prerequisites

- dftb+
- numpy and scipy
- phonopy (`conda install -c conda-forge phonopy`)
- setuptools (`pip install setuptools`)
- cython (`pip install Cython`)

# Installation

- Pull / download latest version from github.
- Run `python setup.py build_ext --inplace` in the terminal to build faster routines.

# Running calculations

- Starting point for all dftBephy calculations is a finished phonopy calculation of the force constants (e.g. FORCE_SETS and phonopy_disp.yaml).
- The working directory should contain a dftb_in.hsd file, which reads the geometry from geo.gen (will be written by dftBephy) and contains the option `WriteHS = Yes` (to be used by dftBephy). The directory may also contain charges.bin from a previous SCC run.
- See the examples/ directory for more details. (It's recommended to copy one of the examples and adapt it to your needs.)


# What you can get
The main purpose of dftBephy is the calculation of electron-phonon couplings. Apart from that, the package also allows the calculation of the electronic band-structure and the electron relaxation-time (at the moment only within SERTA). The latter can be used as an input for BoltzTrap2 to calculate transport properties. The scripts/ directory contains some templates for computing
- EPCs and relaxation times along a band path (q and k paths, respectively) -- graphene-ephline.py and graphene-lws.py -- the results are stored in as json.
- EPCs at k-point on a fine q-mesh -- graphene_epc-py -- the result is stored in a hdf5 file. 
- Relaxation times and input for Boltztrap2 -- graphene_mobility.py -- the results are stored in a hdf5 file. For this script, ase and spglib are required.
See jupyter notebooks in notebooks/ for how to read, use and visualize the output.