# dftBephy
Calculating electron-phonon couplings (EPCs) with DFTB.

![dftbephy](logo.png)

# Prerequisites

- numpy and scipy
- [dftb+](https://github.com/dftbplus/dftbplus)
- [phonopy](https://github.com/phonopy/phonopy) (`conda install -c conda-forge phonopy`)
- setuptools (`pip install setuptools`)
- cython (`pip install Cython`)

# Installation

- Pull / download latest version from github.
- Run `python setup.py build_ext --inplace` in the terminal to build faster routines.

# Running calculations

- Starting point for all dftBephy calculations is a finished phonopy calculation of the force constants (e.g. FORCE_SETS and phonopy_disp.yaml).
- The working directory should contain a dftb_in.hsd file, which reads the geometry from geo.gen (will be written by dftBephy) and contains the option `WriteHS = Yes` (to be used by dftBephy). The directory may also contain charges.bin from a previous SCC run.
- See the examples/ directory for more details. (It's recommended to copy one of the examples and adapt it to your needs.)
- A detailed information about DFTBephy input (dftbephy_in.hsd) can be found [here][dftbephyinput].

# What you can get
The main purpose of dftBephy is the calculation of electron-phonon couplings. Apart from that, the package also allows the calculation of the electronic band-structure and the electron relaxation-time (at the moment only within SERTA). The latter can be used as an input for BoltzTrap2 to calculate transport properties. The scripts/ directory contains some templates for computing
- EPCs and relaxation times along a band path (q and k paths, respectively) -- graphene-ephline.py and graphene-lws.py -- the results are stored in a json file;
- EPCs at k-point on a (fine) q-mesh -- graphene-epc.py -- the result is stored in a hdf5 file;
- Relaxation times on a (fine) k-mesh -- graphene-relaxationtimes-fine.py -- the results are stored in a hdf5 file;
- Relaxation times as input for Boltztrap2 -- graphene-mobility-bt2.py -- the results are stored in a hdf5 file. For this script, ase and spglib are required;
- Conductivity tensor -- graphene-mobility.py -- the results are stored in a json file. For this script, ase and spglib are required.

See jupyter notebooks in notebooks/ for how to read, use and visualize the output.

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job.- http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

[dftbephyinput]: <https://github.com/CoMeT4MatSci/dftbephy/blob/master/Input_for_DFTBephy.md>
