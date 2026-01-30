# DFTBephy Usage

In preparation for DFTBephy calculations, it is essential first to carry out tight geometry optimization of the system with [DFTB+](https://dftbplus.org/about/index.html) and then to perform phonon calculations using [phonopy][phonopyDFTB+].


## Workflow

First, a tight DFTB+ geometry optimization is performed, followed by phonon calculations with phonopy. Based on the optimized geometry and the information in `phonopy_disp.yaml`, a supercell is constructed. DFTBephy uses DFTB+ to obtain the Hamiltonian and the overlap matrices in real space and DFTBephy then takes the Fourier transform of them. In the next step, each atom in the central cell is displaced and the real-space gradients of the Hamiltonian and the overlap matrices are found from a finite-difference scheme. Then, the gradients are also Fourier transformed. Now, together with the output of phonopy, i.e. phonon frequencies and polarizations, the electron-phonon coupling (EPC) constants $g^\lambda_{nm}(\vec{k},\vec{q})$ can be computed. From the EPCs, one can then calculate the scattering rates and the electrical conductivities.

:::{figure} workflow.png
:width: 100%
:name: fig-workflow
Workflow for DFTBephy calculations.
:::

At present, only the self-energy relaxation-time approximation (SERTA) is implemented for the scattering rates. For the conductivity calculations, a constant relaxation-time approximation (CRTA) option is also available. Mobilities are calculated as a ratio of conductivity to carrier concentration in the post-process.The detailed information about the formulation can be found in the [paper](https://link.springer.com/article/10.1007/s10825-023-02033-9).

## Preparing the workspace

On GitHub, all input files for graphene are provided in the [/examples/Graphene](https://github.com/CoMeT4MatSci/dftbephy/tree/master/examples/Graphene) directory. It's recommended that you copy and adapt one of the examples to suit your requirements.

Including the DFTBephy input file and the corresponding Python scripts, `/examples/Graphene` can be used directly as workspace. The directory structure is shown schematically below.


```
.
├── opt/                  # not used by DFTBephy
├── phonons/              # Only the files used by DFTBephy are shown
│   ├── phonopy_disp.yaml
│   └── FORCE_SETS
├── el-ph/
│   ├── dftb_in.hsd
│   └── geo.gen           # optimized unit-cell geometry
├── dftbephy_in.hsd       # DFTBephy input file
└── dftbephy-*.py         # Python driver scripts
```

The important parameters to be added to `dftb_in.hsd` in `el-ph/` directory:

```
Geometry = GenFormat {
  <<< geo.gen
}

Hamiltonian = DFTB {
  SCC = No
# ReadInitialCharges = Yes              # the charge file must be provided
  MaxAngularMomentum = {}

  SlaterKosterFiles = Type2FileNames {
    Prefix = "../../../"
    Separator = "-"
    Suffix = ".skf"
  }
  KPointsAndWeights = SuperCellFolding {
    1 0 0
    0 1 0
    0 0 1
    0.5 0.5 0.0
  }
}

Options {
   WriteResultsTag = Yes
   WriteAutotestTag = Yes
   WriteHS = Yes                        # to print out H and S matrices
}

Analysis = {
  PrintForces = Yes
}
```

Hamiltonian and overlap matrices can also be built with previously converged charge distribution (frozen charge method). In this case, the directory should also contain charges.bin, and the option `ReadInitialCharges = Yes` should be added in `dftb_in.hsd` file as shown above. The corresponding charge data file (`charges.bin` or `charges.dat`) therefore has to be provided as well.


A minimal template of `dftbephy_in.hsd` file is:

```
DFTBephy {
  base_dir     =                     # path to workspace (contains phonons/, el-ph/, ...)
  phonopy_dir  =                     # relative to base_dir
  working_dir  =                     # relative to base_dir
  name         =

  angularmomenta = { }               # angular momentum for each element

  Bands { }                          # band settings
  EPCs { }                           # electron-phonon coupling settings
  RelaxationTimes = SERTA { }        # scattering settings (only SERTA)
  Conductivities { }                 # transport settings
}
```

One can find the detailed info about the input file `dftbephy_in.hsd` in {ref}`inputfile-section` section.

## Scripts

We provide scripts for the calculations in [/scripts][dftbephy-scripts] directory.
The Python scripts named `dftbephy-*.py` can be used with the DFTBephy input (`dftbephy_in.hsd`) which allows the certain parameters to be easily modified. Detailed information about `dftbephy_in.hsd` can be found {ref}`inputfile-section` section.


If a conda environment has been created, it should first be activated before running the Python script. For example:
```
conda activate dftbephy
python dftbephy-bands.py
```

The parallel scripts (`dftbephy-*-mpi.py`) distribute the workload over k-points and can be run with MPI, e.g.:

```
mpirun -np 16 python dftbephy-relaxationtimes-mpi.py |& tee -a run-scatterings.log
```

Depending on your system, it can be useful to set the thread limits (e.g. `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, etc.) in the `dftbephy-*-mpi.py` scripts so that the total concurrency (ranks × threads per rank) matches the available CPU cores and avoids oversubscription.


## Post-processing

We also provide Jupyter Notebooks in [/notebooks][dftbephy-notebooks] directory on GitHub. These notebooks can be used for post-processing, like plotting band structures, EPC contour plots, inverse lifetimes, and calculating charge carrier mobilities.


[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job.- http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [phonopyDFTB+]: <http://phonopy.github.io/phonopy/dftb%2B.html#dftbp-interface>

   [dftbephy-examples]: <https://github.com/CoMeT4MatSci/dftbephy/tree/master/examples/Graphene>
   [dftbephy-scripts]: <https://github.com/CoMeT4MatSci/dftbephy/tree/master/scripts>
   [dftbephy-notebooks]: <https://github.com/CoMeT4MatSci/dftbephy/tree/master/notebooks>

   [dftbephy-input]: documentation/inputfile

