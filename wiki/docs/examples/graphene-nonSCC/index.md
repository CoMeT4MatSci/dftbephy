# Graphene

We consider graphene as an example material to demonstrate our approach and provide the inputs for the calculations in [examples/Graphene][dftbephy-examples] directory.
Inputs for geometry optimization and phonon calculations can be found in */opt* and */phonons*, respectively.
The working directory is */el-ph*, in which electron-phonon coupling calculations will be performed. It should contain a dftb_in.hsd and unit cell geo.gen.

In this example, we used the [matsci-0-3](https://github.com/dftbparams/matsci/releases) parameterization set without self-consistent charges (SCC). We additionally tested other parameter sets and the impact of including SCC; however, these variations did not lead to significant changes in the results. For the detailed information about the results for graphene, please see our [paper](https://link.springer.com/article/10.1007/s10825-023-02033-9).


## Geometry Optimization

The structure must be relaxed before starting DFTBephy calculations. For comprehensive guidance and instructions regarding DFTB+ calculations for graphene, one can refer to the [DFTB+ recipes][DFTB+recipes-graphene] page. The input file `dftb_in.hsd`:

```
Geometry = GenFormat {
  <<< graphene.gen
}

Driver = ConjugateGradient {
    LatticeOpt = Yes
    Isotropic = Yes
    MaxForceComponent = 1e-5
}


Hamiltonian = DFTB {
  SCC = No
  MaxAngularMomentum = {
    C = "p"
  }
  Filling = Fermi {
    Temperature [Kelvin] = 100
  }
  SlaterKosterFiles = Type2FileNames {
    Prefix = "/home/service/slako/matsci/matsci-0-3/"
    Separator = "-"
    Suffix = ".skf"
  }
  KPointsAndWeights = SuperCellFolding {
    48 0 0
    0 48 0
    0 0 1
    0.5 0.5 0.0
  }
}
```


The geometry file `graphene.gen`:

```
 2  F
 C
    1  1    0.3333333333E+00    0.3333333333E+00    0.5000000000E+00
    2  1    0.0000000000E+00    0.0000000000E+00    0.5000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.2466605027E+01    0.0000000000E+00    0.0000000000E+00
    0.1233302513E+01    0.2136142614E+01    0.0000000000E+00
    0.0000000000E+00    0.0000000000E+00    0.1399730088E+02
```


## Phonon Calculations

The optimized geometry (`/opt/geo_end.gen`) is copied as `geo.gen` in to the `/phonons` folder and the displacements can then be created via using the following command:

```
phonopy -d --dim="7 7 1" --amplitude=0.005 --tolerance=1e-05 --dftb+
```

Info will be written `phonopy_disp.yaml` file. For each displacement (which is one for graphene-case), calculate the forces.
DFTB+ input for phonon calculations is:

```
Geometry = GenFormat {
  <<< geo.gen
}

Hamiltonian = DFTB {
  SCC = No
  MaxAngularMomentum = {
    C = "p"
  }
  Filling = Fermi {
    Temperature [Kelvin] = 100
  }
  SlaterKosterFiles = Type2FileNames {
    Prefix = "/home/service/slako/matsci/matsci-0-3/"
    Separator = "-"
    Suffix = ".skf"
  }
  KPointsAndWeights = SuperCellFolding {
    12 0 0
    0 12 0
    0 0 1
    0.5 0.5 0.0
  }
}

Options {
   WriteResultsTag = Yes
   WriteAutotestTag = Yes
   WriteHS = No
}

Analysis = {
  PrintForces = Yes
}
```

We also provided `make_dirs.sh` bash script to run the calculations for (each) displacement.

```
./make_dirs.sh
```

Once all displacement calculations have finished, collect the resulting forces and generate the force constants:

```
phonopy -f disp-0*/results.tag --dftb+
```

If it's successfull one should see `"FORCE_SETS" has been created.`.


## Preparing `/el-ph` directory

One must has DFTB+ input and optimized geometry (as a unit cell) in `/el-ph` directory.
Copy the optimized unit cell GEN file used for creating the displacements in `/phonons` directory.
DFTB+ input `dftb_in.hsd` is:

```
Geometry = GenFormat {
  <<< geo.gen
}

Hamiltonian = DFTB {
  SCC = No
  MaxAngularMomentum = {
    C = "p"
  }
  Filling = Fermi {
    Temperature [Kelvin] = 100
  }
  SlaterKosterFiles = Type2FileNames {
    Prefix = "/home/service/slako/matsci/matsci-0-3/"
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
   WriteHS = Yes
}

Analysis = {
  PrintForces = Yes
}
```

In this example we will use `dftbephy-*.py` scripts and therefore we need DFTBephy input file. A minimal input template for graphene is:


```
DFTBephy {
    base_dir = /home/service/workspace/graphene/
    phonopy_dir = phonons/
    working_dir = el-ph/

    name = graphene

    angularmomenta = {
        C = {s p}
    }

    Bands {    }

    EPCs {    }

    RelaxationTimes = SERTA {   }

    Conductivities = SERTA {    }

}
```

A detailed description of each block is provided in the corresponding Graphene exercise sections.


[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job.- http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [DFTB+recipes-graphene]: <https://dftbplus-recipes.readthedocs.io/en/latest/defect/carbon2d-elec.html#perfect-graphene>
   [dftbephy-examples]: <https://github.com/CoMeT4MatSci/dftbephy/tree/master/examples/Graphene>
   [phonopyDFTB+]: <http://phonopy.github.io/phonopy/dftb%2B.html#dftbp-interface>
   [DFTB+recipes]: <https://dftbplus-recipes.readthedocs.io/en/latest/index.html>
   [DFTB+recipes-phonopy]: <https://dftbplus-recipes.readthedocs.io/en/latest/properties/phonopy.html>


