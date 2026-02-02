# Input File

The main input file, `dftbephy_in.hsd`, is written in Human-friendly Structured Data (HSD) format.
Detailed information about the HSD syntax can be found on [Read the Docs][HSDReadthedocs] or on [GitHub][HSDgithub].


## DFTBephy {}

All parameters are defined under a single section *DFTBephy={}*. This section must contain the following definitions:

**base_dir** Adjusts path to the base directory where python scripts and dftbephy_in.hsd are located. Directories of phonon and electron-phonon coupling calculations are located here as well.

**phonopy_dir** Adjusts path to the directory of phonon calculations.

**working_dir** Adjusts path to the directory of electron-phonon coupling calculations.

**name** Specifies the name of the system.

**angularmomenta** Specifies the maximum angular momentum for each type of atom in the system.

**Phonopy** This is an optional block to define settings related to the phonopy calculation. It is not required if the default phonopy parameters are used. Specify this block only when you rename the `yaml` file or change the symmetry tolerance `symprec`.


Add the following blocks depending on the calculation to be performed:

**Bands** Includes the parameters for band structure calculations.

**EPCs** Includes the parameters for electron-phonon coupling (EPC) calculations. Output is written in HDF file.

**RelaxationTimes** Includes the parameters for relaxation time calculations. Output is written in HDF file. 

**Conductivities** Includes the parameters for conductivity calculations. Output is written in JSON file. 

The name of the output file for each calculation can be changed in the Python [scripts](https://github.com/CoMeT4MatSci/dftbephy/tree/master/scripts).


*Example for input:*
```
DFTBephy {
    base_dir = /home/user/Templates/Graphene/
    phonopy_dir = phonons/
    working_dir = el-ph/

    name = graphene

    angularmomenta = {
        C = {s p}
    }

   Phonopy = {
        yaml_file = phonopy_disp.yaml  #default
        symprec   = 1e-5               #default
    }

    Bands = {}
    EPCs = {}
    RelaxationTimes = SERTA {}
    Conductivities = SERTA {}
}
```

Detailed information about the parameters to be defined in *Bands {}*, *EPCs {}*, *RelaxationTimes {}* and *Conductivities {}* blocks are given in the following subsections.

### Bands {}


| Name    | Definition                                                                 | Default |
|-------- | -------------------------------------------------------------------------- | ------- |
| path    | List of high symmetry points in fractional coordinates |         |
| labels  | Labels for the high symmetry points in `path`  |         |
| npoints | Number of *k*-points between each consecutive high symmetry points in `path`|   21    |


*Example for Bands block:*
```
   Bands {
        path = {
             0.0       0.0       0.0
             0.0       0.5       0.0
             0.333333  0.666667  0.0
             0.0       0.0       0.0
        }

        labels = { G M K G }
        npoints = 50

    }
```



### EPCs {}
This section includes the parameters for electron-phonon coupling calculations:

| Name | Definition  | Default |
| ------ | ------ | ------ |
| qpoints | *q*-point mesh for Brillouin-zone integration|   |
| npoints | size of *q*-point mesh | [1, 1, 1]  |
| refinement | factor by which the *q*-points are scaled (optional) | 1 |
| kvec0 | reference *k*-point for the electronic states | [0., 0., 0.] |
| bands | selection of electronic bands (optional) | all |
| velocities | group velocity calculation of charge carriers (optional) | no |

The ```bands``` option only controls which bands are saved in the output. The calculations are always performed including all bands. It can be defined in `EPCs {}`, `RelaxationTimes {}` and `Conductivities {}` blocks.

For EPC calculations along a path (ephline), the ```path```, ```labels```, and ```npoints``` parameters are defined in the same way as in the ```Bands{}``` block.

*Example for EPCs block:*
```
    EPCs {
        kvec0 = 0.32283333 0.64566667 0.0

        ##### only for EPC calc. on mesh ######
        qpoints = Mesh {
            npoints = 200 200 1
            refinement = 1
        }
        bands = 3 4
        velocities = yes


        ########################################
        #### only for EPC calc. along path ####
        path = {
             0.0       0.0       0.0
             0.0       0.5       0.0
             0.333333  0.666667  0.0
             0.0       0.0       0.0
        }
        labels = { G M K G }
        npoints = 51
        ########################################
    }
```

### RelaxationTimes = {}
Relaxation times are calculated within self-energy relaxation time approximation (SERTA). Additional properties for relaxation time calculations are as follows:

| Name | Definition  | Default |
| ------ | ------ | ------ |
| kpoints | *k*-point mesh for Brillouin-zone integration|   |
| npoints | size of *k*-point mesh | [1, 1, 1]  |
| refinement | factor by which the *k*-points are scaled | 1 |
| shift | shift for reference *k*-point for the electronic states | [0., 0., 0.]
| Efermi [eV] | Fermi energy in eV | 0.0 |
| mu [eV] | chemical potential relative to Efermi in eV  |  |
| temperature [eV] | kT energy in eV | 0.0259 |
| sigma [eV] | width of Gaussian smearing function in eV | 0.003 |


*Example for RelaxationTimes block:*
```
    RelaxationTimes = SERTA {
        qpoints = Mesh {
            npoints = 200 200 1                   # size of q-point mesh (nq x nq x 1)
            refinement = 10                       # factor by which the q-points are scaled
        }
        kpoints = Mesh {
            npoints = 16 16 1
            refinement = 20                       # k-point scaling value is always larger than q-point scaling value
            shift = 0.32283333 0.64566667 0.0
        }
        Efermi [eV] = -4.6585
        mu [eV] = { 0.1 }                         # as a list
        sigma [eV] = 0.003
        temperature [eV] = 0.0259
   }

```


*k*-point refinement should be larger than *q*-point refinement.

Chemical potential `mu` can be specified either as a range, `Range { min max number_of_values }`, or as an explicit list, `{ value1 value2 value3 }`.


### Conductivities = {}
One of the ```CRTA``` or ```SERTA``` methods must be selected. The parameters in the `Relaxationtimes` block are defined in the same way here.
When the ```SERTA``` method is selected, ```qpoints``` must be specified using a ```Mesh``` block (see the example below).
Additional parameters are as follows:

| Name | Definition  | Default |
| ------ | ------ | ------ |
| Ecut [eV] | a certain energy range for transport properties in eV |  1.0 |
| SpinDegeneracy | degeneracy of the corresponding energy level | 1 |

*Example for Conductivities block:*
```
     Conductivities = SERTA {
        kpoints = Mesh {
            npoints = 400 400 1
            shift = 0.32283333 0.64566667 0.0
            refinement = 20
        }
        ############## only for SERTA ###########
        qpoints = Mesh {
            npoints = 200 200 1
            refinement = 10
        }
        #########################################

        Efermi [eV] = -4.6585
        bands = 3 4                                 # optional
        mu [eV] = Range { 0.0 0.40 40 }             # as a range
        temperature [eV] = 0.0259

        sigma [eV] = 0.003
        Ecut [eV] = 1.0
        SpinDegeneracy = 2
    }
```


[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job.- http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [HSDgithub]: <https://github.com/dftbplus/hsd-python>
   [HSDReadthedocs]: <https://hsd-python.readthedocs.io/en/latest/>
