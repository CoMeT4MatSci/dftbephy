# Input for DFTBephy

Input is in Human-friendly Structured Data (HSD) format (dftbephy_in.hsd). Detailed information about HSD format can be found on [Read the Docs][HSDReadthedocs] or on [github][HSDgithub].

Before running the calculations, the path to the base directory of DFTBephy package must be added to the [python scripts][scripts]:

```
sys.path.insert(0, '/home/service/dftbephy-master/') # adjust path to the base directory of DFTBephy package.
```

## **DFTBephy{}** 

All parameters are defined under a single section *DFTBephy{}*. This section must contain the following definitions:

**base_dir** Adjusts path to the base directory where python scripts and dftbephy_in.hsd are located. Directories of phonon and electron-phonon coupling calculations are located here as well.

**phonopy_dir** Adjusts path to the directory of phonon calculations.

**working_dir** Adjusts path to the directory of electron-phonon coupling calculations.

**name** Specifies the name of the system.

**angularmomenta** Specifies the maximum angular momentum for each type of atom in the system.

**EPCs** Includes the parameters for electron-phonon coupling (EPC) calculations. Output is written in el-ph-Nq-K-bandsel.hdf5 file. 

**RelaxationTimes** Includes the parameters for relaxation time calculations. Output is written in relaxation-times-fine-mpi-bandsel.hdf5 file.

**Conductivities** Includes the parameters for conductivity calculations. Output is written in transport-mu-mpi.json file.

*Example for input:*
```
DFTBephy = {
    base_dir = /home/user/Templates/Graphyne/
    phonopy_dir = phonons/
    working_dir = el-ph/
    
    name = graphyne
    
    angularmomenta = {
        C = {s p}
    }
    EPCs = {}    
    RelaxationTimes = SERTA {}
    Conductivities = SERTA {}
}
```

Detailed information about the parameters to be defined in *RelaxationTimes*, *Conductivities* and *EPCs* blocks are given in the following subsections.


### **EPCs{}** 
This section includes the parameters for electron-phonon coupling calculations:

| Name | Definition  | Default |
| ------ | ------ | ------ |
| qpoints | *q*-point mesh for Brillouin-zone integration|   |
| npoints | size of *q*-point mesh | [1, 1, 1]  |
| refinement | factor by which the *q*-points are scaled (optional) | 1 |
| kvec0 | reference *k*-point for the electronic states | [0., 0., 0.] |
| bands | selection of electronic bands (optional) | all |
| velocities | group velocity calculation of charge carriers (optional) | no |


*Example for EPCs block:*
```
    EPCs = {
        qpoints = Mesh {
            npoints = 150 150 1
            refinement = 2
        }
        kvec0 = 0.0 0.5 0.0
        bands = 21 26
        velocities = yes
    }
```    

### **RelaxationTimes{}** 
Scattering rates are calculated within the self-energy relaxation time approximation (SERTA). *k*-point mesh must be defined in this section. Band selection is optional for both methods. Additional properties are as follows:

| Name | Definition  | Default | 
| ------ | ------ | ------ |
| kpoints | *k*-point mesh for Brillouin-zone integration|   | 
| npoints | size of *k*-point mesh | [1, 1, 1]  |
| refinement | factor by which the *k*-points are scaled | 1 |
| shift | shift for reference *k*-point for the electronic states | [0., 0., 0.] 
| Efermi [eV] | Fermi level in eV | 0.0 |
| mu [eV] | chemical potential relative to Efermi in eV  |  |
| temperature [eV] | kT energy in eV | 0.0259 |
| sigma [eV] | width of Gaussian smearing function in eV | 0.003 |


*Example for RelaxationTimes block:*
```
    RelaxationTimes = SERTA {
        qpoints = Mesh {
            npoints = 100 100 1
            refinement = 2
        }
        kpoints = Mesh {
            npoints = 16 16 1
            refinement = 4
            shift = 0.0 0.5 0.0
        }
        Efermi [eV] = -4.805057
        bands = 21 26         # optional
        mu [eV] = { -0.78890955 -0.68890955 0.0 0.68890955 0.78890955 }  # either range, list or float
        temperature [eV] = 0.0259
        sigma [eV] = 0.003        
    }
```

### **Conductivities{}** 
Conductivities can be calculated within the constant relaxation time approximation (CRTA) or SERTA. The constant relaxation time value is 1 ps. The parameters in the *Relaxationtimes* section are defined in the same way in this section. *q*-points block is only necessary for SERTA{} method. Additional parameters are as follows:

| Name | Definition  | Default |
| ------ | ------ | ------ |
| Ecut [eV] | a certain energy range for transport properties in eV |  1.0 |
| SpinDegeneracy | degeneracy of the corresponding energy level | 1 |

*Example for Conductivities block:*
```
   Conductivities = CRTA {
        kpoints = Mesh {
            npoints = 400 400 1
            shift = 0.0 0.5 0.0
        }
# only for SERTA
#        qpoints = Mesh {
#            npoints = 100 100 1
#            refinement = 2
#        }
        Efermi [eV] = -4.805057
        bands = 21 26         # optional
        # mu is relativ to Efermi        
        mu [eV] = Range { 0.58890955 0.88890955 20 }
        temperature [eV] = 0.0259
        sigma [eV] = 0.003
        # Ecut is relativ to Efermi 
        Ecut [eV] = 1.0
        SpinDegeneracy = 2
    }
```


[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job.- http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [HSDgithub]: <https://github.com/dftbplus/hsd-python>
   [HSDReadthedocs]: <https://hsd-python.readthedocs.io/en/latest/>
   [mobility-mpi]:<https://github.com/CoMeT4MatSci/dftbephy/blob/master/scripts/dftbephy-mobility-mpi.py>
   [scripts]: <https://github.com/CoMeT4MatSci/dftbephy/tree/master/scripts>
