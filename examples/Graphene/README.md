See the guide on our wiki for how to set up and run calculations: [Graphene (non-SCC) example ](https://comet4matsci.github.io/dftbephy/examples/graphene-nonSCC/notebooks/bands.html)

Before running DFTBephy, first perform the geometry optimization and phonon calculations to generate the required input files.

Be sure to adjust the path to your later–Koster files (SlaKos) in each dftb_in.hsd.

## Geometry optimization
This step is just the usual (tight) geometry optimization of the respective unit cell. 

## Phonon calculations
Phonopy is used to obtain a supercell, the force constants and eventually the dynamical matrix. For this example, running the follwing code in the phonons directory
```
cp ../opt/geo_end.gen geo.gen
phonopy -d --dim="7 7 1" --amplitude=0.005 --tolerance=1e-05 --dftb+

# the last step should give files geo.genS, geo.genS-001 and phonopy_disp.yaml
# now we run dftb+ to get the forces
./make_dirs.sh

# collect forces and create force-constants
phonopy -f disp-0*/results.tag --dftb+
```
will give the neccesary input for DFTBephy.

## DFTBephy calculations

Be sure to adjust the path to the base directory (`base_dir`) in dftbephy_in.hsd.

This workflow supports:
- **Electronic band structures**
- **Phonon band structures**
- **Electron–phonon coupling (EPC)**
- **Scattering rates**
- **Electrical conductivity**
