Each electron-phonon coupling calculation starts with two preparatory steps. Be sure to adjust the path to your SlaKos in each dftb_in.hsd.

# Geometry optimization
This step is just the usual (tight) geometry optimization of the respective unit cell. 

# Phonon calculations
In the next step phonopy is used to obtain a supercell, the force constants and eventually the dynamical matrix. For this example, running the follwing code in the phonons directory
```
cp ../geo_end.gen geo.gen
phonopy -d --dim="5 5 1" --amplitude=0.0005 --tolerance=1e-4 --dftb+

# the last step should give files geo.genS, geo.genS-001 and phonopy_disp.yaml
# now we run dftb+ to get the forces
./make_dirs.sh

# collect forces and create force-constants
phonopy -f disp-0*/results.tag --dftb+
```
will give the neccesary input for dftBephy.