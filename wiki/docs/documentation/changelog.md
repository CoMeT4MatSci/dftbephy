# Change Log

## **[v0.2] - April 2026**

### Code
- Unified scripts for serial calculations
- Unified scripts for parallel calculations
- Added `rtline` command for relaxation times along a path
- Fixed issue with spglib and irreducible k-points in `mobility-mpi` script

### Input file (`dftbephy_in.hsd`)
- Added `DFTB` block
- Updated `EPC` block to use `Path {}` for `ephline`
- Updated `RelaxationTimes` block to use `Path {}` for `rtline`

### Documentation
- Updated *About* section
- Updated *DFTBephy Usage*
- Updated examples: *Bands*, *EPCs*
- Updated installation instructions (phonopy v3 and `spglib < 2.7`)
- Updated *Input File* section
- Added *Change Log*

### Notebooks
- Removed *Graphene-BT2*

### Other
- Updated `pyproject.toml`
- Updated `README.md`
- Removed legacy scripts


## **[v0.1] - March 2026**

### Input file (`dftbephy_in.hsd`)
- Added `phonopy` block
- Added band path definition

### Documentation
- Added graphene example to tutorial

### Code
- Cleaned and updated scripts on GitHub
- Refactored phonopy class (`get` method)
