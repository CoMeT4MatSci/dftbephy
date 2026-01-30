# About

DFTBephy is a computational framework based on the density functional tight binding (DFTB) method for evaluating electron–phonon couplings and relaxation times. It interfaces with [phonopy](https://phonopy.github.io/phonopy/reference.html) to obtain vibrational modes and with [DFTB+](https://dftbplus.org/about/index.html) to compute the corresponding transport properties. Our implementation is openly available on [GitHub](https://github.com/CoMeT4MatSci/dftbephy).

## What you can get

The main purpose of DFTBephy is the calculation of electron-phonon couplings. Apart from that, the package also allows the calculation of the electronic band-structure and the charge carrier scattering rates (at the moment only within SERTA).


The [scripts/](https://github.com/CoMeT4MatSci/dftbephy/tree/master/scripts) directory contains some python scripts for computing:
- EPCs at k-point on a (fine) q-mesh and along a given q-path
- Scattering rates on a (fine) k-mesh and along a given k-path
- Conductivity tensor



For more information see our articles:
> Croy, A., Unsal, E., Biele, R., Pecchia A. *DFTBephy: A DFTB-based approach for electron–phonon coupling calculations.* J Comput Electron, 22, 1231 (2023).
> doi:10.1007/s10825-023-02033-9
> [full-text access](https://link.springer.com/article/10.1007/s10825-023-02033-9)

> Unsal, E., Pecchia A., Croy, A., Cuniberti, G. *Charge carrier mobilities in γ-graphynes: a computational approach.* Nanoscale, 17, 24591 (2025).
> doi:10.1039/D5NR02989A
> [full-text access](https://pubs.rsc.org/en/content/articlelanding/2025/nr/d5nr02989a)

The BibTeX entries are available in [`references.bib`](https://github.com/CoMeT4MatSci/dftbephy/blob/master/references.bib).


[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job.- http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

