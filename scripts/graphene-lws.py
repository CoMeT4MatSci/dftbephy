import sys, os
sys.path.insert(0, '../') # adjust path to the base directory of the package

import numpy as np
import matplotlib.pyplot as plt

import phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections

from timeit import default_timer as timer

from dftbephy import DftbSuperCellCalc
from dftbephy.fileio import read_dftb_bands, get_lumo
from dftbephy.analysis import inv_tau_nk
from dftbephy.units import *
from dftbephy.tools import printProgressBar

# this is needed for writing hdf files
#import h5py

# this is needed for writing json files
import json

def convert(x):
    if hasattr(x, "tolist"):  # numpy arrays have this
        return x.tolist()
    raise TypeError(x)

###############################################################################
# 0 Set base directory in put parameters
basedir = '../examples/Graphene'

nq = 100 # size of q-point mesh (nq x nq x 1)

# other parameters below under #5

###############################################################################
# 1 Load phonopy
os.chdir(basedir + '/phonons/')
ph = phonopy.load('phonopy_disp.yaml')

###############################################################################
# 2 Initialize DFTB calculator
os.chdir(basedir + '/el-ph/')

dftb = DftbSuperCellCalc({'C': ['s', 'p']})
dftb.load_phonopy(ph)

if os.path.isfile('reference.npz'):
    print('-- loading reference calculation')
    npzfile = np.load('reference.npz')
    dftb.H0 = npzfile['H0']
    dftb.S0 = npzfile['S0']
else:
    print('-- starting reference calculation ...')
    start = timer()
    dftb.calculate_reference()
    end = timer()
    print('-- finished (%4.1f s).' % (end-start))

    np.savez('reference.npz', H0=dftb.H0, S0=dftb.S0)

###############################################################################
# 3 Run band-structure calculation (optional)

###############################################################################
# 4 Run calculations of electron-phonon couplings

if os.path.isfile('derivatives.npz'):
    print('-- loading derivatives calculation')
    npzfile = np.load('derivatives.npz')
    dftb.H_derivs = npzfile['H_derivs']
    dftb.S_derivs = npzfile['S_derivs']
else:
    print('-- starting derivatives calculation ...')
    start = timer()
    dftb.calculate_derivatives()
    end = timer()
    print('-- finished (%4.1f s).' % (end-start))

    np.savez('derivatives.npz', H_derivs=dftb.H_derivs, S_derivs=dftb.S_derivs)

print('-- starting phonon calculations on mesh ...')
start = timer()
ph.run_mesh([nq,nq,1], with_eigenvectors=True, is_mesh_symmetry=False)
mesh = ph.get_mesh_dict()
end = timer()
print('-- finished (%4.1f s).' % (end-start))

mesh_qpoints= mesh['qpoints']
mesh_frequencies = mesh['frequencies'] 
mesh_eigenvectors = mesh['eigenvectors']

kvec0 = np.array([1/3, 2/3, 0])
print('-- starting el-ph calculation ...')
start = timer()
eps_k, mesh_epskq, mesh_g2 = dftb.calculate_g2(kvec0, mesh_qpoints, mesh_frequencies, mesh_eigenvectors)
end = timer()
print('-- finished (%4.1f s).' % (end-start))


###############################################################################
# 5 Calculate electronic linewidths along band-path
kBT = 0.0259  # temperature in eV
mu = (eps_k[4] + eps_k[3])/2  + 0.1  # chemical potential in eV
sigma_0 = 0.005 # smearing in eV

# path for electron k-points
path=[[[0, 0, 0], [0, 1/2, 0]], [[0, 1/2, 0], [1/3, 2/3, 0]], [[1/3, 2/3, 0], [0, 0, 0]]]
#path=[[[0, 0, 0], [0, 1/2, 0]], [[0, 1/2, 0], [1/3, 1/3, 0]], [[1/3, 1/3, 0], [0, 0, 0]]]
path_labels = [['G', 'M'], ['M', 'K'], ['K', 'G']]



kpoints, connections = get_band_qpoints_and_path_connections(path, npoints=51)

npaths = len(kpoints)
nbands = mesh_g2.shape[2]

bands = []
linewidths = []

for i in range(npaths):
    print('-- processing path %s to %s' % (path[i][0], path[i][1]))
    nkpoints = kpoints[i].shape[0]
    energies = np.zeros((nkpoints, nbands), float)
    taus = np.zeros((nkpoints, nbands), float)

    printProgressBar(0, nkpoints, prefix='k-point', suffix='complete')
    for ik, kvec in enumerate(kpoints[i]):
        eps_k, mesh_epskq, mesh_g2 = dftb.calculate_g2(kvec, mesh_qpoints, mesh_frequencies, mesh_eigenvectors)
        
        for n in range(nbands):
            energies[ik,n] = eps_k[n]
            taus[ik,n] = inv_tau_nk( n, eps_k[n], mu, kBT, mesh_g2, mesh_epskq, mesh_frequencies*THZ__EV, sigma=sigma_0)[0]

        printProgressBar(ik+1, nkpoints, prefix='k-point', suffix='complete')

    bands.append(energies)
    linewidths.append(taus)

# generate output as json
bs_dict = { 'particleType': 'electron', 'numBands': bands[0].shape[1], 'energies': np.vstack(bands), 'energyUnit': 'eV', 
  'coordsType': 'lattice', 'highSymCoordinates': np.vstack(path), 'highSymLabels': np.array(path_labels).flatten(), 
  'highSymIndices': np.cumsum(np.array([[0, kps.shape[0]] for kps in kpoints]).flatten()),
  'wavevectorCoordinates': np.vstack(kpoints), 'wavevectorIndices': list(range(np.vstack(kpoints).shape[0]))}

with open('path_el_bandstructure.json', 'w') as outfile:
    json.dump(json.dumps(bs_dict, default=convert), outfile)

# generate output as json
lw_dict = { 'chemicalPotentialUnit': 'eV', 'chemicalPotentials': [mu,], 'particleType': 'electron', 'numBands': bands[0].shape[1], 
  'energies': np.vstack(bands), 'energyUnit': 'eV', 
  'coordsType': 'lattice', 'linewidths': np.vstack(linewidths), 'linewidthsUnit': 'eV', 'temperatureUnit': 'eV', 'temperatures': [kBT,]}

with open('path_el_relaxation_times.json', 'w') as outfile:
    json.dump(json.dumps(lw_dict, default=convert), outfile)