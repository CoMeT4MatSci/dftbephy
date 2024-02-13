import sys, os

import numpy as np
import matplotlib.pyplot as plt

import phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections

from timeit import default_timer as timer

from dftbephy import DftbSuperCellCalc
from dftbephy.units import *
from dftbephy.tools import printProgressBar

# this is needed for writing json files
import json

def convert(x):
    if hasattr(x, "tolist"):  # numpy arrays have this
        return x.tolist()
    raise TypeError(x)

###############################################################################
# 0 Set base directory in put parameters
basedir = '../examples/Graphene'

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
# 3 Run band-structure calculation

# path for electron k-points
path=[[[0, 0, 0], [0, 1/2, 0]], [[0, 1/2, 0], [1/3, 2/3, 0]], [[1/3, 2/3, 0], [0, 0, 0]]]
path_labels = [['G', 'M'], ['M', 'K'], ['K', 'G']]

# get k-points
kpoints, connections = get_band_qpoints_and_path_connections(path, npoints=51)

npaths = len(kpoints)
nbands = int(dftb.H0.shape[0]/(7*7))
print(nbands)

bands = []
vels = []

for i in range(npaths):
    print('-- processing path %s to %s' % (path[i][0], path[i][1]))
    nkpoints = kpoints[i].shape[0]
    energies = np.zeros((nkpoints, nbands), float)
    velocities = np.zeros((nkpoints, nbands, 3), float)

    printProgressBar(0, nkpoints, prefix='k-point', suffix='complete')
    for ik, kvec in enumerate(kpoints[i]):
        eps_k, vel_k = dftb.calculate_velocity(kvec)
                
        for n in range(nbands):
            energies[ik,n] = eps_k[n]
            velocities[ik,n,:] = vel_k[:,n]

        printProgressBar(ik+1, nkpoints, prefix='k-point', suffix='complete')

    bands.append(energies)
    vels.append(velocities)

# convert q-points to cartesian coordinates
primitive_cell = dftb.primitive.get_cell() * BOHR__AA
a0 = np.linalg.norm(primitive_cell[0,:])

reciprocal_lattice = np.linalg.inv(primitive_cell)
kvecs = a0*np.vstack(kpoints) @reciprocal_lattice.T # in units of 2*np.pi/a0
    
# generate output as json
bs_dict = { 'particleType': 'electron', 'numBands': bands[0].shape[1], 'energies': np.vstack(bands), 'energyUnit': 'eV', 
  'coordsType': 'lattice', 'highSymCoordinates': np.vstack(path), 'highSymLabels': np.array(path_labels).flatten(), 
  'highSymIndices': np.cumsum(np.array([[0, kps.shape[0]] for kps in kpoints]).flatten()),
  'wavevectorCoordinates': np.vstack(kpoints), 
  'wavevectorCoordinatesCart': kvecs,
  'wavevectorIndices': list(range(np.vstack(kpoints).shape[0]))}

with open('path_el_bandstructure.json', 'w') as outfile:
    json.dump(json.dumps(bs_dict, default=convert), outfile)

###############################################################################
# 4 Run phonon dispersion

qpoints = kpoints # use the same paths as for electrons
qvecs = kvecs

print('-- starting phonon calculations on paths ...')
start = timer()
ph.run_band_structure(qpoints, is_band_connection=False, path_connections=connections, with_eigenvectors=True)
bs_dict = ph.get_band_structure_dict()
end = timer()
print('-- finished (%4.1f s).' % (end-start))

######## reorder phonon modes ###########
## based on the idea presented in https://quantumtinkerer.tudelft.nl/blog/connecting-the-dots/
##
from scipy.optimize import linear_sum_assignment

ph_evs = bs_dict['eigenvectors']
ph_oms = bs_dict['frequencies']

npaths = len(ph_oms)

ev_ref = ph_evs[0][0]
for lam in range(ev_ref.shape[1]):
    ev_ref[:,lam] = ev_ref[:,lam] * np.exp( -1j * np.angle(ev_ref[0,lam]) ) 

sorted_freqs = []
sorted_evs = []
for pi in range(npaths):

    sorted_freqs_path = []
    sorted_evs_path = []
    for i in range(ph_evs[pi].shape[0]):
        eigvecs = ph_evs[pi][i]
        for lam in range(eigvecs.shape[1]):
            eigvecs[:,lam] = eigvecs[:,lam] * np.exp( -1j * np.angle(eigvecs[0,lam]) ) 
        metric = np.abs(np.dot(ev_ref.conj().T, eigvecs))
        
        assignment = linear_sum_assignment(-metric)[1]
        
        sorted_freqs_path.append(ph_oms[pi][i][assignment])
        sorted_evs_path.append(eigvecs[:, assignment])
        ev_ref = eigvecs[:, assignment]
        
    sorted_freqs.append(np.array(sorted_freqs_path))
    sorted_evs.append(np.array(sorted_evs_path))
######################################### 

qpoints= bs_dict['qpoints']
frequencies = sorted_freqs
eigenvectors = sorted_evs


# generate output as json
bs_dict = { 'particleType': 'phonon', 'numModes': frequencies[0].shape[1], 'frequencies': np.vstack(frequencies)*THZ__EV, 'frequencyUnit': 'eV', 
  'coordsType': 'lattice', 'highSymCoordinates': np.vstack(path), 'highSymLabels': np.array(path_labels).flatten(), 
  'highSymIndices': np.cumsum(np.array([[0, qps.shape[0]] for qps in qpoints]).flatten()),
  'wavevectorCoordinates': np.vstack(qpoints),
  'wavevectorCoordinatesCart': qvecs,
  'wavevectorIndices': list(range(np.vstack(qpoints).shape[0]))}

with open('path_ph_bandstructure.json', 'w') as outfile:
    json.dump(json.dumps(bs_dict, default=convert), outfile)
