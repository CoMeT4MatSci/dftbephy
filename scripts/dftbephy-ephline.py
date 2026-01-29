import sys, os

import numpy as np

import phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
from phonopy.phonon.dos import TotalDos, Dos

from timeit import default_timer as timer

from dftbephy import DftbSuperCellCalc
from dftbephy.fileio import read_dftb_bands, get_lumo
from dftbephy.analysis import phononDOS, alphaFnk, inv_tau_nk
from dftbephy.units import *
from dftbephy.tools import printProgressBar

import json

# for reading hsd files
import hsd
###############################################################################

def convert(x):
    if hasattr(x, "tolist"):  # numpy arrays have this
        return x.tolist()
    raise TypeError(x)

###############################################################################
# 0 Set base directory and load parameters
def check_hsd_input(inp_dict, name):
    ret_dict = None
    if name in inp_dict.keys():
        ret_dict = inp_dict[name]
    else:
        print('-- ERROR: %s not found in input file.')
    return ret_dict

hsdinput = hsd.load("dftbephy_in.hsd")

inp_dict = check_hsd_input(hsdinput, 'DFTBephy')
assert(type(inp_dict) is dict)

basedir = inp_dict.get('base_dir', './')

phonopy_dir = inp_dict.get('phonopy_dir', 'phonons/')
working_dir = inp_dict.get('working_dir', 'el-ph/')
results_dir = inp_dict.get('results_dir', working_dir)

name = inp_dict.get('name', '')

angular_momenta = inp_dict.get('angularmomenta', {})
if len(angular_momenta) == 0:
    if (rank==0):
        print('-- angular momenta per element have to be specified.')
    quit()

# read section for epc calculations
epc_dict = check_hsd_input(inp_dict, 'EPCs')

default_mesh = {'Mesh': {'npoints': [1, 1, 1], 'refinement': 1, 'shift': [0.0,0.0,0.0]}}
qp_dict = epc_dict.get('qpoints', default_mesh)
if 'Mesh' in qp_dict.keys():
    # size of q-point mesh (eg. nq x nq x 1)
    q_mesh = qp_dict['Mesh'].get('npoints', default_mesh['Mesh']['npoints'])
    # factor by which the q-points are scaled
    q_mesh_refinement = qp_dict['Mesh'].get('refinement', default_mesh['Mesh']['refinement'])
    q_mesh_shift = qp_dict['Mesh'].get('shift', default_mesh['Mesh']['shift'])

k_mesh_shift = epc_dict.get('kvec0', [0., 0., 0.])
kvec0 = np.array(k_mesh_shift) # reference k-point for electrons
#kvec0 = -np.array([2/3,1/3,0])*(1 - 0.0315) # k-point for electrons


# read section for band-path for ephline
path_nodes = epc_dict.get('path', [])
labels = epc_dict.get('labels', [])
npoints = epc_dict.get('npoints', 21)  #21 is the default


band_sel = epc_dict.get('bands', None)
if type(band_sel) is list:
    band_sel[1] = band_sel[1]+1

calculate_velocities = epc_dict.get('velocities', False)

# read section for Phonopy
phonopy_dict = inp_dict.get('Phonopy', {})
phonopy_yaml    = phonopy_dict.get('yaml_file', 'phonopy_disp.yaml')
phonopy_symprec = phonopy_dict.get('symprec', None)

###############################################################################
# 0 Set base directory and band path

# "read section for band-path for ephline" above  does this
###############################################################################
# 1 Load phonopy
print('-- looking for phonopy results in %s' % (basedir + phonopy_dir))
os.chdir(basedir + phonopy_dir)

if phonopy_symprec is None:
    ph = phonopy.load(phonopy_yaml)
    print("-- using Phonopy's default symmetry precision.")
else:
    ph = phonopy.load(phonopy_yaml, symprec=phonopy_symprec)
    print("-- symmetry precision set to", phonopy_symprec)

###############################################################################
# 2 Initialize DFTB calculator
print('-- working in %s' % (basedir + working_dir))
os.chdir(basedir + working_dir)

dftb = DftbSuperCellCalc(angular_momenta)
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

# 3 Run band-structure calculation -- not required

###############################################################################
# 4 Run calculations of electron-phonon couplings along band-path
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
###############################################################################
# to get the path as in the previous version
path = [
    [path_nodes[i], path_nodes[i + 1]]
    for i in range(len(path_nodes) - 1)
]
path_labels = [
    [labels[i], labels[i + 1]]
    for i in range(len(labels) - 1)
]
###############################################################################


qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=npoints)
npaths = len(qpoints)   # number of paths

print('-- starting phonon calculations on path ...')
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


bands = []
g_kq = []
for i in range(npaths):
    print('-- processing path %s to %s' % (path[i][0], path[i][1]))
    nqpoints = qpoints[i].shape[0]

    eps_k, epskq, mesh_epskmq, g2 = dftb.calculate_g2(kvec0, qpoints[i], frequencies[i], eigenvectors[i])

    bands.append(epskq)
    g_kq.append(np.sqrt(g2))


# convert q-points to cartesian coordinates
primitive_cell = dftb.primitive.get_cell() * BOHR__AA
a0 = np.linalg.norm(primitive_cell[0,:])

reciprocal_lattice = np.linalg.inv(primitive_cell)
qvecs = a0*np.vstack(qpoints) @reciprocal_lattice.T # in units of 2*np.pi/a0


# generate output as json
ephl_dict = { 'numBands': bands[0].shape[1], 'energies': np.vstack(bands), 'energyUnit': 'eV',
            'numModes': frequencies[0].shape[1], 'frequencies': np.vstack(frequencies)*THZ__EV, 'frequencyUnit': 'eV',
             'eigenvectors_real': np.real(np.vstack(eigenvectors)),
             'eigenvectors_imaginary': np.imag(np.vstack(eigenvectors)),
            'epcs': np.vstack(g_kq), 'epcUnit': 'eV',
            'coordsType': 'lattice', 'highSymCoordinates': np.vstack(path), 'highSymLabels': np.array(path_labels).flatten(),
            'highSymIndices': np.cumsum(np.array([[0, qps.shape[0]] for qps in qpoints]).flatten()),
            'wavevectorCoordinates': np.vstack(qpoints),
            'wavevectorCoordinatesCart': qvecs,
            'wavevectorIndices': list(range(np.vstack(qpoints).shape[0]))}

with open('ephline.json', 'w') as outfile:
    json.dump(json.dumps(ephl_dict, default=convert), outfile)
