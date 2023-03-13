import sys, os
sys.path.insert(0, '../') # adjust path to the base directory of the package

import numpy as np

import phonopy

from timeit import default_timer as timer

from dftbephy import DftbSuperCellCalc
from dftbephy.units import *
from dftbephy.tools import printProgressBar

# this is needed for writing hdf files
import h5py

# for reading hsd files
import hsd

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

band_sel = epc_dict.get('bands', None)
if type(band_sel) is list:
    band_sel[1] = band_sel[1]+1

calculate_velocities = epc_dict.get('velocities', False)

###############################################################################
# 1 Load phonopy
print('-- looking for phonopy results in %s' % (basedir + phonopy_dir))
os.chdir(basedir + phonopy_dir)
ph = phonopy.load('phonopy_disp.yaml')

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
ph.run_mesh(q_mesh, with_eigenvectors=False, is_mesh_symmetry=False, is_gamma_center=False)
mesh = ph.get_mesh_dict()
end = timer()
print('-- finished (%4.1f s).' % (end-start))

mesh_qpoints= mesh['qpoints']

# we are only interested in a small region around kvec0
qps = mesh['qpoints']/q_mesh_refinement
ph.run_qpoints(qps, with_eigenvectors=True)

mesh = ph.get_qpoints_dict()
mesh_qpoints = qps
mesh_frequencies = mesh['frequencies'] 
mesh_eigenvectors = mesh['eigenvectors']


print('-- starting el-ph calculation ...')
start = timer()
eps_k, mesh_epskq, mesh_epskmq, mesh_g2 = dftb.calculate_g2(kvec0, mesh_qpoints, mesh_frequencies, mesh_eigenvectors, band_sel=band_sel)
end = timer()
print('-- finished (%4.1f s).' % (end-start))

# convert q-points to cartesian coordinates
primitive_cell = dftb.primitive.get_cell() * BOHR__AA
a0 = np.linalg.norm(primitive_cell[0,:])

reciprocal_lattice = np.linalg.inv(primitive_cell)
qvecs = a0*np.vstack(mesh_qpoints) @reciprocal_lattice.T # in units of 2*np.pi/a0

# calculate velocities
if calculate_velocities:
    print('-- starting velocity calculation ...')
    start = timer()
    nqpoints = mesh_g2.shape[0]
    nbands = mesh_g2.shape[2]
    energies = np.zeros((nqpoints, nbands), float)
    velocities = np.zeros((nqpoints, nbands, 3), float)
    for iq, qvec in enumerate(mesh_qpoints):
        eps_k, vel_k = dftb.calculate_velocity(kvec0+qvec, band_sel=band_sel)

        for n in range(nbands):
            energies[iq,n] = eps_k[n]
            velocities[iq,n,:] = vel_k[:,n]
    end = timer()
    print('-- finished (%4.1f s).' % (end-start))

# store coupling matrix
nkp = 0
ik = 0
with h5py.File('el-ph-Nq%i-K-bandsel.hdf5' % (q_mesh[0]), 'w') as f:
    ph_grp = f.create_group('ph')
    ph_grp.attrs['mesh'] = ph._mesh.get_mesh_numbers()
    ds = ph_grp.create_dataset('omega', data=mesh_frequencies*THZ__EV)
    ds = ph_grp.create_dataset('qpoints', data=mesh_qpoints)
    ds = ph_grp.create_dataset('qpointsCart', data=qvecs)    
    
    el_grp = f.create_group('el')
    elph_grp = f.create_group('el-ph')

    ds = elph_grp.create_dataset('g2_%i' % (nkp+ik), data=mesh_g2)
    ds.attrs['kvec'] = kvec0
    ds = el_grp.create_dataset('eps_%i' % (nkp+ik), data=eps_k)
    ds.attrs['kvec'] = kvec0
    ds = el_grp.create_dataset('eps_q_%i' % (nkp+ik), data=mesh_epskq)
    ds.attrs['kvec'] = kvec0
    ds = el_grp.create_dataset('eps_mq_%i' % (nkp+ik), data=mesh_epskmq)
    ds.attrs['kvec'] = kvec0
    
    if calculate_velocities:
        ds = el_grp.create_dataset('energies', data=energies)
        ds = el_grp.create_dataset('velocities', data=velocities)
