import sys, os
sys.path.insert(0, '../') # adjust path to the base directory of the package

import numpy as np

import phonopy

from timeit import default_timer as timer

from dftbephy import DftbSuperCellCalc
from dftbephy.fileio import read_dftb_bands, get_lumo
from dftbephy.units import *
from dftbephy.tools import printProgressBar

# this is needed for writing hdf files
import h5py

###############################################################################
# 0 Set base directory in put parameters
basedir = '../examples/Graphene'

nq = 200 # size of q-point mesh (nq x nq x 1)
q_mesh_refinement = 10 # factor by which the q-points are scaled
calculate_veclocities = False

kvec0 = -np.array([2/3,1/3,0])*(1 - 0.0315) # k-point for electrons


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
ph.run_mesh([nq,nq,1], with_eigenvectors=False, is_mesh_symmetry=False, is_gamma_center=False)
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


## very simple filter of ZA and ZO modes
zpol = np.abs(mesh_eigenvectors[:,2,:]) + np.abs(mesh_eigenvectors[:,5,:])
new_frequencies = mesh_frequencies + 100*(zpol>1.1)
srt_idx = np.argsort(new_frequencies, axis=1)
new_frequencies = np.take_along_axis(mesh_frequencies, srt_idx, axis=1)

new_evecs = mesh_eigenvectors
for alpha in range(mesh_eigenvectors.shape[1]):
        new_evecs[:,alpha,:] = np.take_along_axis(mesh_eigenvectors[:,alpha,:], srt_idx, axis=1)
    
mesh_eigenvectors = new_evecs
mesh_frequencies = new_frequencies
##


print('-- starting el-ph calculation ...')
start = timer()
eps_k, mesh_epskq, mesh_epskmq, mesh_g2 = dftb.calculate_g2(kvec0, mesh_qpoints, mesh_frequencies, mesh_eigenvectors)
end = timer()
print('-- finished (%4.1f s).' % (end-start))

# convert q-points to cartesian coordinates
primitive_cell = dftb.primitive.get_cell() * BOHR__AA
a0 = np.linalg.norm(primitive_cell[0,:])

reciprocal_lattice = np.linalg.inv(primitive_cell)
qvecs = a0*np.vstack(mesh_qpoints) @reciprocal_lattice.T # in units of 2*np.pi/a0

# calculate velocities
if calculate_veclocities:
    print('-- starting velocity calculation ...')
    start = timer()
    nqpoints = mesh_g2.shape[0]
    nbands = mesh_g2.shape[2]
    energies = np.zeros((nqpoints, nbands), float)
    velocities = np.zeros((nqpoints, nbands, 3), float)
    for iq, qvec in enumerate(mesh_qpoints):
        eps_k, vel_k = dftb.calculate_velocity(kvec0+qvec)

        for n in range(nbands):
            energies[iq,n] = eps_k[n]
            velocities[iq,n,:] = vel_k[:,n]
    end = timer()
    print('-- finished (%4.1f s).' % (end-start))

# store coupling matrix
nkp = 0
ik = 0
with h5py.File('el-ph-Nq%i-K.hdf5' % (nq), 'w') as f:
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
    
    if calculate_veclocities:
        ds = el_grp.create_dataset('energies', data=energies)
        ds = el_grp.create_dataset('velocities', data=velocities)
