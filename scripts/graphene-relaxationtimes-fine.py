import os, sys

import numpy as np

import phonopy
from phonopy.structure.grid_points import GridPoints
from timeit import default_timer as timer

from dftbephy import DftbSuperCellCalc
from dftbephy.fileio import read_dftb_bands, get_lumo
from dftbephy.analysis import inv_tau_nk, inv_tau_nk_lam
from dftbephy.units import *
from dftbephy.tools import printProgressBar


# this is needed for writing hdf files
import h5py


###############################################################################
# 0 Set base directory
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
# 3 Run band-structure calculation (optional)

###############################################################################
# 4 Run calculations of electron-phonon couplings
nq = 200 # size of q-point mesh (nq x nq x 1)
q_mesh_refinement = 10 # factor by which the q-points are scaled

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

kvec0 = np.array([-2/3, -1/3, 0]) # k-point for electrons

print('-- starting el-ph calculation ...')
start = timer()
eps_k, mesh_epskq, mesh_epskmq, mesh_g2 = dftb.calculate_g2(kvec0, mesh_qpoints, mesh_frequencies, mesh_eigenvectors)
end = timer()
print('-- finished (%4.1f s).' % (end-start))


###############################################################################
# 6 Calculate relaxation times
nk = 16 # size of k-point mesh (nk x nk x 1)


#mus = np.linspace(eps_k[3] - 0.20, eps_k[3] + 0.20, 21) # chemical potentials in eV
#kBTs = 0.0259 * np.ones_like(mus)  # temperatures in eV

sigmas = np.array([0.003, 0.006, 0.009]) # smearing in eV
kBTs = 0.0259 * np.ones_like(sigmas)   # temperatures in eV
mus = (eps_k[4] + 0.10)* np.ones_like(kBTs) # chemical potentials in eV

EF = (eps_k[3]+eps_k[4])/2    # Fermi energy in eV


temps_chems = zip(kBTs, mus)
ncalcs = len(list(temps_chems))

print('-- running %d calculations for mu, kBT, sigma:' % (ncalcs))
for ic, (kBT, mu, sigma_0) in enumerate(zip(kBTs, mus, sigmas)):
    print(ic, mu, kBT, sigma_0)


gp = GridPoints(np.array([nk,nk,1]), np.linalg.inv(dftb.primitive.cell), 
                 rotations=None, is_mesh_symmetry=False, is_gamma_center=False, is_time_reversal=False)
kpoints = gp.qpoints / 20
kpoints -= np.array([2/3,1/3,0] * kpoints.shape[0]).reshape((kpoints.shape[0], -1))

nkpoints = len(kpoints)
nbands = mesh_g2.shape[2]
nmodes = mesh_frequencies.shape[1]
energies = np.zeros((nkpoints, nbands), float)
inv_taus = np.zeros((ncalcs, nkpoints, nbands, nmodes), float)
velocities = np.zeros((nkpoints, nbands, 3), float)


with h5py.File('relaxation-times-fine.hdf5', 'w') as f: 

    struct_grp = f.create_group('struct')
    struct_grp.attrs['name'] = 'graphene'
    el_grp = f.create_group('el')
    elph_grp = f.create_group('el-ph')

    ds = struct_grp.create_dataset('positions', data=dftb.primitive.get_positions()* BOHR__AA)
    ds = struct_grp.create_dataset('numbers', data=dftb.primitive.get_atomic_numbers())
    ds = struct_grp.create_dataset('cell', data=dftb.primitive.get_cell()* BOHR__AA)

    print('-- processing %d k-points' % (nkpoints))
    printProgressBar(0, nkpoints, prefix='k-point', suffix='complete')
    for ik, kvec in enumerate(kpoints):
        eps_k, mesh_epskq, mesh_epskmq, mesh_g2 = dftb.calculate_g2(kvec, mesh_qpoints, mesh_frequencies, mesh_eigenvectors)
        eps_k, vel_k = dftb.calculate_velocity(kvec)

        for n in range(nbands):
            energies[ik,n] = eps_k[n]
            velocities[ik,n,:] = vel_k[:,n]
            
            for ic, (kBT, mu, sigma_0) in enumerate(zip(kBTs, mus, sigmas)):
                inv_taus[ic,ik,n,:] = inv_tau_nk_lam( n, eps_k[n], mu, kBT, mesh_g2, mesh_epskq, mesh_frequencies*THZ__EV, sigma=sigma_0)[0]/q_mesh_refinement**2
            
        printProgressBar(ik+1, nkpoints, prefix='k-point', suffix='complete')
    
    ds = elph_grp.create_dataset('linewidths', data=inv_taus)
    ds.attrs['kBTs'] = kBTs
    ds.attrs['mus'] = mus
    ds.attrs['Nq'] = nq
    ds.attrs['sigma'] = sigmas
    
    ds = el_grp.create_dataset('eps_kn', data=energies)
    ds.attrs['EF'] = EF
    ds = el_grp.create_dataset('vel_kn', data=velocities)
    ds = el_grp.create_dataset('k-points', data=kpoints)
