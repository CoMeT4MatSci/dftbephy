import os, sys
sys.path.insert(0, '../') # adjust path to the base directory of the package

import numpy as np
import matplotlib.pyplot as plt

import phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
from phonopy.structure.grid_points import GridPoints
from timeit import default_timer as timer

from dftbephy import DftbSuperCellCalc
from dftbephy.fileio import read_dftb_bands, get_lumo
from dftbephy.analysis import inv_tau_nk
from dftbephy.units import *
from dftbephy.tools import printProgressBar

import ase
import spglib

# this is needed for writing hdf files
import h5py

# this is needed for writing json files
import json

def convert(x):
    if hasattr(x, "tolist"):  # numpy arrays have this
        return x.tolist()
    raise TypeError(x)

# _FD_XMAX = np.log(np.finfo(float).max)/2
_FD_XMAX = 18.42068 # from BoltzTrap
    
def fermi(x):
    f = np.where(x < 0., 1., 0.)
    norm_ind = np.logical_and(x > -_FD_XMAX, x < _FD_XMAX)
    f[norm_ind] = 1/(np.exp(x[norm_ind]) + 1.)
    return f

def dfermi_deps(x):
    f = np.zeros_like(x)
    norm_ind = np.logical_and(x > -_FD_XMAX, x < _FD_XMAX)
    f[norm_ind] = -np.exp(x[norm_ind])/(np.exp(x[norm_ind]) + 1.)**2
    return f
    
###############################################################################
# 0 Set base directory
basedir = '../examples/Graphene'

q_mesh_refinement = 10

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
ph.run_mesh([nq,nq,1], with_eigenvectors=False, is_mesh_symmetry=False)
mesh = ph.get_mesh_dict()
end = timer()
print('-- finished (%4.1f s).' % (end-start))

#mesh_qpoints= mesh['qpoints']
#mesh_frequencies = mesh['frequencies'] 
#mesh_eigenvectors = mesh['eigenvectors']

# we are only interested in a small region around kvec0
qps = mesh['qpoints']/q_mesh_refinement
ph.run_qpoints(qps, with_eigenvectors=True)

mesh = ph.get_qpoints_dict()
mesh_qpoints = qps
mesh_frequencies = mesh['frequencies'] 
mesh_eigenvectors = mesh['eigenvectors']

kvec0 = np.array([-2/3, -1/3, 0]) # k-point for electrons

print('-- starting el-ph calculation ...')
start = timer()
eps_k, mesh_epskq, mesh_epskmq, mesh_g2 = dftb.calculate_g2(kvec0, mesh_qpoints, mesh_frequencies, mesh_eigenvectors)
end = timer()
print('-- finished (%4.1f s).' % (end-start))


###############################################################################
# 6 Calculate relaxation times
nk = 200 # size of k-point mesh (nk x nk x 1)
mus = np.linspace(eps_k[4], eps_k[4] + 0.40, 40) # chemical potentials in eV
kBTs = 0.0259 * np.ones_like(mus)  # temperatures in eV
EF = (eps_k[3]+eps_k[4])/2    # Fermi energy in eV
sigma_0 = 0.003 # smearing in eV
Ecut = 1. # cutoff in eV

temps_chems = zip(kBTs, mus)
ncalcs = len(list(zip(kBTs, mus)))

print('-- running %d calculations for mu, kBT:' % (ncalcs))
for ic, (kBT, mu) in enumerate(zip(kBTs, mus)):
    print(ic, mu, kBT)

atoms = ase.Atoms(positions=dftb.primitive.get_positions()* BOHR__AA, 
                  numbers=dftb.primitive.get_atomic_numbers(), 
                  cell=dftb.primitive.get_cell()* BOHR__AA, 
                  pbc=[1, 1, 1])
fromspglib = spglib.get_ir_reciprocal_mesh([nk, nk, 1], atoms)

indices = np.unique(fromspglib[0]).tolist()
weights = [list(fromspglib[0]).count(i) for i in indices]
kpoints = fromspglib[1] / float(nk)
kpoints = kpoints[indices,:]
# gp = GridPoints(np.array([nk,nk,1]), np.linalg.inv(dftb.primitive.cell), 
#                 rotations=None, is_mesh_symmetry=True, is_gamma_center=True, is_time_reversal=True)

nkpoints = len(kpoints)
nbands = mesh_g2.shape[2]
energies = np.zeros((nkpoints, nbands), float)
inv_taus = np.zeros((ncalcs, nkpoints, nbands), float)
velocities = np.zeros((nkpoints, nbands, 3), float)

with h5py.File('relaxation-times.hdf5', 'w') as f: 

    struct_grp = f.create_group('struct')
    struct_grp.attrs['name'] = 'gamma-graphyne'
    el_grp = f.create_group('el')
    elph_grp = f.create_group('el-ph')
    ph_grp = f.create_group('ph')    

    ds = struct_grp.create_dataset('positions', data=dftb.primitive.get_positions()* BOHR__AA)
    ds = struct_grp.create_dataset('numbers', data=dftb.primitive.get_atomic_numbers())
    ds = struct_grp.create_dataset('cell', data=dftb.primitive.get_cell()* BOHR__AA)

    ds = ph_grp.create_dataset('q-points', data=mesh_qpoints)
    ds = ph_grp.create_dataset('frequencies', data=mesh_frequencies*THZ__EV)
    
    print('-- processing %d k-points' % (nkpoints))
    printProgressBar(0, nkpoints, prefix='k-point', suffix='complete')
    nepccalcs = 0
    for ik, kvec in enumerate(kpoints):
        eps_k, vel_k = dftb.calculate_velocity(kvec)

        # consider only k-points within a certain energy range
        if np.abs(eps_k - EF).min() < Ecut:
        
            eps_k, mesh_epskq, mesh_epskmq, mesh_g2 = dftb.calculate_g2(kvec, mesh_qpoints, mesh_frequencies, mesh_eigenvectors)
            nepccalcs += 1
            
            for n in range(nbands):
                energies[ik,n] = eps_k[n]
                velocities[ik,n,:] = vel_k[:,n]

                for ic, (kBT, mu) in enumerate(zip(kBTs, mus)):
                    inv_taus[ic,ik,n] = inv_tau_nk( n, eps_k[n], mu, kBT, mesh_g2, mesh_epskq, mesh_frequencies*THZ__EV, sigma=sigma_0)[0]/q_mesh_refinement**2
                
        else:
            for n in range(nbands):
                energies[ik,n] = eps_k[n]
                velocities[ik,n,:] = vel_k[:,n]
                inv_taus[:,ik,n] = 1e10
            
        printProgressBar(ik+1, nkpoints, prefix='k-point', suffix='complete')
    
    ds = elph_grp.create_dataset('linewidths', data=inv_taus)
    ds.attrs['kBTs'] = kBTs
    ds.attrs['mus'] = mus
    ds.attrs['Nq'] = nq
    ds.attrs['sigma'] = sigma_0
    
    ds = el_grp.create_dataset('eps_kn', data=energies)
    ds.attrs['EF'] = EF
    ds = el_grp.create_dataset('vel_kn', data=velocities)
    ds = el_grp.create_dataset('k-points', data=kpoints)

    print('-- %d epc calculations ' % (nepccalcs))