import os, sys

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

import spglib

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
ph.run_mesh([nq,nq,1], with_eigenvectors=False, is_mesh_symmetry=False)
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


kvec0 = np.array([-2/3, -1/3, 0]) # k-point for electrons

print('-- starting el-ph calculation ...')
start = timer()
eps_k, mesh_epskq, mesh_epskmq, mesh_g2 = dftb.calculate_g2(kvec0, mesh_qpoints, mesh_frequencies, mesh_eigenvectors)
end = timer()
print('-- finished (%4.1f s).' % (end-start))

###############################################################################
# 5 Calculate mobility
nk = 400 # size of k-point mesh (nk x nk x 1)
mus = np.linspace(eps_k[4], eps_k[4] + 0.40, 40) # chemical potentials in eV
kBTs = 0.0259 * np.ones_like(mus)  # temperatures in eV
EF = (eps_k[3]+eps_k[4])/2    # Fermi energy in eV
sigma_0 = 0.003 # smearing in eV
Ecut = 1. # cutoff in eV
spin_factor  = 2 # spin degeneracy


temps_chems = zip(kBTs, mus)
ncalcs = len(list(temps_chems))

print('-- running %d calculations for mu, kBT:' % (ncalcs))
for ic, (kBT, mu) in enumerate(zip(kBTs, mus)):
    print(ic, mu, kBT)

print('-- constructing k-mesh')
cell = dftb.primitive.get_cell()* BOHR__AA, \
       dftb.primitive.get_positions()* BOHR__AA, \
       dftb.primitive.get_atomic_numbers()

fromspglib = spglib.get_ir_reciprocal_mesh([nk, nk, 1], cell)

indices, weights = np.unique(fromspglib[0], return_counts=True)
weights = np.asarray(weights, dtype='int')
kpoints = fromspglib[1]
kpoints = kpoints[indices,:] / np.array([nk, nk, 1])

nmeshpoints = np.prod([nk, nk, 1])
nkpoints = len(kpoints)
nbands = mesh_g2.shape[2]

energies = np.zeros((nkpoints, nbands), float)
velocities = np.zeros((nkpoints, nbands, 3), float)
conductivities = np.zeros((ncalcs, 3, 3), float)
conductivities0 = np.zeros((ncalcs, 3, 3), float)
densities = np.zeros((ncalcs, ), float)
densities0 = np.zeros((ncalcs, ), float)

    
print('-- processing %d k-points' % (nkpoints))
printProgressBar(0, nkpoints, prefix='k-point', suffix='complete')
nepccalcs = 0
for ik, kvec in enumerate(kpoints):
    eps_k, vel_k = dftb.calculate_velocity(kvec)

    for n in range(nbands):
        energies[ik,n] = eps_k[n]
        velocities[ik,n,:] = vel_k[:,n]

        for ic, (kBT, mu) in enumerate(zip(kBTs, mus)):
            densities[ic]  += (weights[ik]/nmeshpoints) * ( fermi((eps_k[n]-mu)/kBT) - fermi((eps_k[n]-EF)/kBT) )
            densities0[ic] += (weights[ik]/nmeshpoints) * ( fermi((eps_k[n]-EF)/kBT) )                    
    
    # consider only k-points within a certain energy range for transport properties
    if np.abs(eps_k - EF).min() < Ecut:

        eps_k, mesh_epskq, mesh_epskmq, mesh_g2 = dftb.calculate_g2(kvec, mesh_qpoints, mesh_frequencies, mesh_eigenvectors)
        nepccalcs += 1

        for ic, (kBT, mu) in enumerate(zip(kBTs, mus)):
            for n in range(nbands):

                inv_tau = inv_tau_nk( n, eps_k[n], mu, kBT, mesh_g2, mesh_epskq, mesh_frequencies*THZ__EV, sigma=sigma_0)[0]/q_mesh_refinement**2
                
                conductivities[ic,:,:]  += (AA_EV__m_s**2 * EV__ps * 1e-12)*(weights[ik]/nmeshpoints)*(-dfermi_deps((eps_k[n]-mu)/kBT)/kBT) * np.outer(vel_k[:,n],vel_k[:,n]) * (1/inv_tau)
                conductivities0[ic,:,:] += (AA_EV__m_s**2 * 1e-12)*(weights[ik]/nmeshpoints)*(-dfermi_deps((eps_k[n]-mu)/kBT)/kBT) * np.outer(vel_k[:,n],vel_k[:,n])

    printProgressBar(ik+1, nkpoints, prefix='k-point', suffix='complete')
    
print('-- %d epc calculations done.' % (nepccalcs))

cell = ph.primitive.get_cell() * BOHR__AA
cell_area = np.linalg.norm(np.cross(cell[0],cell[1]))


# generate output as json
tr_dict = { 'kBTs': kBTs, 'mus': mus, 'Nq': nq, 'Nk': nk, 'sigma': sigma_0, 'EF': EF, 'cell_area': cell_area,
           'nepccalcs': nepccalcs,
  'conductivities': spin_factor*conductivities/cell_area,
  'conductivities0': spin_factor*conductivities0/cell_area,
  'densities': spin_factor*densities/cell_area,
  'densities0': spin_factor*densities0/cell_area}

with open('transport.json', 'w') as outfile:
    json.dump(json.dumps(tr_dict, default=convert), outfile)

#np.savez('eband_vband.npz', eband=(energies), vband=(velocities))

