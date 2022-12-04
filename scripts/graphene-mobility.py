import sys, os
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

nq = 100 # size of q-point mesh (nq x nq x 1)

# other parameters below under #6

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
# 5 Calculate mobility directly (still work in progress)
# nk = 2*nq # size of k-point mesh (nk x nk x 1)
# kBT = 0.0259  # temperature in eV
# mu = eps_k[24] + 0.05   # chemical potential in eV
# EF = (eps_k[24]+eps_k[23])/2    # Fermi energy in eV
# sigma_0 = 0.025 # smearing in eV


# print(mu, EF)
# # print((-dfermi_deps((eps_k-mu)/kBT)/kBT))

# # atoms = ase.Atoms(positions=dftb.primitive.get_positions()* BOHR__AA, 
# #                   numbers=dftb.primitive.get_atomic_numbers(), 
# #                   cell=dftb.primitive.get_cell()* BOHR__AA, 
# #                   pbc=[1, 1, 1])
# # fromspglib = spglib.get_ir_reciprocal_mesh([nk, nk, 1], atoms)

# # indices = np.unique(fromspglib[0]).tolist()
# # weights = np.array([list(fromspglib[0]).count(i) for i in indices])
# # kpoints = fromspglib[1] / float(nk)
# # kpoints = kpoints[indices,:]

# gp = GridPoints(np.array([nk,nk,1]), np.linalg.inv(dftb.primitive.cell), 
#                 rotations=None, is_mesh_symmetry=False, is_gamma_center=True, is_time_reversal=False)
# kpoints = gp.qpoints
# weights = gp.weights

# nmeshpoints = np.prod([nk, nk, 1])
# nkpoints = len(kpoints)
# nbands = mesh_g2.shape[2]

# conductivity = np.zeros((3,3), dtype=float)
# conductivity0 = np.zeros((3,3), dtype=float)
# density = 0.0
# density0 = 0.0

# eband = []
# vband = []
# fdband = []

# print('-- processing %d (from %d x %d) k-points' % (nkpoints, nk, nk))
# printProgressBar(0, nkpoints, prefix='k-point', suffix='complete')
# for ik, kvec in enumerate(kpoints):
# #     eps_k, mesh_epskq, mesh_g2 = dftb.calculate_g2(kvec, mesh_qpoints, mesh_frequencies, mesh_eigenvectors)
#     eps_k, vel_k = dftb.calculate_velocity(kvec)

#     for n in range(nbands):        
# #         inv_tau_n = inv_tau_nk( n, eps_k[n], mu, kBT, mesh_g2, mesh_epskq, mesh_frequencies*THZ__EV, sigma=sigma_0)
# #         conductivity += (AA_EV__m_s**2 * EV__ps * 1e-12)*(1/nkpoints)*(-dfermi_deps((eps_k[n]-mu)/kBT)/kBT) * np.outer(vel_k[:,n],vel_k[:,n]) * (1/inv_tau_n)
#         conductivity0 += (AA_EV__m_s**2 * 1e-12)*(weights[ik]/nmeshpoints)*(-dfermi_deps((eps_k[n]-mu)/kBT)/kBT) * np.outer(vel_k[:,n],vel_k[:,n])
#         density += (weights[ik]/nmeshpoints) * ( fermi((eps_k[n]-mu)/kBT) - fermi((eps_k[n]-EF)/kBT) )
#         density0 += (weights[ik]/nmeshpoints) * ( fermi((eps_k[n]-EF)/kBT) )
        
#     eband.append(eps_k)
#     vband.append(vel_k)
#     fdband.append((-dfermi_deps((eps_k-mu)/kBT)/kBT))

#     printProgressBar(ik+1, nkpoints, prefix='k-point', suffix='complete')

# cell = ph.primitive.get_cell() * BOHR__AA
# cell_area = np.linalg.norm(np.cross(cell[0],cell[1]))


# print(2*conductivity0/cell_area)
# print(2*conductivity/cell_area)

# print(2*density/cell_area)
# print(2*density0/cell_area)

# print(conductivity/(density))


# np.savez('eband_vband.npz', eband=np.vstack(eband), vband=np.vstack(vband), fdband=np.vstack(fdband))


# bs_dict = { 'particleType': 'electron', 'numBands': bands[0].shape[1], 'energies': np.vstack(bands), 'energyUnit': 'eV', 
#   'coordsType': 'lattice', 'highSymCoordinates': np.vstack(path), 'highSymLabels': np.array(path_labels).flatten(), 
#   'highSymIndices': np.cumsum(np.array([[0, kps.shape[0]] for kps in kpoints]).flatten()),
#   'wavevectorCoordinates': np.vstack(kpoints), 'wavevectorIndices': list(range(np.vstack(kpoints).shape[0]))}

# with open('path_el_bandstructure.json', 'w') as outfile:
#     json.dump(json.dumps(bs_dict, default=convert), outfile)


# lw_dict = { 'chemicalPotentialUnit': 'eV', 'chemicalPotentials': [mu,], 'particleType': 'electron', 'numBands': bands[0].shape[1], 
#   'energies': np.vstack(bands), 'energyUnit': 'eV', 
#   'coordsType': 'lattice', 'linewidths': np.vstack(linewidths), 'linewidthsUnit': 'eV', 'temperatureUnit': 'eV', 'temperatures': [kBT,]}

# with open('path_el_relaxation_times.json', 'w') as outfile:
#     json.dump(json.dumps(lw_dict, default=convert), outfile)


###############################################################################
# 6 Calculate relaxation times for BT2
nk = nq # size of k-point mesh (nk x nk x 1)
mus = np.linspace(eps_k[3] - 0.150, eps_k[3] + 0.150, 20) # chemical potentials in eV
kBTs = 0.0259 * np.ones_like(mus)  # temperatures in eV
EF = (eps_k[3]+eps_k[4])/2    # Fermi energy in eV
sigma_0 = 0.005 # smearing in eV

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

    ds = struct_grp.create_dataset('positions', data=dftb.primitive.get_positions()* BOHR__AA)
    ds = struct_grp.create_dataset('numbers', data=dftb.primitive.get_atomic_numbers())
    ds = struct_grp.create_dataset('cell', data=dftb.primitive.get_cell()* BOHR__AA)

    print('-- processing %d k-points' % (nkpoints))
    printProgressBar(0, nkpoints, prefix='k-point', suffix='complete')
    for ik, kvec in enumerate(kpoints):
        eps_k, mesh_epskq, mesh_g2 = dftb.calculate_g2(kvec, mesh_qpoints, mesh_frequencies, mesh_eigenvectors)
        eps_k, vel_k = dftb.calculate_velocity(kvec)

        for n in range(nbands):
            energies[ik,n] = eps_k[n]
            velocities[ik,n,:] = vel_k[:,n]
            
            for ic, (kBT, mu) in enumerate(zip(kBTs, mus)):
                inv_taus[ic,ik,n] = inv_tau_nk( n, eps_k[n], mu, kBT, mesh_g2, mesh_epskq, mesh_frequencies*THZ__EV, sigma=sigma_0)[0]
            
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
