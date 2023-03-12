import os, sys
sys.path.insert(0, '../') # adjust path to the base directory of the package

# limit the number of threads for numpy
os.environ["OMP_NUM_THREADS"] = "12" # export OMP_NUM_THREADS=16
os.environ["OPENBLAS_NUM_THREADS"] = "12" # export OPENBLAS_NUM_THREADS=16 
os.environ["MKL_NUM_THREADS"] = "12" # export MKL_NUM_THREADS=16
os.environ["VECLIB_MAXIMUM_THREADS"] = "12" # export VECLIB_MAXIMUM_THREADS=16
os.environ["NUMEXPR_NUM_THREADS"] = "12" # export NUMEXPR_NUM_THREADS=16

###############################################################################
import numpy as np

import phonopy
from phonopy.structure.grid_points import GridPoints
from timeit import default_timer as timer

from dftbephy import DftbSuperCellCalc
from dftbephy.analysis import inv_tau_nk_lam
from dftbephy.units import *
from dftbephy.tools import printProgressBar

# for writing hdf files
import h5py

# for reading hsd files
import hsd

# parallelization
from mpi4py import MPI

###############################################################################
comm = MPI.COMM_WORLD
world_size = comm.Get_size()
rank = comm.Get_rank()
start_time = MPI.Wtime()

if rank==0:
    print('-- working with %i processes' % (world_size))
    sys.stdout.flush()
    
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

# read section for relaxation time calculations
rt_dict = check_hsd_input(inp_dict, 'RelaxationTimes')

if 'SERTA'in rt_dict.keys():
    rta_method = 'SERTA'
    rt_dict = rt_dict['SERTA']
else:
    if rank==0:
        print('-- unsupported method for RelaxationTimes.')
    quit()

default_mesh = {'Mesh': {'npoints': [1, 1, 1], 'refinement': 1, 'shift': [0.0,0.0,0.0]}}
qp_dict = rt_dict.get('qpoints', default_mesh)
if 'Mesh' in qp_dict.keys():
    # size of q-point mesh (eg. nq x nq x 1)
    q_mesh = qp_dict['Mesh'].get('npoints', default_mesh['Mesh']['npoints'])
    # factor by which the q-points are scaled
    q_mesh_refinement = qp_dict['Mesh'].get('refinement', default_mesh['Mesh']['refinement'])
    q_mesh_shift = qp_dict['Mesh'].get('shift', default_mesh['Mesh']['shift'])

kp_dict = rt_dict.get('kpoints', default_mesh)
if 'Mesh' in kp_dict.keys():
    # size of k-point mesh (eg. nk x nk x 1)
    k_mesh = kp_dict['Mesh'].get('npoints', default_mesh['Mesh']['npoints'])
    # factor by which the k-points are scaled
    k_mesh_refinement = kp_dict['Mesh'].get('refinement', default_mesh['Mesh']['refinement'])
    k_mesh_shift = kp_dict['Mesh'].get('shift', default_mesh['Mesh']['shift'])
    
kvec0 = np.array(k_mesh_shift) # reference k-point for electrons

band_sel = rt_dict.get('bands', None)
if type(band_sel) is list:
    band_sel[1] = band_sel[1]+1

if type(rt_dict['mu']) is dict:
    m = rt_dict['mu']['Range']
    mu_list = np.arange(float(m[0]), float(m[1]), int(m[2]))
elif type(rt_dict['mu']) is list:
    mu_list = rt_dict['mu']
elif type(rt_dict['mu']) is float:
    mu_list = [rt_dict['mu']]
else:
    if rank==0:
        print('-- unkown type for mu (neither Range, list or float).')
        quit()

kBT0 = rt_dict.get('temperature', 0.0259)
assert(type(kBT0) is float)

sigma0 = rt_dict.get('sigma', 0.003)
assert(type(sigma0) is float)

EF = rt_dict.get('Efermi', 0.00)
assert(type(EF) is float)


###############################################################################
# 1 Load phonopy
if rank == 0:
    print('-- looking for phonopy results in %s' % (basedir + phonopy_dir))
os.chdir(basedir + phonopy_dir)
ph = phonopy.load('phonopy_disp.yaml')

###############################################################################
# 2 Initialize DFTB calculator
if rank == 0:
    print('-- working in %s' % (basedir + working_dir))
os.chdir(basedir + working_dir)

dftb = DftbSuperCellCalc(angular_momenta)
dftb.load_phonopy(ph)

if os.path.isfile('reference.npz'):
    print('-- [%i] loading reference calculation' % (rank))
    npzfile = np.load('reference.npz')
    dftb.H0 = npzfile['H0']
    dftb.S0 = npzfile['S0']
else:
    print('-- run reference calculation in serial!' )
sys.stdout.flush()
###############################################################################
# 3 Run band-structure calculation (optional)

###############################################################################
# 4 Run calculations of electron-phonon couplings

if os.path.isfile('derivatives.npz'):
    print('-- [%i] loading derivatives calculation' % (rank))
    npzfile = np.load('derivatives.npz')
    dftb.H_derivs = npzfile['H_derivs']
    dftb.S_derivs = npzfile['S_derivs']
else:
    print('-- run reference calculation in serial!' )
sys.stdout.flush()

print('-- [%i] starting phonon calculations on mesh ...' % (rank))
start = MPI.Wtime()
ph.run_mesh(q_mesh, with_eigenvectors=False, is_mesh_symmetry=False, is_gamma_center=False)
mesh = ph.get_mesh_dict()
end = MPI.Wtime()
print('-- [%i] finished (%4.1f s).' % (rank, end-start))
sys.stdout.flush()

# we are only interested in a small region around kvec0
qps = mesh['qpoints']/q_mesh_refinement
ph.run_qpoints(qps, with_eigenvectors=True)

mesh = ph.get_qpoints_dict()
mesh_qpoints = qps
mesh_frequencies = mesh['frequencies'] 
mesh_eigenvectors = mesh['eigenvectors']

# do one calculation of epcs at kvec0
print('-- [%i] starting el-ph calculation ...' % (rank))
start = MPI.Wtime()
eps_k, mesh_epskq, mesh_epskmq, mesh_g2 = dftb.calculate_g2(kvec0, mesh_qpoints, mesh_frequencies, mesh_eigenvectors, band_sel=band_sel)
end = MPI.Wtime()
print('-- [%i] finished (%4.1f s).' % (rank, end-start))
sys.stdout.flush()

###############################################################################
# 6 Calculate relaxation times

#EF = (eps_k[nbv]+eps_k[nbc])/2    # Fermi energy in eV
#mus = np.array([eps_k[nbv]-0.1, eps_k[nbv], EF, eps_k[nbc], eps_k[nbc] + 0.10]) # chemical potentials in eV

mus = EF + np.array(mu_list) # chemical potentials
sigmas = sigma0* np.ones_like(mus) # smearing
kBTs = kBT0 * np.ones_like(mus)   # temperatures

# all calculations
temps_chems = zip(kBTs, mus)
ncalcs = len(list(temps_chems))
if rank == 0:
    print('-- running %d calculations for mu, kBT, sigma:' % (ncalcs))
    for ic, (kBT, mu, sigma_0) in enumerate(zip(kBTs, mus, sigmas)):
        print(ic, mu, kBT, sigma_0)
    sys.stdout.flush()

gp = GridPoints(np.array(k_mesh), np.linalg.inv(dftb.primitive.cell), 
                 rotations=None, is_mesh_symmetry=False, is_gamma_center=False, is_time_reversal=False)
kpoints = gp.qpoints / k_mesh_refinement
kpoints += np.array(k_mesh_shift * kpoints.shape[0]).reshape((kpoints.shape[0], -1))

nkpoints = len(kpoints)
nbands = mesh_g2.shape[2]
nmodes = mesh_frequencies.shape[1]

if rank ==0:
    energies = np.zeros((nkpoints, nbands), dtype='float64')
    inv_taus = np.zeros((nkpoints, ncalcs, nbands, nmodes), dtype='float64')
    velocities = np.zeros((nkpoints, nbands, 3), dtype='float64')
else:
    energies = None
    inv_taus = None
    velocities = None

ks_per_rank = np.array_split(np.arange(nkpoints), world_size)
num_ks_per_rank = np.array([len(rpr) for rpr in ks_per_rank])
energies_displ_per_rank = np.insert(np.cumsum(num_ks_per_rank * nbands),0,0)[0:-1]
inv_taus_displ_per_rank = np.insert(np.cumsum(num_ks_per_rank * ncalcs*nbands*nmodes),0,0)[0:-1]
velocities_displ_per_rank = np.insert(np.cumsum(num_ks_per_rank * nbands * 3),0,0)[0:-1]

loc_energies = np.zeros((num_ks_per_rank[rank], nbands), dtype='float64')
loc_inv_taus = np.zeros((num_ks_per_rank[rank], ncalcs, nbands, nmodes), dtype='float64')
loc_velocities = np.zeros((num_ks_per_rank[rank], nbands, 3), dtype='float64')


print('-- [%i] processing %d k-points' % (rank, num_ks_per_rank[rank]))
sys.stdout.flush()
for loc_ik, ik in enumerate(ks_per_rank[rank]):
    kvec = kpoints[ik]
    
    eps_k, mesh_epskq, mesh_epskmq, mesh_g2 = dftb.calculate_g2(kvec, mesh_qpoints, mesh_frequencies, mesh_eigenvectors, band_sel=band_sel)
    eps_k, vel_k = dftb.calculate_velocity(kvec, band_sel=band_sel)

    for n in range(nbands):
        loc_energies[loc_ik,n] = eps_k[n]
        loc_velocities[loc_ik,n,:] = vel_k[:,n]

        for ic, (kBT, mu, sigma_0) in enumerate(zip(kBTs, mus, sigmas)):
            loc_inv_taus[loc_ik,ic,n,:] = inv_tau_nk_lam( n, eps_k[n], mu, kBT, mesh_g2, mesh_epskq, mesh_frequencies*THZ__EV, sigma=sigma_0)[0]/q_mesh_refinement**2

print('-- [%i] done.' % (rank))
sys.stdout.flush()

comm.Gatherv(loc_energies, [energies, num_ks_per_rank * nbands, energies_displ_per_rank, MPI.DOUBLE], root=0)
comm.Gatherv(loc_inv_taus, [inv_taus, num_ks_per_rank * ncalcs*nbands*nmodes, inv_taus_displ_per_rank, MPI.DOUBLE], root=0)
comm.Gatherv(loc_velocities, [velocities, num_ks_per_rank * nbands * 3, velocities_displ_per_rank, MPI.DOUBLE], root=0)

if rank==0:
    
    inv_taus_ = np.zeros((ncalcs, nkpoints, nbands, nmodes), dtype='float64')
    for nc in range(ncalcs):
        for ik in range(nkpoints):
            for n in range(nbands):
                for lam in range(nmodes):
                    inv_taus_[nc, ik, n, lam] = inv_taus[ik, nc, n, lam]
    
    print('-- writing results to %s' % (basedir + results_dir))
    with h5py.File(basedir + results_dir + 'relaxation-times-fine-mpi-bandsel.hdf5', 'w') as f: 

        struct_grp = f.create_group('struct')
        struct_grp.attrs['name'] = name
        el_grp = f.create_group('el')
        elph_grp = f.create_group('el-ph')

        ds = struct_grp.create_dataset('positions', data=dftb.primitive.get_positions()* BOHR__AA)
        ds = struct_grp.create_dataset('numbers', data=dftb.primitive.get_atomic_numbers())
        ds = struct_grp.create_dataset('cell', data=dftb.primitive.get_cell()* BOHR__AA)
        
        ds = elph_grp.create_dataset('linewidths', data=inv_taus_)
        ds.attrs['kBTs'] = kBTs
        ds.attrs['mus'] = mus
        ds.attrs['Nq'] = np.array(q_mesh)
        ds.attrs['sigma'] = sigmas

        ds = el_grp.create_dataset('eps_kn', data=energies)
        ds.attrs['EF'] = EF
        ds = el_grp.create_dataset('vel_kn', data=velocities)
        ds = el_grp.create_dataset('k-points', data=kpoints)
        
    end_time = MPI.Wtime()

    print('-- total runtime: %f s' % (end_time-start_time))
