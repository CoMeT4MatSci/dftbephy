import os, sys
sys.path.insert(0, '../') # adjust path to the base directory of the package

# limit the number of threads for numpy
os.environ["OMP_NUM_THREADS"] = "16" # export OMP_NUM_THREADS=16
os.environ["OPENBLAS_NUM_THREADS"] = "16" # export OPENBLAS_NUM_THREADS=16 
os.environ["MKL_NUM_THREADS"] = "16" # export MKL_NUM_THREADS=16
os.environ["VECLIB_MAXIMUM_THREADS"] = "16" # export VECLIB_MAXIMUM_THREADS=16
os.environ["NUMEXPR_NUM_THREADS"] = "16" # export NUMEXPR_NUM_THREADS=16

###############################################################################
import numpy as np

import phonopy
from phonopy.structure.grid_points import GridPoints
from timeit import default_timer as timer

from dftbephy import DftbSuperCellCalc
from dftbephy.analysis import inv_tau_nk
from dftbephy.units import *
from dftbephy.tools import printProgressBar

#import ase

# for irreducible mesh
import spglib

# for reading hsd files
import hsd

# this is needed for writing json files
import json

# parallelization
from mpi4py import MPI

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

def step(x):
    f = np.where(x < 0., 1., 0.)
    return f

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

# read section for conductivity calculations
mob_dict = check_hsd_input(inp_dict, 'Conductivities')

if 'CRTA' in mob_dict.keys():
    rta_method = 'CRTA'
    mob_dict = mob_dict['CRTA']
elif 'SERTA' in mob_dict.keys():
    rta_method = 'SERTA'
    mob_dict = mob_dict['SERTA']
else:
    if rank==0:
        print('-- unsupported method for RelaxationTimes.')
    quit()

default_mesh = {'Mesh': {'npoints': [1, 1, 1], 'refinement': 1, 'shift': [0.0,0.0,0.0]}}

qp_dict = mob_dict.get('qpoints', default_mesh)
if 'Mesh' in qp_dict.keys():
    # size of q-point mesh (eg. nq x nq x 1)
    q_mesh = qp_dict['Mesh'].get('npoints', default_mesh['Mesh']['npoints'])
    # factor by which the q-points are scaled
    q_mesh_refinement = qp_dict['Mesh'].get('refinement', default_mesh['Mesh']['refinement'])
    q_mesh_shift = qp_dict['Mesh'].get('shift', default_mesh['Mesh']['shift'])

kp_dict = mob_dict.get('kpoints', default_mesh)
if 'Mesh' in kp_dict.keys():
    # size of k-point mesh (eg. nk x nk x 1)
    k_mesh = kp_dict['Mesh'].get('npoints', default_mesh['Mesh']['npoints'])
    # factor by which the k-points are scaled
    k_mesh_refinement = kp_dict['Mesh'].get('refinement', default_mesh['Mesh']['refinement'])
    k_mesh_shift = kp_dict['Mesh'].get('shift', default_mesh['Mesh']['shift'])
    
kvec0 = np.array(k_mesh_shift) # reference k-point for electrons

band_sel = mob_dict.get('bands', None)
if type(band_sel) is list:
    band_sel[1] = band_sel[1]+1

if type(mob_dict['mu']) is dict:
    m = mob_dict['mu']['Range']
    mu_list = np.linspace(float(m[0]), float(m[1]), int(m[2]))
elif type(mob_dict['mu']) is list:
    mu_list = mob_dict['mu']
elif type(mob_dict['mu']) is float:
    mu_list = [mob_dict['mu']]
else:
    if rank==0:
        print('-- unkown type for mu (neither Range, list or float).')
    quit()

kBT0 = mob_dict.get('temperature', 0.0259)
assert(type(kBT0) is float)

sigma0 = mob_dict.get('sigma', 0.003)
assert(type(sigma0) is float)

EF = mob_dict.get('Efermi', 0.00)
assert(type(EF) is float)

Ecut = mob_dict.get('Ecut', 1.0)
assert(type(Ecut) is float)

spin_factor = mob_dict.get('SpinDegeneracy', 1)
assert(type(spin_factor) is int)

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
#
# only for SERTA

if rta_method == 'SERTA':
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


    print('-- [%i] starting el-ph calculation ...' % (rank))
    start = MPI.Wtime()
    eps_k, mesh_epskq, mesh_epskmq, mesh_g2 = dftb.calculate_g2(kvec0, mesh_qpoints, mesh_frequencies, mesh_eigenvectors, band_sel=band_sel)
    end = MPI.Wtime()
    print('-- [%i] finished (%4.1f s).' % (rank, end-start))
    sys.stdout.flush()

###############################################################################
# 5 Calculate conductivity

mus = EF + np.array(mu_list)  # chemical potentials
kBTs = kBT0*np.ones_like(mus) # temperatures
sigma_0 = sigma0 # smearing

temps_chems = zip(kBTs, mus)
ncalcs = len(list(temps_chems))
if rank == 0:
    print('-- running %d calculations for mu, kBT:' % (ncalcs))
    for ic, (kBT, mu) in enumerate(zip(kBTs, mus)):
        print(ic, mu, kBT)
    sys.stdout.flush()

    print('-- constructing k-mesh')
#    atoms = ase.Atoms(positions=dftb.primitive.get_positions()* BOHR__AA, 
#                      numbers=dftb.primitive.get_atomic_numbers(), 
#                      cell=dftb.primitive.get_cell()* BOHR__AA, 
#                      pbc=[1, 1, 1])
    fromspglib = spglib.get_ir_reciprocal_mesh(k_mesh, dftb.primitive)
    
    indices, weights = np.unique(fromspglib[0], return_counts=True)
    weights = np.asarray(weights, dtype='int')
    
#    indices = np.unique(fromspglib[0]).tolist()
#    weights = np.array([list(fromspglib[0]).count(i) for i in indices], dtype='int')
    kpoints = fromspglib[1]
    kpoints = kpoints[indices,:] / np.array(k_mesh)
    
    nkpoints = len(kpoints)
    
    print('-- number of irreducible k-points: %i' % (nkpoints))
else:
    nkpoints = None
    
    print('-- [%i] waiting for Godot ...' % (rank))
    sys.stdout.flush()

nkpoints = comm.bcast(nkpoints, root=0)

if rank != 0:
    kpoints = np.empty((nkpoints,3), dtype='float64')
    weights = np.empty(nkpoints, dtype='int')

comm.Bcast([kpoints, MPI.DOUBLE], root=0)
comm.Bcast([weights, MPI.INT], root=0)

    
nmeshpoints = np.prod(k_mesh)

if rta_method == 'SERTA':
    nbands = mesh_g2.shape[2]
elif (rta_method == 'CRTA') and (band_sel is None):
    nbands = dftb.get_num_bands()
elif (rta_method == 'CRTA') and (type(band_sel) is list):
    nbands = band_sel[1] - band_sel[0]

conductivities = np.zeros((ncalcs, 3, 3), float)
conductivities0 = np.zeros((ncalcs, 3, 3), float)
densities = np.zeros((ncalcs, ), float)
densities0 = np.zeros((ncalcs, ), float)

loc_conductivities = np.zeros((ncalcs, 3, 3), float)
loc_conductivities0 = np.zeros((ncalcs, 3, 3), float)
loc_densities = np.zeros((ncalcs, ), float)
loc_densities0 = np.zeros((ncalcs, ), float)


if rank ==0:
    energies = np.zeros((nkpoints, nbands), dtype='float64')
    velocities = np.zeros((nkpoints, nbands, 3), dtype='float64')
else:
    energies = None
    velocities = None

ks_per_rank = np.array_split(np.arange(nkpoints), world_size)
num_ks_per_rank = np.array([len(rpr) for rpr in ks_per_rank])
energies_displ_per_rank = np.insert(np.cumsum(num_ks_per_rank * nbands),0,0)[0:-1]
velocities_displ_per_rank = np.insert(np.cumsum(num_ks_per_rank * nbands * 3),0,0)[0:-1]

loc_energies = np.zeros((num_ks_per_rank[rank], nbands), dtype='float64')
loc_velocities = np.zeros((num_ks_per_rank[rank], nbands, 3), dtype='float64')

    
print('-- [%i] processing %d k-points' % (rank, num_ks_per_rank[rank]))
sys.stdout.flush()
nepccalcs = 0
for loc_ik, ik in enumerate(ks_per_rank[rank]):
    kvec = kpoints[ik]

    eps_k, vel_k = dftb.calculate_velocity(kvec, band_sel=band_sel)

    for n in range(nbands):
        loc_energies[loc_ik,n] = eps_k[n]
        loc_velocities[loc_ik,n,:] = vel_k[:,n]

        for ic, (kBT, mu) in enumerate(zip(kBTs, mus)):
            loc_densities[ic]  += (weights[ik]/nmeshpoints) * ( fermi((eps_k[n]-mu)/kBT) - step(eps_k[n]-EF) )
            loc_densities0[ic] += (weights[ik]/nmeshpoints) * ( fermi((eps_k[n]-EF)/kBT) )                    
    
    # consider only k-points within a certain energy range for transport properties
    if np.abs(eps_k - EF).min() < Ecut:

        if rta_method == 'SERTA':
            eps_k, mesh_epskq, mesh_epskmq, mesh_g2 = dftb.calculate_g2(kvec, mesh_qpoints, mesh_frequencies, mesh_eigenvectors, band_sel=band_sel)
            nepccalcs += 1

        for ic, (kBT, mu) in enumerate(zip(kBTs, mus)):
            for n in range(nbands):

                if rta_method == 'SERTA':
                    inv_tau = inv_tau_nk( n, eps_k[n], mu, kBT, mesh_g2, mesh_epskq, mesh_frequencies*THZ__EV, sigma=sigma_0)[0]/q_mesh_refinement**2
                
                    loc_conductivities[ic,:,:]  += (AA_EV__m_s**2 * EV__ps * 1e-12)*(weights[ik]/nmeshpoints)*(-dfermi_deps((eps_k[n]-mu)/kBT)/kBT) * np.outer(vel_k[:,n],vel_k[:,n]) * (1/inv_tau)
                    
                loc_conductivities0[ic,:,:] += (AA_EV__m_s**2 * 1e-12)*(weights[ik]/nmeshpoints)*(-dfermi_deps((eps_k[n]-mu)/kBT)/kBT) * np.outer(vel_k[:,n],vel_k[:,n])

print('-- [%i] %d epc calculations done.' % (rank, nepccalcs))
sys.stdout.flush()

comm.Gatherv(loc_energies, [energies, num_ks_per_rank * nbands, energies_displ_per_rank, MPI.DOUBLE], root=0)
comm.Gatherv(loc_velocities, [velocities, num_ks_per_rank * nbands * 3, velocities_displ_per_rank, MPI.DOUBLE], root=0)

comm.Reduce( [loc_densities, MPI.DOUBLE], [densities, MPI.DOUBLE], op = MPI.SUM, root = 0 )
comm.Reduce( [loc_densities0, MPI.DOUBLE], [densities0, MPI.DOUBLE], op = MPI.SUM, root = 0 )
comm.Reduce( [loc_conductivities, MPI.DOUBLE], [conductivities, MPI.DOUBLE], op = MPI.SUM, root = 0 )
comm.Reduce( [loc_conductivities0, MPI.DOUBLE], [conductivities0, MPI.DOUBLE], op = MPI.SUM, root = 0 )


if rank==0:
    cell = ph.primitive.get_cell() * BOHR__AA
    cell_area = np.linalg.norm(np.cross(cell[0],cell[1]))

    # generate output as json
    tr_dict = { 'name': name, 'kBTs': kBTs, 'mus': mus, 'Nq': np.array(q_mesh), 'Nk': np.array(k_mesh), 
               'sigma': sigma_0, 'EF': EF, 'cell_area': cell_area,
               'nepccalcs': nepccalcs,
      'conductivities': spin_factor*conductivities/cell_area,
      'conductivities0': spin_factor*conductivities0/cell_area,
      'densities': spin_factor*densities/cell_area,
      'densities0': spin_factor*densities0/cell_area}
    
    print('-- writing results to %s' % (basedir + results_dir))
    with open(basedir + results_dir + 'transport-mu-mpi.json', 'w') as outfile:
        json.dump(json.dumps(tr_dict, default=convert), outfile)
    
    end_time = MPI.Wtime()
    print('-- total runtime: %f s' % (end_time-start_time))