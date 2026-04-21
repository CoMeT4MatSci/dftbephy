import os, sys
import argparse

# limit the number of threads for numpy
# TODO: use threadpoolctl to manage threads
os.environ["OMP_NUM_THREADS"] = "16" # export OMP_NUM_THREADS=16
os.environ["OPENBLAS_NUM_THREADS"] = "16" # export OPENBLAS_NUM_THREADS=16 
os.environ["MKL_NUM_THREADS"] = "16" # export MKL_NUM_THREADS=16
os.environ["VECLIB_MAXIMUM_THREADS"] = "16" # export VECLIB_MAXIMUM_THREADS=16
os.environ["NUMEXPR_NUM_THREADS"] = "16" # export NUMEXPR_NUM_THREADS=16

###############################################################################
import numpy as np

import phonopy

import dftbephy
from dftbephy import DftbSuperCellCalc
from dftbephy.analysis import inv_tau_nk, inv_tau_nk_lam
from dftbephy.units import *
from dftbephy.tools import printProgressBar
from dftbephy.dftbephy_cli import check_hsd_input, convert

# for reading hsd files
import hsd

# parallelization
from mpi4py import MPI

###############################################################################
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
# commands and descritions for argparse
commands = [('relaxationtimes', 'Calculate relaxation times.'),
            ('mobility', 'Calculate transport properties.')]


def main(arguments=None):
    ############# PARSE COMMAND LINE ARGUMENTS #############
    parser = argparse.ArgumentParser(
                        prog='dftbephy-mpi',
                        description='Calculating electron-phonon couplings (EPCs) with DFTB.')

    parser.add_argument('-v', '--verbose', action='store_true')  # on/off flag

    subparsers = parser.add_subparsers(title='Commands', dest='command')

    for (cmd, desc) in commands:
        p = subparsers.add_parser(cmd, description=desc)

    args = parser.parse_args()


    ############# START #############
    comm = MPI.COMM_WORLD
    world_size = comm.Get_size()
    rank = comm.Get_rank()

    if rank==0:
        #
        print('     _______  __            __       ')
        print(' ___/ / _/ /_/ /  ___ ___  / /  __ __')
        print('/ _  / _/ __/ _ \\/ -_) _ \\/ _ \\/ // /')
        print('\\_,_/_/ \\__/_.__/\\__/ .__/_//_/\\_, / ')
        print('                   /_/        /___/  ')
        print('v%s' % dftbephy.__version__)
        print('')

        print('-- working with %i processes' % (world_size))
        sys.stdout.flush()

    try:
        hsdinput = hsd.load("dftbephy_in.hsd")
    except:
        if rank==0:
            print('ERROR: You have to provide an input file dftbephy_in.hsd')
        return

    inp_dict = check_hsd_input(hsdinput, 'DFTBephy')
    assert(type(inp_dict) is dict)

    basedir = inp_dict.get('base_dir', './')

    phonopy_dir = inp_dict.get('phonopy_dir', 'phonons/')
    working_dir = inp_dict.get('working_dir', 'el-ph/')
    results_dir = inp_dict.get('results_dir', working_dir)

    name = inp_dict.get('name', '')

    # read section for DFTB
    dftb_dict = inp_dict.get('DFTB', {})
    angular_momenta = dftb_dict.get('angularmomenta', {})
    if len(angular_momenta) == 0:
        if rank==0:
            print('ERROR: angular momenta per element have to be specified.')
        return
    dftb_cmd = dftb_dict.get('cmd', None)
    if not ((dftb_cmd is None) or (isinstance(dftb_cmd, str))):
        if rank==0:
            print('ERROR: invalid cmd.')
        return
    
    # read section for Phonopy
    phonopy_dict = inp_dict.get('Phonopy', {})
    phonopy_yaml    = phonopy_dict.get('yaml_file', 'phonopy_disp.yaml')
    phonopy_symprec = phonopy_dict.get('symprec', None)

    # 1 Load phonopy
    if rank == 0:
        print('-- looking for phonopy results in %s' % (basedir + phonopy_dir))
    os.chdir(basedir + phonopy_dir)

    if phonopy_symprec is None:
        ph = phonopy.load(phonopy_yaml)
        if rank == 0:
            print(" - using phonopy's default symmetry precision.")
    else:
        ph = phonopy.load(phonopy_yaml, symprec=phonopy_symprec)
        if rank == 0:
            print(" - phonopy's symmetry precision set to", phonopy_symprec)

    # 2 Initialize DFTB calculator
    if rank == 0:
        print('-- working in %s' % (basedir + working_dir))
    os.chdir(basedir + working_dir)

    if dftb_cmd is None:
        dftb = DftbSuperCellCalc(angular_momenta)
    else:
        dftb = DftbSuperCellCalc(angular_momenta, cmd=dftb_cmd)
        
    dftb.load_phonopy(ph)

    if os.path.isfile('reference.npz'):
        print('-- [%i] loading reference calculation' % (rank))
        npzfile = np.load('reference.npz')
        dftb.H0 = npzfile['H0']
        dftb.S0 = npzfile['S0']
    else:
        print('-- run reference calculation in serial!' )
        return
    sys.stdout.flush()


    if args.command == 'relaxationtimes':
        run_calc_relaxationtimes(ph, dftb, inp_dict, basedir, phonopy_dir, working_dir, results_dir)

    elif args.command == 'mobility':
        run_calc_transport(ph, dftb, inp_dict, basedir, phonopy_dir, working_dir, results_dir)


def run_calc_relaxationtimes(ph, dftb, inp_dict, basedir, phonopy_dir, working_dir, results_dir):
    from phonopy.structure.grid_points import GridPoints
    # for writing hdf files
    import h5py
    
    comm = MPI.COMM_WORLD
    world_size = comm.Get_size()
    rank = comm.Get_rank()
    start_time = MPI.Wtime()
    
    name = inp_dict.get('name', '')

    # read section for relaxation time calculations
    rt_dict = check_hsd_input(inp_dict, 'RelaxationTimes')

    if 'SERTA'in rt_dict.keys():
        rta_method = 'SERTA'
        rt_dict = rt_dict['SERTA']
    else:
        if rank==0:
            print('ERROR: unsupported method for RelaxationTimes.')
        return

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
        mu_list = np.linspace(float(m[0]), float(m[1]), int(m[2]))
    elif type(rt_dict['mu']) is list:
        mu_list = rt_dict['mu']
    elif type(rt_dict['mu']) is float:
        mu_list = [rt_dict['mu']]
    else:
        if rank==0:
            print('ERROR: unkown type for mu (neither Range, list or float).')
        return

    kBT0 = rt_dict.get('temperature', 0.0259)
    assert(type(kBT0) is float)

    sigma0 = rt_dict.get('sigma', 0.003)
    assert(type(sigma0) is float)

    EF = rt_dict.get('Efermi', 0.00)
    assert(type(EF) is float)

    # 4 Run calculations of electron-phonon couplings

    if os.path.isfile('derivatives.npz'):
        print('-- [%i] loading derivatives calculation' % (rank))
        npzfile = np.load('derivatives.npz')
        dftb.H_derivs = npzfile['H_derivs']
        dftb.S_derivs = npzfile['S_derivs']
    else:
        print('-- run derivatives calculation in serial!' )
        return
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

    mus = EF + np.array(mu_list) # chemical potentials
    sigmas = sigma0* np.ones_like(mus) # smearing
    kBTs = kBT0 * np.ones_like(mus)   # temperatures

    # all calculations
    temps_chems = zip(kBTs, mus)
    ncalcs = len(list(temps_chems))
    if rank == 0:
        print('-- running %d calculations for mu, kBT, sigma:' % (ncalcs))
        for ic, (kBT, mu, sigma_0) in enumerate(zip(kBTs, mus, sigmas)):
            print(' - run %d with mu = %4.3f eV, kBT = %4.3f eV, and sigma = %4.3f' % (ic, mu, kBT, sigma_0))
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

            ds = struct_grp.create_dataset('positions', data=dftb.primitive.positions* BOHR__AA)
            ds = struct_grp.create_dataset('numbers', data=dftb.primitive.numbers)
            ds = struct_grp.create_dataset('cell', data=dftb.primitive.cell* BOHR__AA)
        
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


def run_calc_transport(ph, dftb, inp_dict, basedir, phonopy_dir, working_dir, results_dir):
    # this is needed for writing json files
    import json
    # for irreducible mesh
    import spglib

    comm = MPI.COMM_WORLD
    world_size = comm.Get_size()
    rank = comm.Get_rank()
    start_time = MPI.Wtime()

    name = inp_dict.get('name', '')

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
            print('ERROR: unsupported method for RelaxationTimes.')
        return

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
            print('ERROR: unkown type for mu (neither Range, list or float).')
        return

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
            return
        sys.stdout.flush()

        print('-- [%i] starting phonon calculations on mesh ...' % (rank))
        start = MPI.Wtime()
        ph.run_mesh(q_mesh, with_eigenvectors=False, is_mesh_symmetry=False, is_gamma_center=False)
        mesh = ph.get_mesh_dict()
        end = MPI.Wtime()
        print(' - [%i] finished (%4.1f s).' % (rank, end-start))
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
        print(' - [%i] finished (%4.1f s).' % (rank, end-start))
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
            print(' - run %d with mu = %4.3f eV and kBT = %4.3f eV' % (ic, mu, kBT))
        sys.stdout.flush()

        print(' - constructing k-mesh')
        cell = (dftb.primitive.cell * BOHR__AA, \
           dftb.primitive.scaled_positions, \
           dftb.primitive.numbers)

        fromspglib = spglib.get_ir_reciprocal_mesh(k_mesh, cell)
    
        indices, weights = np.unique(fromspglib[0], return_counts=True)
        weights = np.asarray(weights, dtype='int')    
        kpoints = fromspglib[1]
        kpoints = kpoints[indices,:] / np.array(k_mesh)
    
        nkpoints = len(kpoints)
    
        print(' - number of irreducible k-points: %i' % (nkpoints))
    else:
        nkpoints = None
    
        print(' - [%i] waiting for Godot ...' % (rank))
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

    print(' - [%i] %d epc calculations done.' % (rank, nepccalcs))
    sys.stdout.flush()

    comm.Gatherv(loc_energies, [energies, num_ks_per_rank * nbands, energies_displ_per_rank, MPI.DOUBLE], root=0)
    comm.Gatherv(loc_velocities, [velocities, num_ks_per_rank * nbands * 3, velocities_displ_per_rank, MPI.DOUBLE], root=0)

    comm.Reduce( [loc_densities, MPI.DOUBLE], [densities, MPI.DOUBLE], op = MPI.SUM, root = 0 )
    comm.Reduce( [loc_densities0, MPI.DOUBLE], [densities0, MPI.DOUBLE], op = MPI.SUM, root = 0 )
    comm.Reduce( [loc_conductivities, MPI.DOUBLE], [conductivities, MPI.DOUBLE], op = MPI.SUM, root = 0 )
    comm.Reduce( [loc_conductivities0, MPI.DOUBLE], [conductivities0, MPI.DOUBLE], op = MPI.SUM, root = 0 )

    if rank==0:
        cell = ph.primitive.cell * BOHR__AA
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

if __name__ == '__main__':
    main(arguments=None)

