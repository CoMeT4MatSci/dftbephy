import os
import argparse
import numpy as np

import phonopy

from timeit import default_timer as timer

import dftbephy
from dftbephy import DftbSuperCellCalc
from dftbephy.units import *
from dftbephy.tools import printProgressBar

import hsd

def check_hsd_input(inp_dict, name):
    ret_dict = None
    if name in inp_dict.keys():
        ret_dict = inp_dict[name]
    else:
        print('ERROR: %s not found in input file.')
    return ret_dict

def convert(x):
    if hasattr(x, "tolist"):  # numpy arrays have this
        return x.tolist()
    raise TypeError(x)

def path_dict_to_paths(path_dict):
    path_nodes = path_dict.get('path', [])
    labels = path_dict.get('labels', [])
    npoints = path_dict.get('npoints', 21)  # 21 is the default

    ############################################
    # to get the path as in the previous version
    path = [
        [path_nodes[i], path_nodes[i + 1]]
        for i in range(len(path_nodes) - 1)
    ]
    path_labels = [
        [labels[i], labels[i + 1]]
        for i in range(len(labels) - 1)
    ]
    ############################################
    
    return path, path_labels, npoints

def reorder_modes(ph_oms, ph_evs):
    ######## reorder phonon modes ###########
    ## based on the idea presented in https://quantumtinkerer.tudelft.nl/blog/connecting-the-dots/
    ##
    from scipy.optimize import linear_sum_assignment

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

    return sorted_freqs, sorted_evs

def cart_vecs(points, cell):
    a0 = np.linalg.norm(cell[0,:])

    reciprocal_lattice = np.linalg.inv(cell)
    # cartesian vectors in units of 2*np.pi/a0
    return a0*np.vstack(points) @reciprocal_lattice.T 

# commands and descritions for argparse
commands = [('init',  'Prepare calculations.'),
            ('bands', 'Calculate electronic band-structure and phonon dispersions along a path.'),
            ('epc',   'Calculate electron-phonon coupling matrix on a mesh.'),
            ('ephline', 'Calculate electron-phonon coupling matrix along a qp path.'),
            ('rtline', 'Calculate electronic relaxation-times along a kp path.')]
            
#            ('relaxationtimes', 'Calculate relaxation times.'),
#            ('mobility'), 'Calculate transport properties.']

def main(arguments=None):
    ############# PARSE COMMAND LINE ARGUMENTS #############
    parser = argparse.ArgumentParser(
                        prog='dftbephy',
                        description='Calculating electron-phonon couplings (EPCs) with DFTB.')

    parser.add_argument('-v', '--verbose', action='store_true')  # on/off flag

    subparsers = parser.add_subparsers(title='Commands', dest='command')

    for (cmd, desc) in commands:
        p = subparsers.add_parser(cmd, description=desc)

    args = parser.parse_args()


    ############# START #############
    #
    print('     _______  __            __       ')
    print(' ___/ / _/ /_/ /  ___ ___  / /  __ __')
    print('/ _  / _/ __/ _ \\/ -_) _ \\/ _ \\/ // /')
    print('\\_,_/_/ \\__/_.__/\\__/ .__/_//_/\\_, / ')
    print('                   /_/        /___/  ')
    print('v%s' % dftbephy.__version__)
    print('')

    try:
        hsdinput = hsd.load("dftbephy_in.hsd")
    except:
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
        print('ERROR: angular momenta per element have to be specified.')
        return
    dftb_cmd = dftb_dict.get('cmd', None)
    if not ((dftb_cmd is None) or (isinstance(dftb_cmd, str))):
        print('ERROR: invalid cmd.')
        return
        
    
    # read section for Phonopy
    phonopy_dict = inp_dict.get('Phonopy', {})
    phonopy_yaml    = phonopy_dict.get('yaml_file', 'phonopy_disp.yaml')
    phonopy_symprec = phonopy_dict.get('symprec', None)

    # 1 Load phonopy
    print('-- looking for phonopy results in %s' % (basedir + phonopy_dir))
    os.chdir(basedir + phonopy_dir)

    if phonopy_symprec is None:
        ph = phonopy.load(phonopy_yaml)
        print(" - using phonopy's default symmetry precision.")
    else:
        ph = phonopy.load(phonopy_yaml, symprec=phonopy_symprec)
        print(" - phonopy's symmetry precision set to", phonopy_symprec)


    # 2 Initialize DFTB calculator
    dftb = run_init(ph, angular_momenta, dftb_cmd, basedir, phonopy_dir, working_dir, results_dir)

    if args.command == 'bands':
        run_calc_el_bands(ph, dftb, inp_dict, basedir, phonopy_dir, working_dir, results_dir)
        run_calc_ph_bands(ph, dftb, inp_dict, basedir, phonopy_dir, working_dir, results_dir)

    elif args.command == 'epc':
        run_calc_epc(ph, dftb, inp_dict, basedir, phonopy_dir, working_dir, results_dir)

    elif args.command == 'ephline':
        run_calc_ephline(ph, dftb, inp_dict, basedir, phonopy_dir, working_dir, results_dir)

    elif args.command == 'rtline':
        run_calc_rtline(ph, dftb, inp_dict, basedir, phonopy_dir, working_dir, results_dir)


def run_init(ph, angular_momenta, dftb_cmd, basedir, phonopy_dir, working_dir, results_dir):
    '''
        Initialize DFTB calculator and run reference calculations
    '''
    print('-- working in %s' % (basedir + working_dir))
    os.chdir(basedir + working_dir)

    if dftb_cmd is None:
        dftb = DftbSuperCellCalc(angular_momenta)
    else:
        dftb = DftbSuperCellCalc(angular_momenta, cmd=dftb_cmd)
        
    dftb.load_phonopy(ph)

    # Store Hamiltonian and overlap matrices
    if os.path.isfile('reference.npz'):
        print(' - loading reference calculation')
        npzfile = np.load('reference.npz')
        dftb.H0 = npzfile['H0']
        dftb.S0 = npzfile['S0']
    else:
        print(' - starting reference calculation ...')
        start = timer()
        dftb.calculate_reference()
        end = timer()
        print(' - finished (%4.1f s).' % (end-start))

        np.savez('reference.npz', H0=dftb.H0, S0=dftb.S0)

    # Run preparation calculations for electron-phonon couplings
    # and store gradients
    if os.path.isfile('derivatives.npz'):
        print(' - loading derivatives calculation')
        npzfile = np.load('derivatives.npz')
        dftb.H_derivs = npzfile['H_derivs']
        dftb.S_derivs = npzfile['S_derivs']
    else:
        print(' - starting derivatives calculation ...')
        start = timer()
        dftb.calculate_derivatives()
        end = timer()
        print(' - finished (%4.1f s).' % (end-start))

        np.savez('derivatives.npz', H_derivs=dftb.H_derivs, S_derivs=dftb.S_derivs)

    print('-- all set for running dftbephy!')
    
    return dftb


def run_calc_el_bands(ph, dftb, inp_dict, basedir, phonopy_dir, working_dir, results_dir):
    import json
    from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
    
    # read section for band calculations
    bands_dict = check_hsd_input(inp_dict, 'Bands')
    path, path_labels, npoints = path_dict_to_paths(bands_dict)

    # get k-points
    kpoints, connections = get_band_qpoints_and_path_connections(path, npoints=npoints)

    # call dftbephy to use get_num_bands
    npaths = len(kpoints)
    nbands = dftb.get_num_bands()

    # rund electronic band-structure calculation
    print('-- starting electronic calculations on paths ...')
    bands = []
    vels = []
    for i in range(npaths):
        print(' - processing path %s to %s' % (path[i][0], path[i][1]))
        nkpoints = kpoints[i].shape[0]
        energies = np.zeros((nkpoints, nbands), float)
        velocities = np.zeros((nkpoints, nbands, 3), float)

        printProgressBar(0, nkpoints, prefix='k-point', suffix='complete')
        for ik, kvec in enumerate(kpoints[i]):
            eps_k, vel_k = dftb.calculate_velocity(kvec)
                
            for n in range(nbands):
                energies[ik,n] = eps_k[n]
                velocities[ik,n,:] = vel_k[:,n]

            printProgressBar(ik+1, nkpoints, prefix='k-point', suffix='complete')

        bands.append(energies)
        vels.append(velocities)

    # convert k-points to cartesian coordinates
    primitive_cell = dftb.primitive.cell * BOHR__AA
    kvecs = cart_vecs(kpoints, primitive_cell)  # in units of 2*np.pi/a0
    
    # generate output as json
    bs_dict = { 'particleType': 'electron', 'numBands': bands[0].shape[1], 'energies': np.vstack(bands), 'energyUnit': 'eV', 
      'coordsType': 'lattice', 'highSymCoordinates': np.vstack(path), 'highSymLabels': np.array(path_labels).flatten(), 
      'highSymIndices': np.cumsum(np.array([[0, kps.shape[0]] for kps in kpoints]).flatten()),
      'wavevectorCoordinates': np.vstack(kpoints), 
      'wavevectorCoordinatesCart': kvecs,
      'wavevectorIndices': list(range(np.vstack(kpoints).shape[0]))}

    with open('path_el_bandstructure.json', 'w') as outfile:
        json.dump(json.dumps(bs_dict, default=convert), outfile)


def run_calc_ph_bands(ph, dftb, inp_dict, basedir, phonopy_dir, working_dir, results_dir):
    import json
    from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
    
    # read section for band calculations
    bands_dict = check_hsd_input(inp_dict, 'Bands')
    path, path_labels, npoints = path_dict_to_paths(bands_dict)

    # get q-points
    qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=npoints)
    npaths = len(qpoints)

    print('-- starting phonon calculations on paths ...')
    start = timer()
    ph.run_band_structure(qpoints, is_band_connection=False, path_connections=connections, with_eigenvectors=True)
    bs_dict = ph.get_band_structure_dict()
    end = timer()
    print(' - finished (%4.1f s).' % (end-start))

    # reorder phonon modes
    ph_evs = bs_dict['eigenvectors']
    ph_oms = bs_dict['frequencies']
    qpoints= bs_dict['qpoints']
    frequencies, eigenvectors = reorder_modes(ph_oms, ph_evs)

    # convert q-points to cartesian coordinates
    primitive_cell = dftb.primitive.cell * BOHR__AA
    qvecs = cart_vecs(qpoints, primitive_cell)  # in units of 2*np.pi/a0
    
    # generate output as json
    bs_dict = { 'particleType': 'phonon', 'numModes': frequencies[0].shape[1], 'frequencies': np.vstack(frequencies)*THZ__EV, 'frequencyUnit': 'eV',
      'coordsType': 'lattice', 'highSymCoordinates': np.vstack(path), 'highSymLabels': np.array(path_labels).flatten(),
      'highSymIndices': np.cumsum(np.array([[0, qps.shape[0]] for qps in qpoints]).flatten()),
      'wavevectorCoordinates': np.vstack(qpoints),
      'wavevectorCoordinatesCart': qvecs,
      'wavevectorIndices': list(range(np.vstack(qpoints).shape[0]))}

    with open('path_ph_bandstructure.json', 'w') as outfile:
        json.dump(json.dumps(bs_dict, default=convert), outfile)
    

def run_calc_epc(ph, dftb, inp_dict, basedir, phonopy_dir, working_dir, results_dir):
    import h5py
    
    # read section for epc calculations from input file
    epc_dict = check_hsd_input(inp_dict, 'EPCs')

    default_mesh = {'Mesh': {'npoints': [1, 1, 1], 'refinement': 1, 'shift': [0.0,0.0,0.0]}}
    qp_dict = epc_dict.get('qpoints', default_mesh)
    if 'Mesh' in qp_dict.keys():
        # size of q-point mesh (eg. nq x nq x 1)
        q_mesh = qp_dict['Mesh'].get('npoints', default_mesh['Mesh']['npoints'])
        # factor by which the q-points are scaled
        q_mesh_refinement = qp_dict['Mesh'].get('refinement', default_mesh['Mesh']['refinement'])
        q_mesh_shift = qp_dict['Mesh'].get('shift', default_mesh['Mesh']['shift'])
    else:
        print('ERROR: no q-point mesh specified in EPC section of the input file')
        return

    k_mesh_shift = epc_dict.get('kvec0', [0., 0., 0.])
    kvec0 = np.array(k_mesh_shift) # reference k-point for electrons

    band_sel = epc_dict.get('bands', None)
    if type(band_sel) is list:
        band_sel[1] = band_sel[1]+1

    calculate_velocities = epc_dict.get('velocities', False)
    
    print('-- starting phonon calculations on mesh ...')
    start = timer()
    ph.run_mesh(q_mesh, with_eigenvectors=False, is_mesh_symmetry=False, is_gamma_center=False)
    mesh = ph.get_mesh_dict()
    end = timer()
    print(' - finished (%4.1f s).' % (end-start))

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
    print(' - finished (%4.1f s).' % (end-start))

    # convert q-points to cartesian coordinates
    primitive_cell = dftb.primitive.cell * BOHR__AA
    qvecs = cart_vecs(mesh_qpoints, primitive_cell)  # in units of 2*np.pi/a0

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
        print(' - finished (%4.1f s).' % (end-start))

    # store coupling matrix
    nkp = 0
    ik = 0
    with h5py.File('el-ph-Nq%i-K-bandsel.hdf5' % (q_mesh[0]), 'w') as f:
        ph_grp = f.create_group('ph')
        ph_grp.attrs['mesh'] = ph._mesh.mesh_numbers
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


def run_calc_ephline(ph, dftb, inp_dict, basedir, phonopy_dir, working_dir, results_dir):
    import json
    from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections

    # read section for epc calculations
    epc_dict = check_hsd_input(inp_dict, 'EPCs')

    default_mesh = {'Mesh': {'npoints': [1, 1, 1], 'refinement': 1, 'shift': [0.0,0.0,0.0]}}
    qp_dict = epc_dict.get('qpoints', default_mesh)
    if 'Path' in qp_dict.keys():
        # read section for band-path for ephline
        path, path_labels, npoints = path_dict_to_paths(qp_dict['Path'])
    else:
        print('ERROR: no band-path specified in EPC section of the input file')
        return

    k_mesh_shift = epc_dict.get('kvec0', [0., 0., 0.])
    kvec0 = np.array(k_mesh_shift) # reference k-point for electrons

    # selection of a subset of bands
    band_sel = epc_dict.get('bands', None)
    if type(band_sel) is list:
        band_sel[1] = band_sel[1]+1

    # TODO: currently ignored for ephline
    calculate_velocities = epc_dict.get('velocities', False)
  
    qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=npoints)
    npaths = len(qpoints)   # number of paths

    print('-- starting phonon calculations on path ...')
    start = timer()
    ph.run_band_structure(qpoints, is_band_connection=False, path_connections=connections, with_eigenvectors=True)
    bs_dict = ph.get_band_structure_dict()
    end = timer()
    print(' - finished (%4.1f s).' % (end-start))

    # reorder phonon modes
    ph_evs = bs_dict['eigenvectors']
    ph_oms = bs_dict['frequencies']
    qpoints= bs_dict['qpoints']
    frequencies, eigenvectors = reorder_modes(ph_oms, ph_evs)

    print('-- starting el-ph calculation ...')
    bands = []
    g_kq = []
    for i in range(npaths):
        print(' - processing path %s to %s' % (path[i][0], path[i][1]))
        nqpoints = qpoints[i].shape[0]

        eps_k, epskq, mesh_epskmq, g2 = dftb.calculate_g2(kvec0, qpoints[i], frequencies[i], eigenvectors[i], band_sel=band_sel)

        bands.append(epskq)
        g_kq.append(np.sqrt(g2))

    # convert q-points to cartesian coordinates
    primitive_cell = dftb.primitive.cell * BOHR__AA
    qvecs = cart_vecs(qpoints, primitive_cell)  # in units of 2*np.pi/a0

    # generate output as json
    ephl_dict = {'numBands': bands[0].shape[1], 'energies': np.vstack(bands), 'energyUnit': 'eV',
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


def run_calc_rtline(ph, dftb, inp_dict, basedir, phonopy_dir, working_dir, results_dir):
    import json
    from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
    from dftbephy.analysis import inv_tau_nk_lam

    # read section for relaxation time calculations
    rt_dict = check_hsd_input(inp_dict, 'RelaxationTimes')

    if 'SERTA'in rt_dict.keys():
        rta_method = 'SERTA'
        rt_dict = rt_dict['SERTA']
    else:
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
    if 'Path' in kp_dict.keys():
        # read section for band-path for ephline
        path, path_labels, npoints = path_dict_to_paths(qp_dict['Path'])
    else:
        print('ERROR: no band-path specified in RelaxationTimes section of the input file')
        return

    k_mesh_shift = rt_dict.get('kvec0', [0., 0., 0.])
    kvec0 = np.array(k_mesh_shift) # reference k-point for electrons

    # selection of a subset of bands
    band_sel = rt_dict.get('bands', None)
    if type(band_sel) is list:
        band_sel[1] = band_sel[1]+1

    # chemical potential(s)
    if type(rt_dict['mu']) is dict:
        m = rt_dict['mu']['Range']
        mu_list = np.linspace(float(m[0]), float(m[1]), int(m[2]))
    elif type(rt_dict['mu']) is list:
        mu_list = rt_dict['mu']
    elif type(rt_dict['mu']) is float:
        mu_list = [rt_dict['mu']]
    else:
        print('ERROR: unkown type for mu (neither Range, list or float).')
        return

    kBT0 = rt_dict.get('temperature', 0.0259)
    assert(type(kBT0) is float)

    sigma0 = rt_dict.get('sigma', 0.003)
    assert(type(sigma0) is float)

    EF = rt_dict.get('Efermi', 0.00)
    assert(type(EF) is float)


    print('-- starting phonon calculations on mesh ...')
    start = timer()
    ph.run_mesh(q_mesh, with_eigenvectors=False, is_mesh_symmetry=False, is_gamma_center=False)
    mesh = ph.get_mesh_dict()
    end = timer()
    print(' - finished (%4.1f s).' % (end-start))

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
    print(' - finished (%4.1f s).' % (end-start))


    # get k-points along paths
    kpoints, connections = get_band_qpoints_and_path_connections(path, npoints=npoints)
    npaths = len(kpoints)   # number of paths
    nbands = mesh_g2.shape[2]

    print('-- starting relaxation-times calculation ...')
    bands = []
    linewidths = []
    for i in range(npaths):
        print(' - processing path %s to %s' % (path[i][0], path[i][1]))
        nkpoints = kpoints[i].shape[0]
        energies = np.zeros((nkpoints, nbands), float)
        taus = np.zeros((nkpoints, nbands), float)

        printProgressBar(0, nkpoints, prefix='k-point', suffix='complete')
        for ik, kvec in enumerate(kpoints[i]):
            eps_k, mesh_epskq, mesh_epskmq, mesh_g2 = dftb.calculate_g2(kvec, mesh_qpoints, mesh_frequencies, mesh_eigenvectors, band_sel=band_sel)
        
            for n in range(nbands):
                energies[ik,n] = eps_k[n]
                taus[ik,n] = inv_tau_nk_lam( n, eps_k[n], mu, kBT, mesh_g2, mesh_epskq, mesh_frequencies*THZ__EV, sigma=sigma_0)[0]/q_mesh_refinement**2

            printProgressBar(ik+1, nkpoints, prefix='k-point', suffix='complete')

        bands.append(energies)
        linewidths.append(taus)

    # convert q-points to cartesian coordinates
    primitive_cell = dftb.primitive.cell * BOHR__AA
    kvecs = cart_vecs(kpoints, primitive_cell)  # in units of 2*np.pi/a0

    # generate output as json
    bs_dict = { 'particleType': 'electron', 'numBands': bands[0].shape[1], 'energies': np.vstack(bands), 'energyUnit': 'eV', 
      'coordsType': 'lattice', 'highSymCoordinates': np.vstack(path), 'highSymLabels': np.array(path_labels).flatten(), 
      'highSymIndices': np.cumsum(np.array([[0, kps.shape[0]] for kps in kpoints]).flatten()),
      'wavevectorCoordinates': np.vstack(kpoints), 'wavevectorIndices': list(range(np.vstack(kpoints).shape[0]))}

    with open('path_el_bandstructure.json', 'w') as outfile:
        json.dump(json.dumps(bs_dict, default=convert), outfile)

    # generate output as json
    lw_dict = { 'chemicalPotentialUnit': 'eV', 'chemicalPotentials': [mu,], 'particleType': 'electron', 'numBands': bands[0].shape[1], 
      'energies': np.vstack(bands), 'energyUnit': 'eV', 
      'coordsType': 'lattice', 'linewidths': np.vstack(linewidths), 'linewidthsUnit': 'eV', 'temperatureUnit': 'eV', 'temperatures': [kBT,]}

    with open('path_el_relaxation_times.json', 'w') as outfile:
        json.dump(json.dumps(lw_dict, default=convert), outfile)


if __name__ == '__main__':
    main(arguments=None)

