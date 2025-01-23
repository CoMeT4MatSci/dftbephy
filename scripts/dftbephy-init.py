import sys, os

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
# 4 Run preparation calculations for electron-phonon couplings

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

###############################################################################
print('-- all set for running dftbephy!')
