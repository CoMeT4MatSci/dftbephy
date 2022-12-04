import numpy as np
from scipy import linalg

import multiprocessing as mp
#from multiprocessing import shared_memory

from phonopy.structure.cells import get_smallest_vectors

from .units import *
from .dftb import std_orbital_order

print('*** USING NEW EPC EQUATIONS ***')

try:
    from .extensions import calculate_lattice_ft, calculate_lattice_ft_derivative, calc_g_loc_new, calc_g2_inner_new
except ModuleNotFoundError:
    print('Warning: C-extensions not available')
    from .fourier import calculate_lattice_ft, calculate_lattice_ft_derivative

    def calc_g_loc_new(g_H_loc, g_S_k_loc, g_S_kq_loc, dHdR_k, dSdR_k, dHdR_kq, dSdR_kq, ph_ev, uc2idx, uc_masses, pos_qvec):

        nuc = ph_ev.shape[0]
        g_H_loc.fill(0. + 0j)
        g_S_k_loc.fill(0. + 0j)
        g_S_kq_loc.fill(0. + 0j)

        for s in range(nuc): # atoms in unit cell
            m_s = np.sqrt(uc_masses[s])
            phase_s = np.exp(2j*np.pi*pos_qvec[s])/m_s

            for sp in range(nuc): # atoms in unit cell
                m_sp = np.sqrt(uc_masses[sp])
                phase_sp = np.exp(2j*np.pi*pos_qvec[sp])/m_sp

                for alpha in range(3): # cartesian coordinates
                    g_H_loc[uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]] += phase_s*ph_ev[s,alpha] * dHdR_k[alpha,uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]] - phase_sp*ph_ev[sp,alpha] * dHdR_kq[alpha,uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]]
                    g_S_k_loc[uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]] += phase_s*ph_ev[s,alpha] * dSdR_k[alpha,uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]] 
                    g_S_kq_loc[uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]] += phase_sp*ph_ev[sp,alpha] * dSdR_kq[alpha,uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]]

    
    def calc_g2_inner_new(mesh_g2, iq, lam, gH, gS_k, gS_kq, eps_k, eps_kq, prefactor, band0, band1):
        nbands = mesh_g2.shape[2]
        for j, jband in enumerate(range(band0,band1)):
            for k, kband in enumerate(range(band0,band1)):
                g2 = prefactor * np.abs( gH[jband,kband] - (eps_k[kband]*gS_k[jband,kband] - eps_kq[jband]*gS_kq[jband,kband]) )**2
                mesh_g2[iq, lam, j, k] = g2


def calc_g2(mesh_g2, iq, lam, g_H_loc, g_S_k_loc, g_S_kq_loc, eps_k, eps_kq, U_k, U_kq_H, prefactor, band0, band1):
    gH = U_kq_H @ g_H_loc @ U_k
    gS_k = U_kq_H @ g_S_k_loc @ U_k
    gS_kq = U_kq_H @ g_S_kq_loc @ U_k

    calc_g2_inner_new(mesh_g2, iq, lam, gH, gS_k, gS_kq, eps_k, eps_kq, prefactor, band0, band1)



def calculate_g2(kvec0, band_sel, mesh_qpoints, mesh_frequencies, mesh_eigenvectors, uc2sc, sc2uc, sc2c, supercell, primitive, angular_momenta, ham0, S0, ham_derivs, ovr_derivs):
    """Calculate electron-phonon coupling matrix g_{mn}^lambda(kvec0, q) for the q-points in mesh_qpoints. 
       The respective phonon frequencies and eigenvectors are provided in mesh_frequencies, mesh_eigenvectors.
       Information about the supercell and the primitive cell are obtained from phonopy's supercell and primitive.
       
       returns: eps(kvec) = electronic band-energies at kvec (shape: (nbands,))
                eps(kvec + q) = electronic band-energies at kvec + qvec (shape: (nqpoints, nbands))
                mesh_g2 = absolute square of g_{mn}^lambda(kvec0, q) (shape: (nqpoints, nmodes, nbands, nbands))
    """
    orbitals = [sum([std_orbital_order[am] for am in angular_momenta[cs]], []) for cs in supercell.get_chemical_symbols()]
    norbitals = [len(o) for o in orbitals]
    sc2idx = np.insert(np.cumsum(norbitals), 0, 0) # map of atom in supercell to first orbital in hamiltonian/overlap

    uc_orbitals = [orbitals[i] for i in uc2sc]
    norbitals = [len(o) for o in uc_orbitals]
    uc2idx = np.insert(np.cumsum(norbitals), 0, 0) # map of atom in unitcell to first orbital in hamiltonian/overlap

    uc_pos = primitive.scaled_positions - primitive.scaled_positions[0]
    svecs, multi = get_smallest_vectors(supercell.get_cell(), supercell.scaled_positions, supercell.scaled_positions[uc2sc] )
    trans_mat_float = np.dot(supercell.get_cell(), np.linalg.inv(primitive.get_cell()))
    trans_mat = np.rint(trans_mat_float).astype(int)
    assert (np.abs(trans_mat_float - trans_mat) < 1e-8).all()
    svecs = np.array(np.dot(svecs, trans_mat), dtype="double", order="C")

    # get matrices for initial state at k
    h0_uc = calculate_lattice_ft(ham0, kvec0, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
    s0_uc = calculate_lattice_ft(S0, kvec0, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
    eps_k, U_k = linalg.eigh(h0_uc, b=s0_uc) # diagonalize Hamiltonian
    for iband in range(U_k.shape[1]): # fix gauge by choosing first element of each vector to be real
        U_k[:,iband] = U_k[:,iband] * np.exp( -1j * np.angle(U_k[0,iband]) )
    dHdR_k = calculate_lattice_ft_derivative(ham_derivs, kvec0, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
    dSdR_k = calculate_lattice_ft_derivative(ovr_derivs, kvec0, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)

    uc_masses = supercell.get_masses()[uc2sc] # masses in primitive cell
    Nsc = supercell.get_supercell_matrix().diagonal().prod() # number of primitive cells in super cell

    nqpoints = mesh_frequencies.shape[0] # no of qpoints
    nmodes = mesh_frequencies.shape[1] # no of phonon branches
    nbands = eps_k.shape[0] # no of electronic bands
    
    if (type(band_sel) is list) and (len(band_sel) == 2):
        band0 = max(band_sel[0], 0)
        band1 = min(band_sel[1], nbands)
        nbands = band1 - band0
    else:
        band0 = 0
        band1 = nbands

    # temporary storage for el-ph coupling
    g_H_loc = np.zeros((sum(norbitals), sum(norbitals)), complex)
    g_S_k_loc = np.zeros((sum(norbitals), sum(norbitals)), complex)
    g_S_kq_loc = np.zeros((sum(norbitals), sum(norbitals)), complex)

    # output array for g2 and eps_kq
    mesh_g2 = np.zeros((nqpoints, nmodes, nbands, nbands), float)
    mesh_epskq = np.zeros((nqpoints, nbands), float)

    # print('--[calculate_g2] start iteration')
    for iq, qvec in enumerate(mesh_qpoints):

        # get matrices for state at k + q
        h0_uc = calculate_lattice_ft(ham0, kvec0 + qvec, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
        s0_uc = calculate_lattice_ft(S0, kvec0 + qvec, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
        eps_kq, U_kq = linalg.eigh(h0_uc, b=s0_uc) # diagonalize Hamiltonian
        for iband in range(U_kq.shape[1]): # fix gauge by choosing first element of each vector to be real
            U_kq[:,iband] = U_kq[:,iband] * np.exp( -1j * np.angle(U_kq[0,iband]) )        
        dHdR_kq = calculate_lattice_ft_derivative(ham_derivs, kvec0 + qvec, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
        dSdR_kq = calculate_lattice_ft_derivative(ovr_derivs, kvec0 + qvec, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
        
        mesh_epskq[iq, :] = eps_kq[band0:band1]
        pos_qvec = uc_pos @ qvec

        # calculate g for each phonon branch
        for lam in range(nmodes):
            ph_om = mesh_frequencies[iq][lam]
            ph_ev = mesh_eigenvectors[iq][:,lam].reshape((-1,3))

            # get partial el-ph coupling matrices in local basis
            calc_g_loc_new(g_H_loc, g_S_k_loc, g_S_kq_loc, dHdR_k, dSdR_k, dHdR_kq, dSdR_kq, ph_ev, uc2idx, uc_masses, pos_qvec)

            # transform local el-ph coupling matrix to eigenbasis
            # note: the prefactor does not contain Nsc
            calc_g2(mesh_g2, iq, lam, g_H_loc, g_S_k_loc, g_S_kq_loc, eps_k, eps_kq, U_k, U_kq.conj().T, hbar_amu_THz/(np.maximum(2*ph_om, 1e-5)), band0, band1)
                    
    return eps_k, mesh_epskq, mesh_g2