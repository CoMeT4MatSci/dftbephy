import numpy as np
from scipy import linalg

#import multiprocessing as mp
#from multiprocessing import shared_memory

from phonopy.structure.cells import get_smallest_vectors

from .units import *
from .dftb import std_orbital_order


try:
    from .extensions import calculate_lattice_ft, calculate_lattice_ft_derivative, calculate_lattice_double_ft_derivative, calc_g_loc
except ModuleNotFoundError:
    print('Warning: C-extensions not available')
    from .fourier import calculate_lattice_ft, calculate_lattice_ft_derivative, calculate_lattice_double_ft_derivative

    def calc_g_loc(g_H_loc, g_S_k_loc, g_S_kq_loc, dHdR_k, dSdR_k, dHdR_kq, dSdR_kq, ph_ev, uc2idx, uc_masses, pos_qvec):

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


def calc_g2_inner(mesh_g2, iq, lam, gH, gS_k, gS_kq, eps_k, eps_kq, prefactor, band0, band1):
    nbands = mesh_g2.shape[2]
    for j, jband in enumerate(range(band0,band1)):
        for k, kband in enumerate(range(band0,band1)):
            g = gH[jband,kband] - (eps_k[kband]*gS_k[jband,kband] - eps_kq[jband]*gS_kq[jband,kband])

            if np.isrealobj(mesh_g2):
                mesh_g2[iq, lam, j, k] = prefactor * np.abs( g )**2
            else:
                mesh_g2[iq, lam, j, k] = np.sqrt(prefactor) * g
                    

def calc_g2(mesh_g2, iq, lam, g_H_loc, g_S_k_loc, g_S_kq_loc, eps_k, eps_kq, U_k, U_kq_H, prefactor, band0, band1):
    gH = U_kq_H @ g_H_loc @ U_k
    gS_k = U_kq_H @ g_S_k_loc @ U_k
    gS_kq = U_kq_H @ g_S_kq_loc @ U_k

    calc_g2_inner(mesh_g2, iq, lam, gH, gS_k, gS_kq, eps_k, eps_kq, prefactor, band0, band1)



def calculate_g2(kvec0, band_sel, mesh_qpoints, mesh_frequencies, mesh_eigenvectors, uc2sc, sc2uc, sc2c, supercell, primitive, angular_momenta, ham0, S0, ham_derivs, ovr_derivs):
    """Calculate electron-phonon coupling matrix g_{mn}^lambda(kvec0, q) for the q-points in mesh_qpoints. 
       The respective phonon frequencies and eigenvectors are provided in mesh_frequencies, mesh_eigenvectors.
       Information about the supercell and the primitive cell are obtained from phonopy's supercell and primitive.
       
       returns: eps(kvec) = electronic band-energies at kvec (shape: (nbands,))
                eps(kvec + q) = electronic band-energies at kvec + qvec (shape: (nqpoints, nbands))
                mesh_g2 = absolute square of g_{mn}^lambda(kvec0, q) (shape: (nqpoints, nmodes, nbands, nbands))
    """
    orbitals = [sum([std_orbital_order[am] for am in angular_momenta[cs]], []) for cs in supercell.symbols]
    norbitals = [len(o) for o in orbitals]
    sc2idx = np.insert(np.cumsum(norbitals), 0, 0) # map of atom in supercell to first orbital in hamiltonian/overlap

    uc_orbitals = [orbitals[i] for i in uc2sc]
    norbitals = [len(o) for o in uc_orbitals]
    uc2idx = np.insert(np.cumsum(norbitals), 0, 0) # map of atom in unitcell to first orbital in hamiltonian/overlap

    uc_pos = primitive.scaled_positions - primitive.scaled_positions[0]
    svecs, multi = get_smallest_vectors(supercell.cell, supercell.scaled_positions, supercell.scaled_positions[uc2sc] )
    trans_mat_float = np.dot(supercell.cell, np.linalg.inv(primitive.cell))
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

    uc_masses = supercell.masses[uc2sc] # masses in primitive cell
    Nsc = supercell.supercell_matrix.diagonal().prod() # number of primitive cells in super cell

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

    # output array for g2, eps_kq, eps_kmq
    mesh_g2 = np.zeros((nqpoints, nmodes, nbands, nbands), complex)
    mesh_epskq = np.zeros((nqpoints, nbands), float)
    mesh_epskmq = np.zeros((nqpoints, nbands), float)

    # print('--[calculate_g2] start iteration')
    for iq, qvec in enumerate(mesh_qpoints):

        # get matrices for state at k - q
        h0_uc = calculate_lattice_ft(ham0, kvec0 - qvec, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
        s0_uc = calculate_lattice_ft(S0, kvec0 - qvec, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
        eps_kmq, U_kmq = linalg.eigh(h0_uc, b=s0_uc) # diagonalize Hamiltonian
        mesh_epskmq[iq, :] = eps_kmq[band0:band1]

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
            calc_g_loc(g_H_loc, g_S_k_loc, g_S_kq_loc, dHdR_k, dSdR_k, dHdR_kq, dSdR_kq, ph_ev, uc2idx, uc_masses, pos_qvec)

            # transform local el-ph coupling matrix to eigenbasis
            # note: the prefactor does not contain Nsc
            calc_g2(mesh_g2, iq, lam, g_H_loc, g_S_k_loc, g_S_kq_loc, eps_k, eps_kq, U_k, U_kq.conj().T, hbar_amu_THz/(np.maximum(2*ph_om, 1e-5)), band0, band1)
                    
    return eps_k[band0:band1], mesh_epskq, mesh_epskmq, mesh_g2


def calc_g_loc_scc(g_H_loc, g_S_k_loc, g_S_kq_loc, dHdR_kq, dSdR_k, dSdR_kq, ph_ev, uc2idx, uc_masses, pos_qvec):
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
                g_S_k_loc[uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]]  += phase_s*ph_ev[s,alpha] * dSdR_k[alpha,uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]]
                g_S_kq_loc[uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]] += phase_sp*ph_ev[sp,alpha] * dSdR_kq[alpha,uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]]

            for spp in range(nuc): # atoms in unit cell
                m_spp = np.sqrt(uc_masses[spp])
                phase_spp = np.exp(2j*np.pi*pos_qvec[spp])/m_spp

                for alpha in range(3): # cartesian coordinates
                    g_H_loc[uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]] += phase_spp*ph_ev[spp,alpha] * dHdR_kq[spp,alpha,uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]]


def calculate_scc_g2(kvec0, band_sel, mesh_qpoints, mesh_frequencies, mesh_eigenvectors, uc2sc, sc2uc, sc2c, supercell, primitive, angular_momenta, ham0, S0, ham_derivs, ovr_derivs):
    """Calculate self-consistent electron-phonon coupling matrix g_{mn}^lambda(kvec0, q) for the q-points in mesh_qpoints. 
       The respective phonon frequencies and eigenvectors are provided in mesh_frequencies, mesh_eigenvectors.
       Information about the supercell and the primitive cell are obtained from phonopy's supercell and primitive.
       
       returns: eps(kvec) = electronic band-energies at kvec (shape: (nbands,))
                eps(kvec + q) = electronic band-energies at kvec + qvec (shape: (nqpoints, nbands))
                mesh_g2 = absolute square of g_{mn}^lambda(kvec0, q) (shape: (nqpoints, nmodes, nbands, nbands))
    """
    print('** INFO: using scc EPC calculation v2 **')
    
    orbitals = [sum([std_orbital_order[am] for am in angular_momenta[cs]], []) for cs in supercell.symbols]
    norbitals = [len(o) for o in orbitals]
    sc2idx = np.insert(np.cumsum(norbitals), 0, 0) # map of atom in supercell to first orbital in hamiltonian/overlap

    uc_orbitals = [orbitals[i] for i in uc2sc]
    norbitals = [len(o) for o in uc_orbitals]
    uc2idx = np.insert(np.cumsum(norbitals), 0, 0) # map of atom in unitcell to first orbital in hamiltonian/overlap

    uc_pos = primitive.scaled_positions - primitive.scaled_positions[0]
    svecs, multi = get_smallest_vectors(supercell.cell, supercell.scaled_positions, supercell.scaled_positions[uc2sc] )
    trans_mat_float = np.dot(supercell.cell, np.linalg.inv(primitive.cell))
    trans_mat = np.rint(trans_mat_float).astype(int)
    assert (np.abs(trans_mat_float - trans_mat) < 1e-8).all()
    svecs = np.array(np.dot(svecs, trans_mat), dtype="double", order="C")

    # get matrices for initial state at k
    h0_uc = calculate_lattice_ft(ham0, kvec0, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
    s0_uc = calculate_lattice_ft(S0, kvec0, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
    eps_k, U_k = linalg.eigh(h0_uc, b=s0_uc) # diagonalize Hamiltonian
    for iband in range(U_k.shape[1]): # fix gauge by choosing first element of each vector to be real
        U_k[:,iband] = U_k[:,iband] * np.exp( -1j * np.angle(U_k[0,iband]) )
    # dHdR_k = calculate_lattice_ft_derivative(ham_derivs, kvec0, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
    dSdR_k = calculate_lattice_ft_derivative(ovr_derivs, kvec0, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)

    uc_masses = supercell.masses[uc2sc] # masses in primitive cell
    Nsc = supercell.supercell_matrix.diagonal().prod() # number of primitive cells in super cell

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

    # output array for g2, eps_kq, eps_kmq
    mesh_g2 = np.zeros((nqpoints, nmodes, nbands, nbands), complex)
    mesh_epskq = np.zeros((nqpoints, nbands), float)
    mesh_epskmq = np.zeros((nqpoints, nbands), float)

    # print('--[calculate_g2] start iteration')
    for iq, qvec in enumerate(mesh_qpoints):

        # get matrices for state at k + q
        h0_uc = calculate_lattice_ft(ham0, kvec0 + qvec, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
        s0_uc = calculate_lattice_ft(S0, kvec0 + qvec, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
        eps_kq, U_kq = linalg.eigh(h0_uc, b=s0_uc) # diagonalize Hamiltonian
        for iband in range(U_kq.shape[1]): # fix gauge by choosing first element of each vector to be real
            U_kq[:,iband] = U_kq[:,iband] * np.exp( -1j * np.angle(U_kq[0,iband]) )
        dHdR_kq = calculate_lattice_double_ft_derivative(ham_derivs, kvec0 + qvec, kvec0, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
        dSdR_kq = calculate_lattice_ft_derivative(ovr_derivs, kvec0 + qvec, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
        
        mesh_epskq[iq, :] = eps_kq[band0:band1]
        pos_qvec = uc_pos @ qvec

        # calculate g for each phonon branch
        for lam in range(nmodes):
            ph_om = mesh_frequencies[iq][lam]
            ph_ev = mesh_eigenvectors[iq][:,lam].reshape((-1,3))

            # get partial el-ph coupling matrices in local basis
            calc_g_loc_scc(g_H_loc, g_S_k_loc, g_S_kq_loc, dHdR_kq, dSdR_k,  dSdR_kq, ph_ev, uc2idx, uc_masses, pos_qvec)

            # transform local el-ph coupling matrix to eigenbasis
            # note: the prefactor does not contain Nsc
            calc_g2(mesh_g2, iq, lam, g_H_loc, g_S_k_loc, g_S_kq_loc, eps_k, eps_kq, U_k, U_kq.conj().T, hbar_amu_THz/(np.maximum(2*ph_om, 1e-5)), band0, band1)
                    
    return eps_k[band0:band1], mesh_epskq, mesh_epskmq, mesh_g2


def calc_M_S_loc(at, uc_pos, pos_qvec, s0_k_uc, s0_kq_uc, uc2idx, norbitals):
    nuc = len(pos_qvec)

    z_at = uc_pos[at, 2] # z component of position of atom at
    
    M_kq_loc = np.zeros((sum(norbitals), sum(norbitals)), complex)
    S_kq_loc = np.zeros((sum(norbitals), sum(norbitals)), complex)
    for s in range(nuc): # atoms in unit cell
        phase_s = np.exp(2j*np.pi*pos_qvec[s])
        z_s = uc_pos[s, 2] 

        for sp in range(nuc): # atoms in unit cell
            phase_sp = np.exp(2j*np.pi*pos_qvec[sp])
            z_sp = uc_pos[sp, 2] 

            M_kq_loc[uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]] += phase_s * s0_k_uc[uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]]
            M_kq_loc[uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]] += phase_sp* s0_kq_uc[uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]]
            
            S_kq_loc[uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]] += np.sign(z_s - z_at)  * phase_s * s0_k_uc[uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]]
            S_kq_loc[uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]] += np.sign(z_sp - z_at) * phase_sp* s0_kq_uc[uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]]

    return M_kq_loc, S_kq_loc  


def calculate_long_g2(kvec0, band_sel, mesh_qpoints, mesh_frequencies, mesh_eigenvectors, uc2sc, sc2uc, sc2c, supercell, primitive, angular_momenta, ham0, S0, Zstar, alpha2D):
    """Calculate the long-range electron-phonon coupling matrix g_{mn}^lambda(kvec0, q) for the q-points in mesh_qpoints. 
       The respective phonon frequencies and eigenvectors are provided in mesh_frequencies, mesh_eigenvectors.
       Information about the supercell and the primitive cell are obtained from phonopy's supercell and primitive.
       Additionally, the Born charges (Zstar) and the 2D polarizability tensor (alpha2D) have to be specified.

       The formalism follows from Deng et al, Phys. Rev. B 103, 075410 (2021)
       https://doi.org/10.1103/PhysRevB.103.075410
       
       returns: eps(kvec) = electronic band-energies at kvec (shape: (nbands,))
                eps(kvec + q) = electronic band-energies at kvec + qvec (shape: (nqpoints, nbands))
                mesh_g2 = absolute square of g_{mn}^lambda(kvec0, q) (shape: (nqpoints, nmodes, nbands, nbands))
    """
    print('** INFO: using long-range EPC calculation **')
    
    orbitals = [sum([std_orbital_order[am] for am in angular_momenta[cs]], []) for cs in supercell.symbols]
    norbitals = [len(o) for o in orbitals]
    sc2idx = np.insert(np.cumsum(norbitals), 0, 0) # map of atom in supercell to first orbital in hamiltonian/overlap

    uc_orbitals = [orbitals[i] for i in uc2sc]
    norbitals = [len(o) for o in uc_orbitals]
    uc2idx = np.insert(np.cumsum(norbitals), 0, 0) # map of atom in unitcell to first orbital in hamiltonian/overlap

    uc_pos = primitive.scaled_positions - primitive.scaled_positions[0]
    svecs, multi = get_smallest_vectors(supercell.cell, supercell.scaled_positions, supercell.scaled_positions[uc2sc] )
    trans_mat_float = np.dot(supercell.cell, np.linalg.inv(primitive.cell))
    trans_mat = np.rint(trans_mat_float).astype(int)
    assert (np.abs(trans_mat_float - trans_mat) < 1e-8).all()
    svecs = np.array(np.dot(svecs, trans_mat), dtype="double", order="C")

    # get matrices for initial state at k
    h0_uc = calculate_lattice_ft(ham0, kvec0, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
    s0_k_uc = calculate_lattice_ft(S0, kvec0, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
    eps_k, U_k = linalg.eigh(h0_uc, b=s0_k_uc) # diagonalize Hamiltonian
    for iband in range(U_k.shape[1]): # fix gauge by choosing first element of each vector to be real
        U_k[:,iband] = U_k[:,iband] * np.exp( -1j * np.angle(U_k[0,iband]) )

    uc_masses = supercell.masses[uc2sc] # masses in primitive cell
    Nsc = supercell.supercell_matrix.diagonal().prod() # number of primitive cells in super cell
    nuc = len(uc_masses)

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

    # output array for g2, eps_kq, eps_kmq
    mesh_g2 = np.zeros((nqpoints, nmodes, nbands, nbands), complex)
    mesh_epskq = np.zeros((nqpoints, nbands), float)
    mesh_epskmq = np.zeros((nqpoints, nbands), float)

    # direction along z-axis
    zdir = np.array([0., 0., 1.])

    # print('--[calculate_g2] start iteration')
    for iq, qvec in enumerate(mesh_qpoints):

        qnorm = np.linalg.norm(qvec)
        if qnorm < 1e-6:
            qdir = np.array([0., 0., 0])
        else:
            qdir = qvec/qnorm
            
        # get matrices for state at k + q
        h0_uc = calculate_lattice_ft(ham0, kvec0 + qvec, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
        s0_kq_uc = calculate_lattice_ft(S0, kvec0 + qvec, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi)
        eps_kq, U_kq = linalg.eigh(h0_uc, b=s0_kq_uc) # diagonalize Hamiltonian
        for iband in range(U_kq.shape[1]): # fix gauge by choosing first element of each vector to be real
            U_kq[:,iband] = U_kq[:,iband] * np.exp( -1j * np.angle(U_kq[0,iband]) )
        
        mesh_epskq[iq, :] = eps_kq[band0:band1]
        pos_qvec = uc_pos @ qvec

        denominator = 1 + 2*np.pi * np.dot(qdir, np.dot(alpha2D, qdir)) * qnorm

        # calculate g for each phonon branch
        for lam in range(nmodes):
            ph_om = mesh_frequencies[iq][lam]
            ph_ev = mesh_eigenvectors[iq][:,lam].reshape((-1,3)) # atom, direction

            prefactor = (echrg/(4*np.pi*diel))**2 * hbar_amu_THz/(np.maximum(2*ph_om, 1e-5))

            g = np.zeros((nbands, nbands), complex) 
            for s in range(nuc): # atoms in unit cell
                m_s = np.sqrt(uc_masses[s])
                phase_s = np.exp(-2j*np.pi*pos_qvec[s])/m_s

                Zstar_lam = np.dot(Zstar[s], ph_ev[s])
                
                g_par = 1j * np.dot(qdir, Zstar_lam) * phase_s / denominator
                g_perp = -1 * np.dot(zdir, Zstar_lam) * phase_s / denominator
                
                M_kq_loc, S_kq_loc = calc_M_S_loc(s, uc_pos, pos_qvec, s0_k_uc, s0_kq_uc, uc2idx, norbitals)
                M_kq = 0.5 * ( U_kq.conj().T @ M_kq_loc @ U_k )
                S_kq = 0.5 * ( U_kq.conj().T @ S_kq_loc @ U_k )

                g += g_par * M_kq[band0:band1, band0:band1] + g_perp * S_kq[band0:band1, band0:band1] 

            if np.isrealobj(mesh_g2):
                mesh_g2[iq, lam, :, :] += prefactor * np.abs(g)**2
            else:
                mesh_g2[iq, lam, :, :] += np.sqrt(prefactor) * g
                                   
    return eps_k[band0:band1], mesh_epskq, mesh_epskmq, mesh_g2
