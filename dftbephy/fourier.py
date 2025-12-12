import numpy as np

from phonopy.structure.cells import get_smallest_vectors

from .dftb import std_orbital_order

def calculate_lattice_ft(ham0, kvec, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi):
    """Calculate Fourier transform of supercell hamiltonian/overlap matrix ham0 at k-point kvec. 
       
       returns: complex matrix H(kvec)
    """
    norbitals = uc2idx[len(uc2sc)]
    h0_uc = np.zeros((norbitals, norbitals), dtype=complex)
    for i in uc2sc: # atoms in uc
        s = sc2uc[i]
        l = sc2c[i]

        for j in np.arange(len(sc2uc)): # atoms in sc
            sp = sc2uc[j]
            lp = sc2c[j]
            
            mul = multi[lp,sc2uc[l]]
            vecs = svecs[lp,sc2uc[l],:mul]
            phase = np.exp(2j*np.pi*np.dot( vecs, kvec)).sum()/mul

            h0_uc[uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]] += ham0[sc2idx[i]:sc2idx[i+1],sc2idx[j]:sc2idx[j+1]] * phase
    return h0_uc

def calculate_lattice_ft_derivative(ham_derivs, kvec, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi):
    """Calculate Fourier transform of supercell hamiltonian/overlap matrix derivative ham_derivs at k-point kvec.        
       
       returns: complex matrix dH/dR(kvec)
    """    
    norbitals = uc2idx[len(uc2sc)]
    dHdR = np.zeros((3, norbitals, norbitals), dtype=complex)
    for i in uc2sc:  # atoms in uc
        s = sc2uc[i]
        l = sc2c[i]

        for j in np.arange(len(sc2uc)): # atoms in sc
            sp = sc2uc[j]
            lp = sc2c[j]
            
            mul = multi[lp,sc2uc[l]]
            vecs = svecs[lp,sc2uc[l],:mul]
            phase = np.exp(2j*np.pi*np.dot( vecs, kvec)).sum()/mul
            
            dHdR[:, uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]] += ham_derivs[s, :, sc2idx[i]:sc2idx[i+1],sc2idx[j]:sc2idx[j+1]]* phase
    return dHdR

def calculate_lattice_double_ft_derivative(ham_derivs, kvec, kvec2, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi):
    """Calculate double Fourier transform of supercell hamiltonian/overlap matrix derivative ham_derivs at k-points kvec and kvec2.        
       
       returns: complex matrix dH/dR(s'', dir, kvec, kvec2)
    """    
    norbitals = uc2idx[len(uc2sc)]
    ndisp = ham_derivs.shape[0] # number of displaced atoms in reference cell
    l0 = sc2c[uc2sc[0]]  # cell idx of reference cell

    dHdR = np.zeros((ndisp, 3, norbitals, norbitals), dtype=complex)
    for i in np.arange(len(sc2uc)): # atoms in sc
        s = sc2uc[i] # atom idx in uc
        l = sc2c[i]  # cell idx
        
        mul = multi[l,sc2uc[l0]]
        vecs = svecs[l,sc2uc[l0],:mul]
        phase = np.exp(-2j*np.pi*np.dot( vecs, kvec)).sum()/mul

        for j in np.arange(len(sc2uc)): # atoms in sc
            sp = sc2uc[j] # atom idx in uc
            lp = sc2c[j]  # cell idx
            
            mul2 = multi[lp,sc2uc[l0]]
            vecs2 = svecs[lp,sc2uc[l0],:mul2]
            phase2 = np.exp(2j*np.pi*np.dot( vecs2, kvec2)).sum()/mul2

            dHdR[:, :, uc2idx[s]:uc2idx[s+1],uc2idx[sp]:uc2idx[sp+1]] += ham_derivs[:, :, sc2idx[i]:sc2idx[i+1], sc2idx[j]:sc2idx[j+1]] * phase * phase2
    return dHdR


def calculate_lattice_ft_kderivative(ham0, kvec, uc2sc, sc2uc, sc2c, uc2idx, sc2idx, svecs, multi):
    """Calculate Fourier transform of supercell hamiltonian/overlap matrix ham0 at k-point kvec. 
       
       returns: complex matrix dH(kvec)/dkvec
    """
    norbitals = uc2idx[len(uc2sc)]
    dh0_dk = np.zeros((3, norbitals, norbitals), dtype=complex)
    for i in uc2sc: # atoms in uc
        s = sc2uc[i]
        l = sc2c[i]

        for j in np.arange(len(sc2uc)): # atoms in sc
            sp = sc2uc[j]
            lp = sc2c[j]
            
            mul = multi[lp,sc2uc[l]]
            vecs = svecs[lp,sc2uc[l],:mul]
            phases = np.exp(2j*np.pi*np.dot( vecs, kvec))

            for m in range(mul):
                for alpha in range(3):
                    dh0_dk[alpha, uc2idx[s]:uc2idx[s+1], uc2idx[sp]:uc2idx[sp+1]] += 2j*np.pi * vecs[m,alpha] * ham0[sc2idx[i]:sc2idx[i+1],sc2idx[j]:sc2idx[j+1]] * phases[m]/mul
    return dh0_dk
