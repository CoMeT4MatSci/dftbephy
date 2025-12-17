# distutils: language = c++

# cimport numpy as np
import numpy as np

cimport numpy as np

import cython

from libc.math cimport sqrt, exp

cdef extern from "<complex.h>" namespace "std" nogil:
    double complex exp(double complex z)
    double abs(double complex z)    

cdef double complex I = 1j
cdef double PI = np.pi
cdef double complex TWO_PI_I = 2*np.pi*I

DTYPE = np.float64

# fourier

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calculate_lattice_ft(double[:,:] ham0, double[:] kvec, long[:] uc2sc, long[:] sc2uc, long[:] sc2c, long[:] uc2idx, long[:] sc2idx, double[:,:,:,:] svecs, int[:,:] multi):
    """Calculate Fourier transform of supercell hamiltonian/overlap matrix ham0 at k-point kvec. 
       Information about the supercell and the primitive cell are obtained from phonopy's supercell and primitive.
       
       returns: complex matrix H(kvec)
    """    
    cdef Py_ssize_t nsc = sc2uc.shape[0]
    cdef Py_ssize_t nuc = uc2sc.shape[0]
    cdef Py_ssize_t norbitals = uc2idx[nuc]

    h0_uc = np.zeros((norbitals, norbitals), dtype=np.complex128)    
    cdef double complex[:,:] h0_view = h0_uc

    cdef Py_ssize_t ni, m, i, j, s, l, sp, lp, ns, nsp
    cdef int mul
    cdef double complex phase
    cdef double dot_vk

    for ni in range(nuc): # atoms in uc
        i = uc2sc[ni]
        s = sc2uc[i]
        l = sc2c[i]

        for j in range(nsc): # atoms in sc
            sp = sc2uc[j]
            lp = sc2c[j]
            
            mul = multi[lp,sc2uc[l]]
            phase = 0. + 0.j
            for m in range(mul):
                dot_vk  = svecs[lp,sc2uc[l],m,0] * kvec[0]
                dot_vk += svecs[lp,sc2uc[l],m,1] * kvec[1]
                dot_vk += svecs[lp,sc2uc[l],m,2] * kvec[2]
                phase += exp(TWO_PI_I*dot_vk)
            phase = phase/mul

            for ns in range(uc2idx[s+1]-uc2idx[s]):
                for nsp in range(uc2idx[sp+1]-uc2idx[sp]):
                    h0_view[uc2idx[s]+ns,uc2idx[sp]+nsp] = h0_view[uc2idx[s]+ns,uc2idx[sp]+nsp] + ham0[sc2idx[i]+ns,sc2idx[j]+nsp] * phase
    return h0_uc


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calculate_lattice_ft_derivative(double[:,:,:,:] ham_derivs, double[:] kvec, long[:] uc2sc, long[:] sc2uc, long[:] sc2c, long[:] uc2idx, long[:] sc2idx, double[:,:,:,:] svecs, int[:,:] multi):
    """Calculate Fourier transform of supercell hamiltonian/overlap matrix derivative ham_derivs at k-point kvec. 
       Information about the supercell and the primitive cell are obtained from phonopy's supercell and primitive.
       
       returns: complex matrix dH/dR(kvec)
    """
    cdef Py_ssize_t nsc = sc2uc.shape[0]
    cdef Py_ssize_t nuc = uc2sc.shape[0]
    cdef Py_ssize_t norbitals = uc2idx[nuc]    

    dHdR = np.zeros((3, norbitals, norbitals), dtype=np.complex128)    
    cdef double complex[:,:,:] dHdR_view = dHdR

    cdef Py_ssize_t ni, m, i, j, s, l, sp, lp, ns, nsp, alpha
    cdef int mul
    cdef double complex phase
    cdef double dot_vk

    for ni in range(nuc): # atoms in uc
        i = uc2sc[ni]
        s = sc2uc[i]
        l = sc2c[i]

        for j in range(nsc): # atoms in sc
            sp = sc2uc[j]
            lp = sc2c[j]
            
            mul = multi[lp,sc2uc[l]]
            phase = 0. + 0.j
            for m in range(mul):
                dot_vk  = svecs[lp,sc2uc[l],m,0] * kvec[0]
                dot_vk += svecs[lp,sc2uc[l],m,1] * kvec[1]
                dot_vk += svecs[lp,sc2uc[l],m,2] * kvec[2]
                phase += exp(TWO_PI_I*dot_vk)
            phase = phase/mul

            for alpha in range(3): # cartesian components
                for ns in range(uc2idx[s+1]-uc2idx[s]):
                    for nsp in range(uc2idx[sp+1]-uc2idx[sp]):
                        dHdR_view[alpha,uc2idx[s]+ns,uc2idx[sp]+nsp] = dHdR_view[alpha,uc2idx[s]+ns,uc2idx[sp]+nsp] + ham_derivs[s,alpha,sc2idx[i]+ns,sc2idx[j]+nsp] * phase
    return dHdR

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calculate_lattice_double_ft_derivative(double[:,:,:,:] ham_derivs, double[:] kvec, double[:] kvec2, long[:] uc2sc, long[:] sc2uc, long[:] sc2c, long[:] uc2idx, long[:] sc2idx, double[:,:,:,:] svecs, int[:,:] multi):
    """Calculate double Fourier transform of supercell hamiltonian/overlap matrix derivative ham_derivs at k-points kvec and kvec2.        
       
       returns: complex matrix dH/dR(s'', dir, kvec, kvec2)
    """    
    cdef Py_ssize_t nsc = sc2uc.shape[0]
    cdef Py_ssize_t nuc = uc2sc.shape[0]
    cdef Py_ssize_t norbitals = uc2idx[nuc]
    cdef Py_ssize_t ndisp = ham_derivs.shape[0]

    dHdR = np.zeros((ndisp, 3, norbitals, norbitals), dtype=np.complex128)    
    cdef double complex[:,:,:,:] dHdR_view = dHdR

    cdef Py_ssize_t ni, m, i, j, s, l, l0, sp, lp, ns, nsp, alpha
    cdef int mul, mul2
    cdef double complex phase, phase2
    cdef double dot_vk

    l0 = sc2c[uc2sc[0]] # cell idx of reference cell
    for i in range(nsc): # atoms in sc
        s = sc2uc[i]
        l = sc2c[i]

        mul = multi[l,sc2uc[l0]]
        phase = 0. + 0.j
        for m in range(mul):
            dot_vk  = svecs[l,sc2uc[l0],m,0] * kvec[0]
            dot_vk += svecs[l,sc2uc[l0],m,1] * kvec[1]
            dot_vk += svecs[l,sc2uc[l0],m,2] * kvec[2]
            phase += exp(-TWO_PI_I*dot_vk)
        phase = phase/mul

        for j in range(nsc): # atoms in sc
            sp = sc2uc[j]
            lp = sc2c[j]
            
            mul2 = multi[lp,sc2uc[l0]]
            phase2 = 0. + 0.j
            for m in range(mul):
                dot_vk  = svecs[lp,sc2uc[l0],m,0] * kvec2[0]
                dot_vk += svecs[lp,sc2uc[l0],m,1] * kvec2[1]
                dot_vk += svecs[lp,sc2uc[l0],m,2] * kvec2[2]
                phase2 += exp(TWO_PI_I*dot_vk)
            phase2 = phase2/mul2

            for spp in range(ndisp):
                for alpha in range(3): # cartesian components
                    for ns in range(uc2idx[s+1]-uc2idx[s]):
                        for nsp in range(uc2idx[sp+1]-uc2idx[sp]):
                            dHdR_view[spp,alpha,uc2idx[s]+ns,uc2idx[sp]+nsp] = dHdR_view[spp,alpha,uc2idx[s]+ns,uc2idx[sp]+nsp] + ham_derivs[spp,alpha,sc2idx[i]+ns,sc2idx[j]+nsp] * phase * phase2
    return dHdR

##################################
# epc
          
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calc_g_loc(complex[:,:] g_H_loc, complex[:,:] g_S_k_loc, complex[:,:] g_S_kq_loc, complex[:,:,:] dHdR_k, complex[:,:,:] dSdR_k, complex[:,:,:] dHdR_kq, complex[:,:,:] dSdR_kq, complex[:,:] ph_ev, long [:] uc2idx, double[:] uc_masses, double[:] pos_qvec):

    cdef Py_ssize_t nuc = ph_ev.shape[0]

    cdef double complex[:,:] g_H_view = g_H_loc
    cdef double complex[:,:] g_S_k_view = g_S_k_loc
    cdef double complex[:,:] g_S_kq_view = g_S_kq_loc    
    cdef Py_ssize_t s, sp, ns, nsp, alpha
    cdef double m_s, m_sp
    cdef double complex phase_s, phase_sp

    for s in range(nuc): # atoms in unit cell
        m_s = sqrt(uc_masses[s])
        phase_s = exp(TWO_PI_I*pos_qvec[s])/m_s

        for sp in range(nuc): # atoms in unit cell
            m_sp = sqrt(uc_masses[sp])
            phase_sp = exp(TWO_PI_I*pos_qvec[sp])/m_sp

            for ns in range(uc2idx[s+1]-uc2idx[s]):
                for nsp in range(uc2idx[sp+1]-uc2idx[sp]):
                    g_H_view[uc2idx[s]+ns,uc2idx[sp]+nsp] = 0. + 0.j
                    g_S_k_view[uc2idx[s]+ns,uc2idx[sp]+nsp] = 0. + 0.j
                    g_S_kq_view[uc2idx[s]+ns,uc2idx[sp]+nsp] = 0. + 0.j                    
                    for alpha in range(3): # cartesian coordinates
                        g_H_view[uc2idx[s]+ns,uc2idx[sp]+nsp] += phase_s*ph_ev[s,alpha] * dHdR_k[alpha,uc2idx[s]+ns,uc2idx[sp]+nsp] - phase_sp*ph_ev[sp,alpha] * dHdR_kq[alpha,uc2idx[s]+ns,uc2idx[sp]+nsp]
                        g_S_k_view[uc2idx[s]+ns,uc2idx[sp]+nsp]  += phase_s*ph_ev[s,alpha] * dSdR_k[alpha,uc2idx[s]+ns,uc2idx[sp]+nsp] 
                        g_S_kq_view[uc2idx[s]+ns,uc2idx[sp]+nsp] += phase_sp*ph_ev[sp,alpha] * dSdR_kq[alpha,uc2idx[s]+ns,uc2idx[sp]+nsp]
            
@cython.boundscheck(False)
@cython.wraparound(False)
def calc_g2_inner(double[:,:,:,:] mesh_g2, long iq, long lam, complex[:,:] gH, complex[:,:] gS_k, complex[:,:] gS_kq, double[:] eps_k, double[:] eps_kq, double prefactor, long band0, long band1):
    cdef double[:,:,:,:] mesh_g2_view = mesh_g2

    cdef Py_ssize_t nbands = mesh_g2.shape[2]
    cdef Py_ssize_t jband, kband, k, j
    cdef double g2
    for j in range(nbands):
        jband = j + band0
        for k in range(nbands):
            kband = k + band0
            g2 = prefactor * abs(gH[jband,kband] - (eps_k[kband]*gS_k[jband,kband] - eps_kq[jband]*gS_kq[jband,kband]))**2
            mesh_g2_view[iq, lam, j, k] = g2
            
# NOT USED                        
@cython.boundscheck(False)
@cython.wraparound(False)
def calc_g2(double[:,:,:,:] mesh_g2, long iq, long lam, complex[:,:] g_H_loc, complex[:,:] g_S_loc, double[:] eps_k, complex[:,:] U_k, complex[:,:] U_kq_H, double prefactor):

    cdef Py_ssize_t nbands = mesh_g2.shape[2]
    cdef Py_ssize_t nuc = g_H_loc.shape[0]
    cdef Py_ssize_t jband, kband
    cdef double[:,:,:,:] mesh_g2_view = mesh_g2
    cdef double complex[:,:] U_kq_H_view = U_kq_H
    cdef double complex[:,:] U_k_view = U_k
    cdef double complex[:,:] g_H_loc_view = g_H_loc
    cdef double complex[:,:] g_S_loc_view = g_S_loc

    cdef double g2
    cdef double complex gH, gS

    for jband in range(nbands):
        for kband in range(nbands):

            gH = 0. + 0.j
            gS = 0. + 0.j
            for s in range(nuc):
                for sp in range(nuc):
                    gH += U_kq_H_view[jband,s] * g_H_loc_view[s,sp] * U_k_view[sp,kband]
                    gS += U_kq_H_view[jband,s] * g_S_loc_view[s,sp] * U_k_view[sp,kband]

            g2 = prefactor * abs(gH - eps_k[kband]*gS)**2
            mesh_g2_view[iq, lam, jband, kband] = g2


##################################
# analysis

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def alphaFnk(long n, double eps, double[:] frequency_points, double[:,:,:,:] mesh_g2, double[:,:] mesh_epskq, double[:,:] mesh_frequencies, double sigma):
    cdef Py_ssize_t nbands = mesh_g2.shape[2]
    cdef Py_ssize_t nmodes = mesh_g2.shape[1]
    cdef Py_ssize_t nqpoints = mesh_g2.shape[0]
    cdef Py_ssize_t nfpoints = frequency_points.shape[0]
    
    alphaF = np.zeros((frequency_points.shape[0],), dtype=np.float64)
    cdef double[:] alphaF_view = alphaF
    cdef double[:,:] mesh_epskq_view = mesh_epskq
    cdef double[:,:] mesh_frequencies_view = mesh_frequencies

    cdef Py_ssize_t iq, lam, m, nf
    cdef double weight_ph, weight_el, f
    for iq in range(nqpoints):
        for lam in range(nmodes):
            for m in range(nbands):
                
                weight_el = exp( -(eps - mesh_epskq_view[iq,m])**2/(2*sigma**2) )/(sqrt(2*PI)*sigma)
                for nf in range(nfpoints):
                    f = frequency_points[nf]
                    weight_ph = exp( -(f - mesh_frequencies_view[iq,lam])**2/(2*sigma**2) )/(sqrt(2*PI)*sigma)

                    alphaF_view[nf] +=   weight_ph * mesh_g2[iq,lam,m,n] * weight_el / nqpoints
    return alphaF

@cython.cdivision(True)
cdef double bose(double x):
    return 1.0/(exp(x) - 1.0)

@cython.cdivision(True)
cdef double fermi(double x):
    return 1.0/(exp(x) + 1.0)

@cython.cdivision(True)
cdef double lorentzian(double x, double sigma):
    return (sigma/2.0)/(x*x + (sigma/2.0)*(sigma/2.0))/PI

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def inv_tau_nk_lor(long n, double eps, double mu, double kBT, double[:,:,:,:] mesh_g2, double[:,:] mesh_epskq, double[:,:] mesh_frequencies, double sigma):
    cdef Py_ssize_t nbands = mesh_g2.shape[2]
    cdef Py_ssize_t nmodes = mesh_g2.shape[1]
    cdef Py_ssize_t nqpoints = mesh_g2.shape[0]

    cdef double tau_p, tau_m, weight_p, weight_m
    cdef Py_ssize_t iq, lam, m

    tau_p = 0.0
    tau_m = 0.0
    for iq in range(nqpoints):
        for lam in range(nmodes):
            for m in range(nbands):

                weight_p = lorentzian( eps + mesh_frequencies[iq,lam] - mesh_epskq[iq,m], sigma )
                weight_m = lorentzian( eps - mesh_frequencies[iq,lam] - mesh_epskq[iq,m], sigma )

                tau_p += mesh_g2[iq,lam,m,n] * (bose(mesh_frequencies[iq,lam]/kBT) + fermi((mesh_epskq[iq,m]-mu)/kBT)) * weight_p
                tau_m += mesh_g2[iq,lam,m,n] * (bose(mesh_frequencies[iq,lam]/kBT) + 1. - fermi((mesh_epskq[iq,m]-mu)/kBT)) * weight_m
    return 2*PI*(tau_p+tau_m)/nqpoints, 2*PI*tau_p/nqpoints, 2*PI*tau_m/nqpoints

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def inv_tau_nk_lam_lor(long n, double eps, double mu, double kBT, double[:,:,:,:] mesh_g2, double[:,:] mesh_epskq, double[:,:] mesh_epskmq, double[:,:] mesh_frequencies, double sigma):
    cdef Py_ssize_t nbands = mesh_g2.shape[2]
    cdef Py_ssize_t nmodes = mesh_g2.shape[1]
    cdef Py_ssize_t nqpoints = mesh_g2.shape[0]

    cdef double weight_p, weight_m
    cdef Py_ssize_t iq, lam, m

    tau = np.zeros(nmodes, dtype=DTYPE)
    tau_p = np.zeros(nmodes, dtype=DTYPE)
    tau_m = np.zeros(nmodes, dtype=DTYPE)
    cdef double[:] tau_view = tau  
    cdef double[:] tau_p_view = tau_p
    cdef double[:] tau_m_view = tau_m
    for lam in range(nmodes):
        for iq in range(nqpoints):
            for m in range(nbands):

                weight_p = lorentzian( eps + mesh_frequencies[iq,lam] - mesh_epskq[iq,m], sigma )
                weight_m = lorentzian( eps - mesh_frequencies[iq,lam] - mesh_epskmq[iq,m], sigma )

                tau_p_view[lam] += mesh_g2[iq,lam,m,n] * (bose(mesh_frequencies[iq,lam]/kBT) + fermi((mesh_epskq[iq,m]-mu)/kBT)) * weight_p
                tau_m_view[lam] += mesh_g2[iq,lam,m,n] * (bose(mesh_frequencies[iq,lam]/kBT) + 1. - fermi((mesh_epskmq[iq,m]-mu)/kBT)) * weight_m
        tau_p_view[lam] = tau_p_view[lam]*2*PI/nqpoints
        tau_m_view[lam] = tau_m_view[lam]*2*PI/nqpoints
        tau_view[lam] = tau_p_view[lam] + tau_m_view[lam]
    return tau, tau_p, tau_m

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def inv_tau_nk(long n, double eps, double mu, double kBT, double[:,:,:,:] mesh_g2, double[:,:] mesh_epskq, double[:,:] mesh_frequencies, double sigma):
    cdef Py_ssize_t nbands = mesh_g2.shape[2]
    cdef Py_ssize_t nmodes = mesh_g2.shape[1]
    cdef Py_ssize_t nqpoints = mesh_g2.shape[0]

    cdef double tau_p, tau_m, weight_p, weight_m
    cdef Py_ssize_t iq, lam, m
    
    tau_p = 0.0
    tau_m = 0.0
    for iq in range(nqpoints):
        for lam in range(nmodes):
            for m in range(nbands):

                weight_p = exp( -(eps + mesh_frequencies[iq,lam] - mesh_epskq[iq,m])**2/(2*sigma**2) )/(sqrt(2*PI)*sigma)
                weight_m = exp( -(eps - mesh_frequencies[iq,lam] - mesh_epskq[iq,m])**2/(2*sigma**2) )/(sqrt(2*PI)*sigma)

                tau_p += mesh_g2[iq,lam,m,n] * (bose(mesh_frequencies[iq,lam]/kBT) + fermi((mesh_epskq[iq,m]-mu)/kBT)) * weight_p
                tau_m += mesh_g2[iq,lam,m,n] * (bose(mesh_frequencies[iq,lam]/kBT) + 1. - fermi((mesh_epskq[iq,m]-mu)/kBT)) * weight_m                
    return 2*PI*(tau_p+tau_m)/nqpoints, 2*PI*tau_p/nqpoints, 2*PI*tau_m/nqpoints

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def inv_tau_nk_lam(long n, double eps, double mu, double kBT, double[:,:,:,:] mesh_g2, double[:,:] mesh_epskq, double[:,:] mesh_frequencies, double sigma):
    cdef Py_ssize_t nbands = mesh_g2.shape[2]
    cdef Py_ssize_t nmodes = mesh_g2.shape[1]
    cdef Py_ssize_t nqpoints = mesh_g2.shape[0]

    cdef double weight_p, weight_m
    cdef Py_ssize_t iq, lam, m
    
    tau = np.zeros(nmodes, dtype=DTYPE)
    tau_p = np.zeros(nmodes, dtype=DTYPE)
    tau_m = np.zeros(nmodes, dtype=DTYPE)
    cdef double[:] tau_view = tau  
    cdef double[:] tau_p_view = tau_p
    cdef double[:] tau_m_view = tau_m
    for lam in range(nmodes):
        for iq in range(nqpoints):
            for m in range(nbands):

                weight_p = exp( -(eps + mesh_frequencies[iq,lam] - mesh_epskq[iq,m])**2/(2*sigma**2) )/(sqrt(2*PI)*sigma)
                weight_m = exp( -(eps - mesh_frequencies[iq,lam] - mesh_epskq[iq,m])**2/(2*sigma**2) )/(sqrt(2*PI)*sigma)

                tau_p_view[lam] += mesh_g2[iq,lam,m,n] * (bose(mesh_frequencies[iq,lam]/kBT) + fermi((mesh_epskq[iq,m]-mu)/kBT)) * weight_p
                tau_m_view[lam] += mesh_g2[iq,lam,m,n] * (bose(mesh_frequencies[iq,lam]/kBT) + 1. - fermi((mesh_epskq[iq,m]-mu)/kBT)) * weight_m                
        tau_p_view[lam] = tau_p_view[lam]*2*PI/nqpoints
        tau_m_view[lam] = tau_m_view[lam]*2*PI/nqpoints
        tau_view[lam] = tau_p_view[lam] + tau_m_view[lam] 
    return tau, tau_p, tau_m

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def epa_inner(double[:] energy_points, double[:,:,:,:] g2_kq, double[:] energies_k, double[:,:] energies_kq, double sigma_0, double[:,:,:] g2Epa, double[:,:] epa_norm):

    cdef Py_ssize_t nbands = g2_kq.shape[2]
    cdef Py_ssize_t nmodes = g2_kq.shape[1]
    cdef Py_ssize_t nqpoints = g2_kq.shape[0]
    cdef Py_ssize_t nenergy_points = energy_points.shape[0]

    cdef double[:,:,:] g2Epa_view = g2Epa
    cdef double[:,:] epa_norm_view = epa_norm    

    cdef double ei, ej, weights_i, weights_j
    cdef Py_ssize_t iq, i, j, lam, n, m
    for iq in range(nqpoints):    
        for i in range(nenergy_points):
            for j in range(nenergy_points):
                ei = energy_points[i]
                ej = energy_points[j]

                for n in range(nbands):
                    weights_i = exp( -(ei - energies_k[n])**2/(2*sigma_0**2) )

                    for m in range(nbands):
                        weights_j = exp( -(ej - energies_kq[iq,m])**2/(2*sigma_0**2) )

                        epa_norm_view[i,j] += weights_i * weights_j

                        for lam in range(nmodes):
                            g2Epa_view[lam,i,j] = g2Epa_view[lam,i,j] + weights_i * g2_kq[iq, lam, m, n] * weights_j
