import numpy as np

try:
    from .extensions import alphaFnk, inv_tau_nk, inv_tau_nk_lam
except ModuleNotFoundError:
    print('Warning: C-extensions not available')

    def alphaFnk(n, eps, frequency_points, mesh_g2, mesh_epskq, mesh_frequencies, sigma):
        nbands = mesh_g2.shape[2]
        nmodes = mesh_g2.shape[1]
        nqpoints = mesh_g2.shape[0]
        
        alphaF = np.zeros((frequency_points.shape[0],))
        for iq in range(nqpoints):
            for lam in range(nmodes):
                for m in range(nbands):
                    
                    weight_el = np.exp( -(eps - mesh_epskq[iq,m])**2/(2*sigma**2) )/(np.sqrt(2*np.pi)*sigma)
                    for nf, f in enumerate(frequency_points):
                        weight_ph = np.exp( -(f - mesh_frequencies[iq,lam])**2/(2*sigma**2) )/(np.sqrt(2*np.pi)*sigma)

                        alphaF[nf] +=   weight_ph * mesh_g2[iq,lam,m,n] * weight_el
        return alphaF/nqpoints

    def bose(x):
        return 1/(np.exp(x) - 1)

    def fermi(x):
        return 1/(np.exp(x) + 1)

    def inv_tau_nk(n, eps, mu, kBT, mesh_g2, mesh_epskq, mesh_frequencies, sigma):
        nbands = mesh_g2.shape[2]
        nmodes = mesh_g2.shape[1]
        nqpoints = mesh_g2.shape[0]

        tau_p = 0.0
        tau_m = 0.0        
        for iq in range(nqpoints):
            for lam in range(nmodes):
                for m in range(nbands):

                    weight_p = np.exp( -(eps + mesh_frequencies[iq,lam] - mesh_epskq[iq,m])**2/(2*sigma**2) )/(np.sqrt(2*np.pi)*sigma)
                    weight_m = np.exp( -(eps - mesh_frequencies[iq,lam] - mesh_epskq[iq,m])**2/(2*sigma**2) )/(np.sqrt(2*np.pi)*sigma)                

                    tau_p += mesh_g2[iq,lam,m,n] * (bose(mesh_frequencies[iq,lam]/kBT) + fermi((mesh_epskq[iq,m]-mu)/kBT)) * weight_p
                    tau_m += mesh_g2[iq,lam,m,n] * (bose(mesh_frequencies[iq,lam]/kBT) + 1. - fermi((mesh_epskq[iq,m]-mu)/kBT)) * weight_m                
        return 2*np.pi*(tau_p+tau_m)/nqpoints, 2*np.pi*(tau_p)/nqpoints, 2*np.pi*(tau_m)/nqpoints

    def inv_tau_nk_lam(n, eps, mu, kBT, mesh_g2, mesh_epskq, mesh_frequencies, sigma):
        nbands = mesh_g2.shape[2]
        nmodes = mesh_g2.shape[1]
        nqpoints = mesh_g2.shape[0]

        tau_p = np.zeros((nmodes,))
        tau_m = np.zeros((nmodes,))
        for iq in range(nqpoints):
            for lam in range(nmodes):
                for m in range(nbands):

                    weight_p = np.exp( -(eps + mesh_frequencies[iq,lam] - mesh_epskq[iq,m])**2/(2*sigma**2) )/(np.sqrt(2*np.pi)*sigma)
                    weight_m = np.exp( -(eps - mesh_frequencies[iq,lam] - mesh_epskq[iq,m])**2/(2*sigma**2) )/(np.sqrt(2*np.pi)*sigma)                

                    tau_p[lam] += mesh_g2[iq,lam,m,n] * (bose(mesh_frequencies[iq,lam]/kBT) + fermi((mesh_epskq[iq,m]-mu)/kBT)) * weight_p
                    tau_m[lam] += mesh_g2[iq,lam,m,n] * (bose(mesh_frequencies[iq,lam]/kBT) + 1. - fermi((mesh_epskq[iq,m]-mu)/kBT)) * weight_m                
        return 2*np.pi*(tau_p+tau_m)/nqpoints, 2*np.pi*(tau_p)/nqpoints, 2*np.pi*(tau_m)/nqpoints
    
def phononDOS(frequency_points, mesh_frequencies, sigma):
    nmodes = mesh_frequencies.shape[1]
    nqpoints = mesh_frequencies.shape[0]

    DOS = np.zeros((frequency_points.shape[0],))
    for iq in range(nqpoints):
        for lam in range(nmodes):
            for nf, f in enumerate(frequency_points):
                weight_ph = np.exp( -(f - mesh_frequencies[iq,lam])**2/(2*sigma**2) )/(np.sqrt(2*np.pi)*sigma)
                DOS[nf] +=   weight_ph
    return DOS/nqpoints
