import numpy as np
from scipy import linalg

from phonopy.structure.cells import get_smallest_vectors

from .units import *
from .atomic import *
from .fileio import read_hamsqr1, read_oversqr, read_tagfile, readgen
from .dftb import std_orbital_order, calculate_reference, calculate_hamiltonian_derivs
from .epc import calculate_g2, calculate_scc_g2
from .fourier import calculate_lattice_ft_kderivative

try:
    from .extensions import calculate_lattice_ft, calculate_lattice_ft_derivative
except ModuleNotFoundError:
    print('Warning: C-extensions not available')
    from .fourier import calculate_lattice_ft, calculate_lattice_ft_derivative

class DftbSuperCellCalc:
    def __init__(self, angular_momenta, cmd='dftb+'):

        self.supercell = None
        self.primitive = None
        self.uc2sc = None
        self.sc2c = None
        self.sc2uc = None

        self.Nsc = None
        
        self.sc2idx = None
        self.uc2idx = None

        self.svecs = None
        self.multi = None
        
        self.angular_momenta = angular_momenta
        self.cmd = cmd
        
        self.H0 = None
        self.S0 = None
        self.H_derivs = None
        self.S_derivs = None
        
    def load_phonopy(self, ph):
        """Initialize calculator from existing phonopy results.

            ph: Phonopy object
        """
        self.primitive = ph.get_primitive()
        self.supercell = ph.get_supercell()
        Nsc = self.supercell.get_supercell_matrix().diagonal().prod() # number of primitive cells in super cell
        self.Nsc = Nsc
        Nuc = len(self.primitive) # number of atoms in primitive cell
        
        # map of atom in super cell to the cell index
        self.sc2c = np.tile(np.arange(0,Nsc), Nuc)

        uc2uc = self.supercell.get_unitcell_to_unitcell_map()
        self.sc2uc = np.array([uc2uc[at] for at in self.supercell.get_supercell_to_unitcell_map()])
        
        # indices of primitive cell atoms in the super cell:
        # ideally, we would use
#        self.uc2sc = self.supercell.get_unitcell_to_supercell_map()
        # but for now we work around this 
        self.uc2sc = np.zeros((Nuc,), int)
        self.uc2sc[0] = (Nsc-1)//2
        for i in np.arange(1,Nuc):
            self.uc2sc[i] = self.uc2sc[i-1] + Nsc

        # obtain a map of atom index to Hamiltonian index
        orbitals = [sum([std_orbital_order[am] for am in self.angular_momenta[cs]], []) for cs in self.supercell.get_chemical_symbols()]
        norbitals = [len(o) for o in orbitals]
        self.sc2idx = np.insert(np.cumsum(norbitals), 0, 0)    

        uc_orbitals = [orbitals[i] for i in self.uc2sc]
        norbitals = [len(o) for o in uc_orbitals]
        self.uc2idx = np.insert(np.cumsum(norbitals), 0, 0)    

        # get vectors between atoms in primitive cell coordinates
        self.svecs, self.multi = get_smallest_vectors(self.supercell.get_cell(), self.supercell.scaled_positions, self.supercell.scaled_positions[self.uc2sc] )
        trans_mat_float = np.dot(self.supercell.get_cell(), np.linalg.inv(self.primitive.get_cell()))
        trans_mat = np.rint(trans_mat_float).astype(int)
        assert (np.abs(trans_mat_float - trans_mat) < 1e-8).all()
        self.svecs = np.array(np.dot(self.svecs, trans_mat), dtype="double", order="C")

        # reset matrices
        self.H0 = None
        self.S0 = None
        self.H_derivs = None
        self.S_derivs = None
            
    def calculate_reference(self, scc=False):
        """Obtain Hamiltonian and overlap matrix for the supercell from DFTB. This system serves as a unperturbed reference.

            scc: Perform SCC calculation to obtain charges before getting the matrices (boolean).
        """
        if (self.supercell is None):
            raise RuntimeError("Phonopy supercell not available. Use load_phonopy() first.")
            
        coords = self.supercell.get_positions() * BOHR__AA
        specienames = list(dict.fromkeys(self.supercell.get_chemical_symbols()))
        species = [specienames.index(sym) for sym in self.supercell.get_chemical_symbols()]
        origin = np.array([0., 0., 0.])
        latvecs = self.supercell.get_cell()*BOHR__AA

        calculate_reference(self.cmd, None, coords, specienames, species, origin, latvecs, scc)
        self.H0 = read_hamsqr1()
        self.S0 = read_oversqr()
        
    def calculate_derivatives(self, disp=0.001, scc=False):
        """Calculate derivatives of Hamiltonian and overlap matrix via finite-differences. All atoms in the primitive cell are displaced
            in all cartesian directions by disp (in Angstrom).

            disp: amount of displacement in Angstrom (scalar)
            scc: Perform SCC calculation to obtain charges for each displacement before getting the matrices (boolean).
        """
        if (self.supercell is None):
            raise RuntimeError("Phonopy supercell not available. Use load_phonopy() first.")
            
        coords = self.supercell.get_positions() * BOHR__AA
        specienames = list(dict.fromkeys(self.supercell.get_chemical_symbols()))
        species = [specienames.index(sym) for sym in self.supercell.get_chemical_symbols()]
        origin = np.array([0., 0., 0.])
        latvecs = self.supercell.get_cell()*BOHR__AA
        
        self.H_derivs, self.S_derivs = calculate_hamiltonian_derivs(
                self.cmd, disp, self.uc2sc, coords, specienames, species, origin, latvecs, scc)
    
    def calculate_band_structure(self, kpoints):
        """Calculate the electronic band-structure at the given k-points.

            kpoints: list containing arrays of k-points (nkpoints, 3), one for each path

            returns: bands = list of arrays of energies (nkpoints, nbands), one for each path
        """
        if (self.supercell is None):
            raise RuntimeError("Phonopy supercell not available. Use load_phonopy() first.")
        
        if (self.H0 is None) or (self.S0 is None):
            raise RuntimeError("Reference Hamiltonian and Overlap Matrices not available. Run calculate_reference() first.")
        
        npaths = len(kpoints)
        orbitals = [sum([std_orbital_order[am] for am in self.angular_momenta[cs]], []) for cs in self.supercell.get_chemical_symbols()]
        uc_orbitals = [orbitals[i] for i in self.uc2sc]
        nbands = sum([len(o) for o in uc_orbitals])
        
        bands = []        
        for i in range(npaths):
            nkpoints = kpoints[i].shape[0]
            energies = np.zeros((nkpoints, nbands), float)
            
            for ik, kvec in enumerate(kpoints[i]):
                h0_uc = calculate_lattice_ft(self.H0, kvec, self.uc2sc, self.sc2uc, self.sc2c, self.uc2idx, self.sc2idx, self.svecs, self.multi)
                s0_uc = calculate_lattice_ft(self.S0, kvec, self.uc2sc, self.sc2uc, self.sc2c, self.uc2idx, self.sc2idx, self.svecs, self.multi)

                energies[ik,:] = linalg.eigvalsh(h0_uc, b=s0_uc)
                
            bands.append(energies)
        
        return bands
    
    def calculate_velocity(self, kvec0, band_sel=None):
        """Calculate the electron velocities at kvec0. The electronic bands can be selected with band_sel=[band0,band1+1].

            returns: eps_k = electronic energies at kvec0 (nbands)
                     velocities_k = electron velocities at kvec0 (3,nbands)
        """
        if (self.supercell is None):
            raise RuntimeError("Phonopy supercell not available. Use load_phonopy() first.")
        
        if (self.H0 is None) or (self.S0 is None):
            raise RuntimeError("Reference Hamiltonian and Overlap Matrices not available. Run calculate_reference() first.")

        h0_uc = calculate_lattice_ft(self.H0, kvec0, self.uc2sc, self.sc2uc, self.sc2c, self.uc2idx, self.sc2idx, self.svecs, self.multi)
        s0_uc = calculate_lattice_ft(self.S0, kvec0, self.uc2sc, self.sc2uc, self.sc2c, self.uc2idx, self.sc2idx, self.svecs, self.multi)

        dh0_dk = calculate_lattice_ft_kderivative(self.H0, kvec0, self.uc2sc, self.sc2uc, self.sc2c, self.uc2idx, self.sc2idx, self.svecs, self.multi)
        ds0_dk = calculate_lattice_ft_kderivative(self.S0, kvec0, self.uc2sc, self.sc2uc, self.sc2c, self.uc2idx, self.sc2idx, self.svecs, self.multi)

        eps_k, U_k = linalg.eigh(h0_uc, b=s0_uc) # diagonalize Hamiltonian
        nbands = eps_k.shape[0] # no of electronic bands

        if (type(band_sel) is list) and (len(band_sel) == 2):
            band0 = max(band_sel[0], 0)
            band1 = min(band_sel[1], nbands)
            nbands = band1 - band0
        else:
            band0 = 0
            band1 = nbands
        
        velocities_k = np.zeros((3, nbands), float)
        for n in range(nbands):
            u = U_k[:,band0+n]
            for alpha in range(3):
                velocities_k[alpha, n] = np.real( np.einsum('i,ij,j', u.conj(), dh0_dk[alpha], u) - eps_k[band0+n]*np.einsum('i,ij,j', u.conj(), ds0_dk[alpha], u) )

        # the derivative is with respect to fractional coordinates
        # thus, we have to convert it back to cartesian using
        #     k_cart = 2*pi*inv(cell)*k_frac
        velocities_k = (self.primitive.cell.T*BOHR__AA) @ velocities_k / (2*np.pi)
        return eps_k[band0:band1], velocities_k

    # WORK IN PROGRESS!!
    def calculate_transition_dipole_moment(self, kvec0):
        """Calculate the transition dipole moment at kvec0.

            returns: eps_k = electronic energies at kvec0 (nbands)
                     velocities_k = electron velocities at kvec0 (3,nbands)
        """
        if (self.supercell is None):
            raise RuntimeError("Phonopy supercell not available. Use load_phonopy() first.")
        
        if (self.H0 is None) or (self.S0 is None):
            raise RuntimeError("Reference Hamiltonian and Overlap Matrices not available. Run calculate_reference() first.")

        if (self.S_derivs is None):
            raise RuntimeError("Derivatives of Hamiltonian and Overlap matrices not available. Run calculate_derivatives() first.")
            

        h0_uc = calculate_lattice_ft(self.H0, kvec0, self.uc2sc, self.sc2uc, self.sc2c, self.uc2idx, self.sc2idx, self.svecs, self.multi)
        s0_uc = calculate_lattice_ft(self.S0, kvec0, self.uc2sc, self.sc2uc, self.sc2c, self.uc2idx, self.sc2idx, self.svecs, self.multi)
        dSdR_k = calculate_lattice_ft_derivative(self.S_derivs, kvec0, self.uc2sc, self.sc2uc, self.sc2c, self.uc2idx, self.sc2idx, self.svecs, self.multi)
        
        eps_k, U_k = linalg.eigh(h0_uc, b=s0_uc) # diagonalize Hamiltonian
        for iband in range(U_k.shape[1]): # fix gauge by choosing first element of each vector to be real
            U_k[:,iband] = U_k[:,iband] * np.exp( -1j * np.angle(U_k[0,iband]) )

        nbands = eps_k.shape[0] # no of electronic bands
        tdm_k = np.zeros((3, nbands, nbands), complex)
        for alpha in range(3):
            tdm_k[alpha, :, :] = -1j *U_k.conj().T @ dSdR_k[alpha, :, :] @ U_k
        
        for n in range(nbands):
            for m in range(nbands):                
                if n == m:
                    tdm_k[:, n, m] = 0.+ 0j
                else:                    
                    tdm_k[:, n, m] = -1j*tdm_k[:, n, m]/(eps_k[n]-eps_k[m])
                    
        # hbar [eV s] * hbar [kg m^2 / s] / ( m_e [kg] * Å ) / Å
        prefactor =  ( 6.582119569e-16 * 1.054571817e-34 / (9.1093837e-31 * 1e-10) ) * 1e10

        return eps_k, U_k, prefactor*tdm_k

    def calculate_g2(self, kvec0, qpoints, ph_frequencies, ph_eigenvectors, band_sel=None, scc=False):
        """Calculate the squared electron-phonon couplings at kvec0 and kvec0+qpoints using the specified phonon frequencies (ph_frequencies in THz) and phonon eigenmodes (ph_eigenvectors). The electronic bands can be selected with band_sel=[band0,band1+1].

            returns: eps_k = electronic energies at kvec0 (nqpoints)
                     mesh_epskq =  electronic energies at kvec0+qpoints (nqpoints, nbands)
                     mesh_g2 = squared electron-phonon coupling matrix (nqpoints, nmodes, nbands, nbands)
        """
        if (self.H0 is None) or (self.S0 is None):
            raise RuntimeError("Reference Hamiltonian and Overlap matrices not available. Run calculate_reference() first.")
        if (self.H_derivs is None) or (self.S_derivs is None):
            raise RuntimeError("Derivatives of Hamiltonian and Overlap matrices not available. Run calculate_derivatives() first.")

        if scc==True:
            eps_k, mesh_epskq, mesh_epskmq, mesh_g2 = calculate_scc_g2(kvec0, band_sel, qpoints, ph_frequencies, ph_eigenvectors, self.uc2sc, self.sc2uc, self.sc2c, self.supercell, self.primitive, self.angular_momenta, self.H0, self.S0, self.H_derivs, self.S_derivs)

        else:
            eps_k, mesh_epskq, mesh_epskmq, mesh_g2 = calculate_g2(kvec0, band_sel, qpoints, ph_frequencies, ph_eigenvectors, self.uc2sc, self.sc2uc, self.sc2c, self.supercell, self.primitive, self.angular_momenta, self.H0, self.S0, self.H_derivs, self.S_derivs)

        return eps_k, mesh_epskq, mesh_epskmq, mesh_g2

    def get_num_bands(self):
        """Return number of electronic bands.
        """
        if (self.supercell is None):
            raise RuntimeError("Phonopy supercell not available. Use load_phonopy() first.")
        
        if (self.H0 is None) or (self.S0 is None):
            raise RuntimeError("Reference Hamiltonian and Overlap Matrices not available. Run calculate_reference() first.")
        return int(self.H0.shape[0]/self.Nsc)
        

# WORK IN PROGRESS!!
class DftbMoleculeCalc:
    def __init__(self, angular_momenta):
                
        self.angular_momenta = angular_momenta

        # geometry of the Molecule
        self.specienames = None
        self.species = None
        self.coords = None
        self.origin = None
        self.latvecs  = None

        self.masses = None

        # Hamiltonian, overlap matrices and their derivatives
        self.H0 = None
        self.S0 = None
        self.H_derivs = None
        self.S_derivs = None

        # vibrational modes
        self.modes_dict = None
        
    def load_modes(self, fname = 'vibrations.tag'):
        self.modes_dict = read_tagfile(fname)

    def calculate_reference(self):
        self.specienames, self.species, self.coords, self.origin, self.latvecs = readgen('geo.gen.template')

        self.masses = [atom_data[symbol_map[self.specienames[s]]][3] for s in self.species]

        calculate_reference('dftb+', 'reference.tag', self.coords, self.specienames, self.species, self.origin, self.latvecs)
        self.H0 = read_hamsqr1()
        self.S0 = read_oversqr()
        
    def calculate_derivatives(self, disp=0.001):
        if (self.specienames is None):
            raise RuntimeError("Geometry information not available. Use calculate_reference() first.")
                
        self.H_derivs, self.S_derivs = calculate_hamiltonian_derivs(
                'dftb+', disp, np.arange(len(self.species)), self.coords, self.specienames, self.species, self.origin, self.latvecs)

    def calculate_g2(self):
        if (self.modes_dict is None):
            raise RuntimeError("Vibrational modes not available. Run load_modes() first.")
        if (self.H0 is None) or (self.S0 is None):
            raise RuntimeError("Reference Hamiltonian and Overlap matrices not available. Run calculate_reference() first.")
        if (self.H_derivs is None) or (self.S_derivs is None):
            raise RuntimeError("Derivatives of Hamiltonian and Overlap matrices not available. Run calculate_derivatives() first.")

        eps_n, mesh_g2 = calculate_g2_eigen(self.modes_dict, self.H0, self.S0, self.H_derivs, self.S_derivs, self.species, self.specienames, self.masses, self.angular_momenta)

        return eps_n, mesh_g2
