import os
import subprocess
import re
import numpy as np
import string

from .units import *
from .fileio import writegen 

ENERGY_PATTERN = re.compile(r"mermin_energy[^:]*:[^:]*:[^:]*:\s*(?P<value>\S+)")

# this is the standard order of orbitals as stated in the dftb+ manual
std_orbital_order = {'s': ['s'], 'p': ['py', 'pz', 'px'], 'd': ['dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2']}

def run_dftb(binary, reffile=None):
    """Run dftb calculations using the given binary.
    """
    with open('out.txt', "w") as outfile:
        with open('err.txt', "w") as errfile:
            subprocess.run([binary], stdout=outfile, stderr=errfile)
    if reffile:
        os.rename("autotest.tag", reffile)
    

def calculate_reference(binary, reffile, coords, specienames, species, origin,
                        latvecs, scc=False):
    """Calculates reference system
       adapted from calcderivs.py (DFTB+ developers group)
    """

    # write input geometry
    writegen("geo.gen", (specienames, species, coords, origin, latvecs))

    if scc:
        # read template file
        with open('dftb_in.template') as t:
            template = string.Template(t.read())

        # setup for scc calculation
        final_output = template.substitute(maxscciterations=100, 
                                           readinitialcharges='No', 
                                           writehs='No')
        # write file
        with open('dftb_in.hsd', 'w') as output:
            output.write(final_output)

        # run dftb+
        run_dftb(binary, reffile)

        # copy charges.bin for later
        os.system('cp charges.bin charges0.bin')
    
        # read template file
        with open('dftb_in.template') as t:
            template = string.Template(t.read())

        # setup for scc calculation
        final_output = template.substitute(maxscciterations=1, 
                                           readinitialcharges='Yes', 
                                           writehs='Yes')
        # write file
        with open('dftb_in.hsd', 'w') as output:
            output.write(final_output)
            
    # run dftb+
    run_dftb(binary, reffile)

def calculate_forces(binary, disp, coords, specienames, species, origin,
                     latvecs):
    """Calculates forces by finite differences
       unmodified from calcderivs.py (DFTB+ developers group)
    """

    energy = np.empty((2,), dtype=float)
    forces = np.empty((len(coords), 3), dtype=float)
    for iat in range(len(coords)):
        for ii in range(3):
            for jj in range(2):
                newcoords = np.array(coords)
                newcoords[iat][ii] += float(2 * jj - 1) * disp
                writegen("geo.gen", (specienames, species, newcoords, origin,
                                     latvecs))

                run_dftb(binary)

                fp = open("autotest.tag", "r")
                txt = fp.read()
                fp.close()
                match = ENERGY_PATTERN.search(txt)
                print("iat: %2d, ii: %2d, jj: %2d" % (iat, ii, jj))
                if match:
                    energy[jj] = float(match.group("value"))
                    print("energy:", energy[jj])
                else:
                    raise "No match found!"
            forces[iat][ii] = (energy[0] - energy[1]) / (2.0 * disp)
    return forces

def calculate_hamiltonian_derivs(binary, disp, atom_ids, coords, specienames, species, origin,
                     latvecs, scc=False):
    """Calculates derivatives of hamiltonian and overlap matrix by finite differences. SELECTED atoms with atom_ids are displaced by disp in ALL directions.
       Only the first k-point is currently used!
       based on calculate_forces from calcderivs.py (DFTB+ developers group)
       
       returns: ham_derivs (shape: (len(atom_ids), 3, NALLORB, NALLORB)) in eV/Å
                ovr_derivs (shape: (len(atom_ids), 3, NALLORB, NALLORB)) in 1/Å
    """

    with open('hamsqr1.dat') as f:
        _ = f.readline()
        second_line = f.readline().strip('#')

    temp_dat = second_line.split()

    real_elem = temp_dat[0] == 'T'
    NALLORB = int(temp_dat[1])
    NKPOINT = int(temp_dat[2])

    if real_elem:
        skip = 1
    else:
        skip = 2
        
    hamiltonian = np.empty((2, NALLORB, NALLORB), dtype=float)
    ham_derivs = np.empty((len(atom_ids), 3, NALLORB, NALLORB), dtype=float)
    
    overlap = np.empty((2, NALLORB, NALLORB), dtype=float)
    ovr_derivs = np.empty((len(atom_ids), 3, NALLORB, NALLORB), dtype=float)

    # check if charges from reference calculation exist
    if os.path.isfile('charges0.bin'):
        read_ref_charges = 'Yes'
    else:
        read_ref_charges = 'No'

    for ino, iat in enumerate(atom_ids): # iterate over atoms
    
        for ii in range(3): # iterate over cartesian coordinates
            
            for jj in range(2): # forward-backward displacement
                newcoords = np.array(coords)
                newcoords[iat][ii] += float(2 * jj - 1) * disp 
                writegen("geo.gen", (specienames, species, newcoords, origin,
                                     latvecs))

                if scc:
                    # read template file
                    with open('dftb_in.template') as t:
                        template = string.Template(t.read())

                    # copy charges0.bin from reference calculation
                    if read_ref_charges[0] == 'Y':
                        os.system('cp charges0.bin charges.bin')

                    # setup for scc calculation
                    final_output = template.substitute(maxscciterations=100, 
                                                       readinitialcharges=read_ref_charges, 
                                                       writehs='No')
                    # write file
                    with open('dftb_in.hsd', 'w') as output:
                        output.write(final_output)

                    # run dftb+
                    run_dftb(binary)

                    # read template file
                    with open('dftb_in.template') as t:
                        template = string.Template(t.read())

                    # setup for scc calculation
                    final_output = template.substitute(maxscciterations=1, 
                                                       readinitialcharges='Yes', 
                                                       writehs='Yes')
                    # write file
                    with open('dftb_in.hsd', 'w') as output:
                        output.write(final_output)

                # run dftb+
                run_dftb(binary)
                                      
                # read hamsqr file
                # TODO read all k-points
                # for n in range (1,NKPOINT+1):
                n = 1
                hamiltonian[jj] = np.loadtxt('hamsqr1.dat', skiprows= 2 + n*3 + (n-1)*NALLORB , max_rows=NALLORB)[:,::skip]
                overlap[jj] = np.loadtxt('oversqr.dat', skiprows= 2 + n*3 + (n-1)*NALLORB , max_rows=NALLORB)[:,::skip]
                    
            ham_derivs[ino][ii] = (hamiltonian[1] - hamiltonian[0]) / (2.0 * disp)
            ovr_derivs[ino][ii] = (overlap[1] - overlap[0]) / (2.0 * disp)
            
    return ham_derivs*HARTREE__EV, ovr_derivs
