import os
import numpy as np

from .units import *

def readgen(fname):
    """Reads in the content of a gen file. Units of length are Angström.
       modified from calcderivs.py (DFTB+ developers group)
    
       fname:   filename
       
       returns: specienames, species, coords, origin, latvecs
    """

    fp = open(fname, "r")
    line = fp.readline()
    words = line.split()
    natom = int(words[0])
    periodic = (words[1] == 'S' or words[1] == 's')
    fractional = (words[1] == 'F' or words[1] == 'f')
    periodic = (periodic or fractional)
    line = fp.readline()
    specienames = line.split()
    coords = np.empty((natom, 3), dtype=float)
    species = np.empty(natom, dtype=int)
    for ii in range(natom):
        line = fp.readline()
        words = line.split()
        species[ii] = int(words[1]) - 1
        coords[ii] = (float(words[2]), float(words[3]), float(words[4]))
    if periodic:
        line = fp.readline()
        origin = np.array([float(s) for s in line.split()], dtype=float)
        latvecs = np.empty((3, 3), dtype=float)
        for ii in range(3):
            line = fp.readline()
            latvecs[ii] = [float(s) for s in line.split()]
        if fractional:
            coords = frac2cart(latvecs, coords)
#         origin *= AA__BOHR
#         latvecs *= AA__BOHR
    else:
        origin = None
        latvecs = None
#     coords *= AA__BOHR
    return specienames, species, coords, origin, latvecs

def writegen(fname, data):
    """Writes the geometry as gen file. Units of length are assumed to be Angström.
       modified from calcderivs.py (DFTB+ developers group)
    
       fname: filename
       data:  tuple (specienames, species, coords, origin, latvecs)
    """

    fp = open(fname, "w")
    specienames, species, coords, origin, latvecs = data
    fp.write("%5d %s\n" % (len(coords), latvecs is None and "C" or "S"))
    fp.write(("%2s "*len(specienames) + "\n") % tuple(specienames))
#     coords = coords * BOHR__AA
    for ii in range(len(coords)):
        fp.write("%5d %5d %23.15E %23.15E %23.15E\n"
                 % (ii + 1, species[ii] + 1, coords[ii, 0], coords[ii, 1],
                    coords[ii, 2]))
    if latvecs is not None:
#         origin = origin * BOHR__AA
#         latvecs = latvecs * BOHR__AA
        fp.write("%23.15E %23.15E %23.15E\n" % tuple(origin))
        for ii in range(3):
            fp.write("%23.15E %23.15E %23.15E\n" % tuple(latvecs[ii]))
    fp.close()

def cart2frac(latvecs, coords):
    """Converts cartesian coordinates to fractional coordinates.
       unmodified from calcderivs.py (DFTB+ developers group)
    """
    invlatvecs = np.empty((3, 3), dtype=float)
    invlatvecs = np.transpose(latvecs)
    newcoords = np.array(coords)
    invlatvecs = np.linalg.inv(invlatvecs)
    for iat, atcoords in enumerate(coords):
        newcoords[iat] = np.dot(invlatvecs, atcoords)
    return newcoords


def frac2cart(latvecs, coords):
    """Converts fractional coordinates to cartesian ones.
       unmodified from calcderivs.py (DFTB+ developers group)
    """
    newcoords = np.array(coords)
    for iat, atcoords in enumerate(coords):
        newcoords[iat] = np.dot(np.transpose(latvecs), atcoords)
    return newcoords

def read_dftb_bands(fname):
    """Read band.out file and extract band energies and occupations at different k-points.
    
       fname: filename (usually band.out)
       
       returns: bands, band_occ
    """ 
    f = open(fname, "r")
    energies = []
    occs = []
    bands = []
    band_occ = []
    for x in f:
        line_vec = x.split()

        if len(line_vec) > 0:                
            if (line_vec[0] == 'KPT'):
                if len(energies) > 0:
                    bands.append(np.array(energies))
                    band_occ.append(np.array(occs))
                    energies = []
                    occs = []
            else:
                # print(line_vec[1])
                energies.append(float(line_vec[1]))
                occs.append(float(line_vec[2]))            

    bands.append(np.array(energies))
    band_occ.append(np.array(occs))
    f.close()
    
    return bands, band_occ

def read_detailed_out(fname):
    """Read detailed.out file and extract total energy in eV.
    
       fname: filename (usually band.out)
       
       returns: etot
    """     
    f = open(fname, "r")

    etot = -1.0    
    for x in f:
        
        if 'Total energy:' in x:
            line_vec = x.split()
            etot = float(line_vec[4])
            break

    f.close()
    
    return etot

def get_lumo(occs):
    i_lumo = 0
    for i, occ in enumerate(occs):
        if occ < 1.0:
            i_lumo = i
            break
    return i_lumo

def read_sqr(filename, nkpt = 1):
    """Read *sqr file and get hamiltonian/overlap matrix at the specified k-point.
    
       fname:   filename
       nkpt:    index of k-point, starting at 1
       
       returns: hamiltonian/overlap matrix (nallorbitals * nallorbitals)
    """         
    with open(filename) as f:
        # firt line is '#      REAL   NALLORB   NKPOINT'
        first_line = f.readline()
        # second line gives the respective values (the # in the beginning is for convenience)
        second_line = f.readline().strip('#')
    
    temp_dat = second_line.split()

    real_elem = temp_dat[0] == 'T'
    NALLORB = int(temp_dat[1])
    NKPOINT = int(temp_dat[2])
    
    if nkpt > NKPOINT:
        print("Requested k-point is not available (max. is %i)." % (NKPOINT))
        return None

    if real_elem:
        skip = 1
    else:
        skip = 2
    
    return np.loadtxt(filename, skiprows= 2 + nkpt*3 + (nkpt-1)*NALLORB , max_rows=NALLORB)[:,::skip]

def read_hamsqr1(nkpt = 1):
    """Read hamsqr1.dat in current directory and get hamiltonian matrix at the specified k-point.
    
       nkpt:    index of k-point, starting at 1
       
       returns: hamiltonian matrix (nallorbitals * nallorbitals)
    """    
    return read_sqr('hamsqr1.dat', nkpt)*HARTREE__EV

def read_oversqr(nkpt = 1):
    """Read oversqr.dat in current directory and get overlap matrix at the specified k-point.
    
       nkpt:    index of k-point, starting at 1
       
       returns: overlap matrix (nallorbitals * nallorbitals)
    """    
    return read_sqr('oversqr.dat', nkpt)

def _readline_ncolumns(lines, ncol):
    line = lines.pop(0).strip()
    columns = line.split()
    
    while len(columns) < ncol:
        line = lines.pop(0).strip()
        columns.extend(line.split())

    return columns

def _parse_def_line(line):
    columns = line.split(':')

    name = columns[0].strip()
    ftype = columns[1].strip()
    dims = int(columns[2])
    shape = [int(el) for el in columns[3].split(',')]
        
    return name, ftype, dims, shape
    
def read_tagfile(filen):
    file = open(filen, 'r')
    lines = file.readlines()

    tagdict = {}
    while len(lines) > 0:
        line = lines.pop(0).strip()
        name, ftype, dims, shape = _parse_def_line(line)
        
        if ftype == 'real':
            dtype = float
        elif ftype == 'integer':
            dtype = int
        
        linev = _readline_ncolumns(lines, np.prod(shape))
        vals = np.array([dtype(el) for el in linev]).reshape(np.flip(shape))
        tagdict[name] = vals

    return tagdict
