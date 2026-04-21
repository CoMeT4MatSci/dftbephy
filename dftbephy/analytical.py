from dataclasses import dataclass
from enum import Enum
import numpy as np

################################################################################
# Functions which might be moved to the other modules.
################################################################################

# --------------------------------------------------------------------
# Constants --> units.py
# --------------------------------------------------------------------
@dataclass
class Constants:
    HBAR_eV: float = 6.582e-4                      # [eV * picosecond]
    ELEMENTARY_CHARGE_Aps: float = 1.602176634e-7  # [Ampere * picosecond]
    BOLTZMANN_eV: float = 8.617333262e-5           # [eV / K]
    # this won't be in the Constants class
    qe_SI: float = ELEMENTARY_CHARGE_Aps * 1e-12   # Electron Charge [C]


# --------------------------------------------------------------------
# Fermi-Dirac --> analysis.py
# --------------------------------------------------------------------

###################   Functions from dftbephy-mpi.py  ################

# _FD_XMAX = np.log(np.finfo(float).max)/2
_FD_XMAX = 18.42068 # from BoltzTrap

def fermi(x):
    """
    Fermi–Dirac occupation (unitless)
    -----------------------
        f(x = 1 / (exp((x)/kBT) + 1)

    """
    f = np.where(x < 0., 1., 0.)
    norm_ind = np.logical_and(x > -_FD_XMAX, x < _FD_XMAX)
    f[norm_ind] = 1/(np.exp(x[norm_ind]) + 1.)
    return f

def dfermi_deps(x):
    f = np.zeros_like(x)
    norm_ind = np.logical_and(x > -_FD_XMAX, x < _FD_XMAX)
    f[norm_ind] = -np.exp(x[norm_ind])/(np.exp(x[norm_ind]) + 1.)**2
    return f


def step(x):
    f = np.where(x < 0., 1., 0.)
    return f

#################  end dftbephy-mpi.py functions  #####################

#######################  fermi func. for the analytical models ########
def fermi_eps(eps, mu, kBT):
    """
    Fermi–Dirac occupation
    -----------------------
        f(eps; mu, kBT) = 1 / (exp((eps - mu)/kBT) + 1)

    """
    eps = np.asarray(eps, dtype=float)

    if kBT == 0.0:
        delta = eps - mu
        f = np.where(delta < 0.0, 1.0, 0.0)
        f[np.isclose(delta, 0.0)] = 0.5
        return f

    x = (eps - mu) / kBT
    return fermi(x)


def dfermi_de(eps, mu, kBT):
    """
    Derivative of Fermi–Dirac with respect to energy, df/deps.
    Relation to df/dx
    ------------------------------
        x = (eps - mu) / kBT
        df/deps = (1/kBT) * df/dx
    """
    eps = np.asarray(eps, dtype=float)

    if kBT == 0.0:
        return np.zeros_like(eps, dtype=float)

    x = (eps - mu) / kBT
    return  dfermi_deps(x) / kBT

##############  end of fermi func. for the analytical models ##########

################################################################################
# End of functions which might be moved to the other modules.
################################################################################


class Band(Enum):
    CONDUCTION = "conduction"
    VALENCE = "valence"

@dataclass  # with keyword argument, ordering is not important, Version >= Python 3.7
class Linear2DConduction:
    """
    2D Dirac / linear dispersion model (graphene-like)

    Parameters:
      vF: Fermi velocity [eV * Å]
      gv, gs: valley and spin degeneracy
    """
    gv: int
    gs: int
    vF: float
#    constants: Constants = Constants() # keep it as defined here

    def density(self, eps: float, mu: float, kBT: float) -> float:
        dos = (self.gs * self.gv) / (2.0 * np.pi) * (np.abs(eps) / (self.vF**2)) * fermi_eps(eps, mu, kBT)
        return dos

    def conductivity_crta(self, eps: float, mu: float, kBT: float) -> float:
        vvdos = (Constants.ELEMENTARY_CHARGE_Aps / (Constants.HBAR_eV**2)) * (self.gs * self.gv) / (2.0 * np.pi) * np.abs(eps) * (-dfermi_de(eps, mu, kBT))
        return vvdos

    def conductivity_serta(self, eps: float, mu: float, kBT: float, dTaudEps: float) -> float:
        # to avoid dividing by zero
        eps_abs = max(float(np.abs(eps)), 1e-30)
        vvdos = (Constants.ELEMENTARY_CHARGE_Aps / (Constants.HBAR_eV**2)) * (self.gs * self.gv) / (2.0 * np.pi) * eps_abs * (-dfermi_de(eps, mu, kBT)) * (1.0 / (eps_abs * dTaudEps))
        return vvdos


@dataclass
class Parabolic2D:
    """
    Parabolic band model

    Parameters:
      m_eff:    effective mass (includes HBAR_eV^2 -> (d^2E/dk^2)^-1). The unit is not in electron rest mass [m0]
      eps_edge: band edge energy [eV]. Conduction band minimum if CONDUCTION, Valence band maximum if VALENCE
      band:     CONDUCTION or VALENCE
      gv, gs:   valley and spin degeneracies
    """
    gv: int
    gs: int
    meff: float
    band_type: Band
    eps_edge: float
#    constants: Constants = Constants()  # keep it as defined here


    def _mask(self, eps: float) -> bool:
        if self.band_type is Band.CONDUCTION:
            return eps >= self.eps_edge
        elif self.band_type is Band.VALENCE:
            return eps <= self.eps_edge
        raise ValueError(f"Unknown band type: {self.band_type}")

    def _x(self, eps: float) -> float:
        if self.band_type is Band.CONDUCTION:
            return eps - self.eps_edge
        elif self.band_type is Band.VALENCE:
            return self.eps_edge - eps
        raise ValueError(f"Unknown band type: {self.band_type}")

    def _occ(self, eps: float, mu: float, kBT: float) -> float:
        fd = fermi_eps(eps, mu, kBT)
        if self.band_type is Band.CONDUCTION:
            return fd
        elif self.band_type is Band.VALENCE:
            return 1.0 - fd
        raise ValueError(f"Unknown band type: {self.band_type}")

    def density(self, eps: float, mu: float, kBT: float) -> float:
        if not self._mask(eps):
            return 0.0

        pref = (self.gs * self.gv) / (2.0 * np.pi) * self.meff

        sign = -1.0 if self.band_type is Band.VALENCE else 1.0
        dos = sign * pref * self._occ(eps, mu, kBT)
        return dos

    def conductivity_crta(self, eps: float, mu: float, kBT: float) -> float:
        if not self._mask(eps):
            return 0.0

        pref = (Constants.ELEMENTARY_CHARGE_Aps / (Constants.HBAR_eV**2)) * (self.gs * self.gv) / (2.0 * np.pi)
        vvdos = pref * 2.0 * self._x(eps) * (-dfermi_de(eps, mu, kBT))
        return vvdos


@dataclass
class Kane2D:
    """
    2D Kane (nonparabolic) band model

    Parameters:
      m_eff, eps_edge, band, gv are gs as explained in Parabolic2D class.
      alpha: Kane nonparabolicity parameter [1/eV]
    """
    gv: int
    gs: int
    meff: float
    alpha: float
    band_type: Band
    eps_edge: float
#    constants: Constants = Constants()    # keep it as defined here

    def _mask(self, eps: float) -> bool:
        if self.band_type is Band.CONDUCTION:
            return eps >= self.eps_edge
        elif self.band_type is Band.VALENCE:
            return eps <= self.eps_edge
        raise ValueError(f"Unknown band type: {self.band_type}")

    def _x(self, eps: float) -> float:
        if self.band_type is Band.CONDUCTION:
            return eps - self.eps_edge
        elif self.band_type is Band.VALENCE:
            return self.eps_edge - eps
        raise ValueError(f"Unknown band type: {self.band_type}")

    def _occ(self, eps: float, mu: float, kBT: float) -> float:
        fd = fermi_eps(eps, mu, kBT)
        if self.band_type is Band.CONDUCTION:
            return fd
        elif self.band_type is Band.VALENCE:
            return 1.0 - fd
        raise ValueError(f"Unknown band type: {self.band_type}")

    def density(self, eps: float, mu: float, kBT: float) -> float:
        if not self._mask(eps):
            return 0.0

        x = self._x(eps)
        pref = (self.gs * self.gv) / (2.0 * np.pi) * self.meff

        sign = -1.0 if self.band_type is Band.VALENCE else 1.0
        dos = sign * pref * (1.0 + 2.0 * self.alpha * x) * self._occ(eps, mu, kBT)
        return dos

    def conductivity_crta(self, eps: float, mu: float, kBT: float) -> float:
        if not self._mask(eps):
            return 0.0

        x = self._x(eps)
        if self.band_type is Band.CONDUCTION:
            kaneE = x + self.alpha * x**2
        elif self.band_type is Band.VALENCE:
            kaneE = x - self.alpha * x**2
        else:
            raise ValueError(f"Unknown band_type: {self.band_type}")
        pref = (Constants.ELEMENTARY_CHARGE_Aps / (Constants.HBAR_eV**2)) * (self.gs * self.gv) / (2.0 * np.pi)
        vvdos = pref * 2.0 * (kaneE / (1.0 + 2.0 * self.alpha * x)) * (-dfermi_de(eps, mu, kBT))
        return vvdos
