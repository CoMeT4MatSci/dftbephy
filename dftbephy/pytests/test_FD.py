import numpy as np
import pytest
import dftbephy.analytical as ana

# -------------------------
# Constants / helpers
# -------------------------
FD_XMAX = 18.42068


def test_fd_xmax_value_is_expected():
    # sanity check that module constant didn't change
    assert pytest.approx(ana._FD_XMAX, rel=0.0, abs=1e-12) == FD_XMAX


# -------------------------
# fermi tests (x-space)
# -------------------------
def test_fermi_basic_values_and_bounds():
    x = np.array([-100.0, -10.0, 0.0, 10.0, 100.0])
    f = ana.fermi(x)

    assert isinstance(f, np.ndarray)
    assert f.shape == x.shape
    assert np.all((f >= 0.0) & (f <= 1.0))
    assert pytest.approx(f[x == 0.0][0], abs=1e-15) == 0.5


def test_fermi_tail_clipping_is_exact():
    # Beyond +FD_XMAX -> exactly 0; below -FD_XMAX -> exactly 1
    x = np.array([-FD_XMAX - 1e-6, -FD_XMAX, 0.0, FD_XMAX, FD_XMAX + 1e-6])
    f = ana.fermi(x)

    # boundaries are outside the strict "mid" mask, so they are clipped
    assert f[0] == 1.0
    assert f[1] == 1.0
    assert pytest.approx(f[2], abs=1e-15) == 0.5
    assert f[3] == 0.0
    assert f[4] == 0.0


def test_fermi_particle_hole_symmetry():
    # For any x: f(x) + f(-x) = 1 (logistic identity; clipping preserves it)
    x = np.linspace(-30.0, 30.0, 301)
    f = ana.fermi(x)
    g = ana.fermi(-x)
    assert np.allclose(f + g, 1.0, atol=1e-14)


# -------------------------
# dfermi_deps tests (x-space derivative df/dx)
# -------------------------
def test_dfermi_deps_sign_and_peak_at_zero():
    x = np.array([-1.0, -0.1, 0.0, 0.1, 1.0])
    df = ana.dfermi_deps(x)

    # df/dx should be <= 0 everywhere (monotone decreasing f)
    assert np.all(df <= 0.0)

    # magnitude peaked at x=0 (symmetric)
    assert abs(df[2]) >= abs(df[1])
    assert abs(df[2]) >= abs(df[3])

    # exact known value: df/dx at x=0 is -1/4
    assert pytest.approx(df[2], abs=1e-15) == -0.25


def test_dfermi_deps_tail_is_exact_zero():
    x = np.array([-FD_XMAX - 1e-3, -FD_XMAX, FD_XMAX, FD_XMAX + 1e-3])
    df = ana.dfermi_deps(x)
    assert np.all(df == 0.0)


def test_dfermi_deps_matches_numeric_derivative_in_mid_region():
    # Avoid clipping region: stay safely inside (-FD_XMAX, FD_XMAX)
    x = np.linspace(-5.0, 5.0, 101)
    h = 1e-7
    num = (ana.fermi(x + h) - ana.fermi(x - h)) / (2.0 * h)
    ana_df = ana.dfermi_deps(x)
    assert np.allclose(ana_df, num, rtol=1e-7, atol=1e-9)


# -------------------------
# fermi_eps tests (energy space)
# -------------------------
def test_fermi_eps_zero_temperature_step_convention():
    mu = 0.1
    kBT = 0.0
    eps = np.array([mu - 1e-3, mu, mu + 1e-3])
    f = ana.fermi_eps(eps, mu, kBT)

    assert f[0] == 1.0
    assert f[2] == 0.0
    assert f[1] == 0.5  # convention


def test_fermi_eps_is_fermi_of_dimensionless_argument():
    mu = 0.2
    kBT = 0.05
    eps = np.linspace(mu - 0.4, mu + 0.4, 41)

    x = (eps - mu) / kBT
    f_eps = ana.fermi_eps(eps, mu, kBT)
    f_x = ana.fermi(x)

    assert np.allclose(f_eps, f_x, atol=0.0, rtol=0.0)


def test_fermi_eps_tail_clipping_energy_window():
    mu = 0.0
    kBT = 0.05
    # Clipping in x-space at +/- FD_XMAX corresponds to eps-mu at +/- FD_XMAX*kBT
    dE = FD_XMAX * kBT
    eps = np.array([mu - dE - 1e-6, mu - dE, mu, mu + dE, mu + dE + 1e-6])

    f = ana.fermi_eps(eps, mu, kBT)
    assert f[0] == 1.0
    assert f[1] == 1.0
    assert pytest.approx(f[2], abs=1e-15) == 0.5
    assert f[3] == 0.0
    assert f[4] == 0.0


# -------------------------
# dfermi_de tests (energy derivative df/deps + "units")
# -------------------------
def test_dfermi_de_zero_temperature_is_zero_everywhere():
    mu = 0.0
    kBT = 0.0
    eps = np.array([-1.0, 0.0, 1.0])
    df = ana.dfermi_de(eps, mu, kBT)
    assert np.all(df == 0.0)


def test_dfermi_de_chain_rule_relation_to_dfermi_deps():
    # df/de = (df/dx)/kBT where x=(eps-mu)/kBT
    mu = 0.1
    kBT = 0.07
    eps = np.linspace(mu - 0.3, mu + 0.3, 61)

    x = (eps - mu) / kBT
    dfde = ana.dfermi_de(eps, mu, kBT)
    dfdx = ana.dfermi_deps(x)

    assert np.allclose(dfde, dfdx / kBT, rtol=0.0, atol=0.0)


def test_dfermi_de_peak_value_has_correct_energy_units_scale():
    # At eps=mu: df/dx = -1/4, so df/de = -(1/4)/kBT
    mu = 0.0
    kBT = 0.05
    eps = np.array([mu])
    dfde = ana.dfermi_de(eps, mu, kBT)[0]

    expected = -0.25 / kBT  # units: 1/energy (eV^-1 if kBT in eV)
    assert pytest.approx(dfde, rel=0.0, abs=1e-12) == expected


def test_dfermi_de_integrates_to_minus_one_over_energy():
    # For the ideal logistic without clipping, ∫ df/de dE = f(+∞)-f(-∞) = -1.
    # With clipping at +/- FD_XMAX, this should be extremely close to -1.
    mu = 0.0
    kBT = 0.05

    Emax = FD_XMAX * kBT
    eps = np.linspace(mu - Emax, mu + Emax, 200001)  # fine grid
    dfde = ana.dfermi_de(eps, mu, kBT)

    integral = np.trapezoid(dfde, eps)
    assert pytest.approx(integral, abs=5e-6) == -1.0


def test_minus_dfermi_de_is_positive_and_normalized():
    # -df/de is the positive "thermal broadening kernel" used in transport
    mu = 0.0
    kBT = 0.05

    Emax = FD_XMAX * kBT
    eps = np.linspace(mu - Emax, mu + Emax, 200001)

    w = -ana.dfermi_de(eps, mu, kBT)

    assert np.all(w >= 0.0)
    assert pytest.approx(np.trapezoid(w, eps), abs=5e-6) == 1.0


def test_dfermi_de_tail_is_exact_zero_outside_window():
    mu = 0.0
    kBT = 0.05
    Emax = FD_XMAX * kBT
    eps = np.array([mu - Emax - 1e-6, mu + Emax + 1e-6])

    dfde = ana.dfermi_de(eps, mu, kBT)
    assert np.all(dfde == 0.0)
