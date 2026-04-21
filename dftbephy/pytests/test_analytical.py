import numpy as np
import pytest
import dftbephy.analytical as ana


# -------------------------
# Model class sanity tests
# -------------------------

def test_parabolic2d_masking_and_sign_convention():
    gv, gs = 2, 2
    meff = 0.5
    eps_edge = 0.1
    mu = 0.2
    kBT = 0.05

    cb = ana.Parabolic2D(gv=gv, gs=gs, meff=meff, band_type=ana.Band.CONDUCTION, eps_edge=eps_edge)
    vb = ana.Parabolic2D(gv=gv, gs=gs, meff=meff, band_type=ana.Band.VALENCE, eps_edge=eps_edge)

    # Conduction: below edge => 0
    assert cb.density(eps_edge - 1e-6, mu, kBT) == 0.0
    # Valence: above edge => 0
    assert vb.density(eps_edge + 1e-6, mu, kBT) == 0.0

    # Valence density uses a negative sign in your implementation
    below = vb.density(eps_edge - 1e-3, mu, kBT)
    assert below <= 0.0


def test_linear2d_conductivity_is_nonnegative_for_finite_T():
    # conductivity integrand uses (-dFDde) which should be >= 0 at finite T
    model = ana.Linear2DConduction(gv=2, gs=2, vF=5.0)
    mu = 0.1
    kBT = 0.05

    eps = np.linspace(-1.0, 1.0, 21)
    vals = np.array([model.conductivity_crta(e, mu, kBT) for e in eps])
    assert np.all(vals >= 0.0)


def test_kane2d_masking():
    model = ana.Kane2D(gv=2, gs=2, meff=0.5, alpha=1.0, band_type=ana.Band.CONDUCTION, eps_edge=0.0)
    assert model.density(-1e-6, mu=0.0, kBT=0.1) == 0.0
    assert model.conductivity_crta(-1e-6, mu=0.0, kBT=0.1) == 0.0


def test_kane2d_reduces_to_parabolic_when_alpha_zero():
    gv, gs = 2, 2
    meff = 0.5
    eps_edge = 0.1
    mu = 0.2
    kBT = 0.05

    para = ana.Parabolic2D(
        gv=gv, gs=gs, meff=meff,
        band_type=ana.Band.CONDUCTION,
        eps_edge=eps_edge,
    )
    kane = ana.Kane2D(
        gv=gv, gs=gs, meff=meff, alpha=0.0,
        band_type=ana.Band.CONDUCTION,
        eps_edge=eps_edge,
    )

    eps_grid = np.linspace(eps_edge, eps_edge + 0.5, 101)

    para_density = np.array([para.density(eps, mu, kBT) for eps in eps_grid])
    kane_density = np.array([kane.density(eps, mu, kBT) for eps in eps_grid])

    para_cond = np.array([para.conductivity_crta(eps, mu, kBT) for eps in eps_grid])
    kane_cond = np.array([kane.conductivity_crta(eps, mu, kBT) for eps in eps_grid])

    assert np.allclose(kane_density, para_density, rtol=1e-12, atol=1e-14)
    assert np.allclose(kane_cond, para_cond, rtol=1e-12, atol=1e-14)


def test_band_models_zero_temperature_occupancy_behavior():
    gv, gs = 2, 2
    meff = 0.5
    eps_edge = 0.1
    kBT = 0.0

    # -------------------------
    # Conduction band
    # -------------------------
    cb = ana.Parabolic2D(
        gv=gv, gs=gs, meff=meff,
        band_type=ana.Band.CONDUCTION,
        eps_edge=eps_edge,
    )

    mu_cb = 0.3
    eps_below_mu = 0.2   # inside conduction band and occupied at T=0
    eps_at_mu = 0.3      # half-occupied by convention
    eps_above_mu = 0.4   # empty at T=0

    pref = (gs * gv) / (2.0 * np.pi) * meff

    assert cb.density(eps_below_mu, mu_cb, kBT) == pytest.approx(pref)
    assert cb.density(eps_at_mu, mu_cb, kBT) == pytest.approx(0.5 * pref)
    assert cb.density(eps_above_mu, mu_cb, kBT) == pytest.approx(0.0)

    # -------------------------
    # Valence band
    # -------------------------
    vb = ana.Parabolic2D(
        gv=gv, gs=gs, meff=meff,
        band_type=ana.Band.VALENCE,
        eps_edge=eps_edge,
    )

    mu_vb = 0.0
    eps_below_mu_v = -0.1   # occupied electron states => hole occupancy 0
    eps_at_mu_v = 0.0       # half-hole occupancy by convention
    eps_above_mu_v = 0.05   # above mu, but still inside valence band (<= eps_edge)

    assert vb.density(eps_below_mu_v, mu_vb, kBT) == pytest.approx(0.0)
    assert vb.density(eps_at_mu_v, mu_vb, kBT) == pytest.approx(-0.5 * pref)
    assert vb.density(eps_above_mu_v, mu_vb, kBT) == pytest.approx(-pref)


def test_parabolic2d_transport_integral_low_temperature_limit():
    gv, gs = 2, 2
    meff = 0.5
    eps_edge = 0.1
    mu = 0.3
    kBT = 0.005  # low but finite T

    model = ana.Parabolic2D(
        gv=gv, gs=gs, meff=meff,
        band_type=ana.Band.CONDUCTION,
        eps_edge=eps_edge,
    )

    eps = np.linspace(eps_edge, mu + 20.0 * kBT, 200001)
    integrand = np.array([model.conductivity_crta(e, mu, kBT) for e in eps])
    integral = np.trapezoid(integrand, eps)

    pref = (
        (ana.Constants.ELEMENTARY_CHARGE_Aps / (ana.Constants.HBAR_eV**2))
        * (gs * gv) / (2.0 * np.pi)
    )

    expected = pref * 2.0 * (mu - eps_edge)

    # small tolerance because finite-T broadening and finite integration window remain
    assert integral == pytest.approx(expected, rel=2e-3, abs=1e-8)


def test_kane2d_conductivity_large_alpha_asymptotic_limit():
    gv, gs = 2, 2
    meff = 0.5
    eps_edge = 0.1
    mu = 0.3
    kBT = 0.05
    alpha = 1.0e6  # effectively large-alpha for the chosen x range

    model = ana.Kane2D(
        gv=gv, gs=gs, meff=meff, alpha=alpha,
        band_type=ana.Band.CONDUCTION,
        eps_edge=eps_edge,
    )
    # away from the band edge so that alpha*x >> 1 really holds
    eps_grid = np.linspace(eps_edge + 1e-2, eps_edge + 0.5, 101)
    model_vals = np.array([model.conductivity_crta(eps, mu, kBT) for eps in eps_grid])

    pref = (
        (ana.Constants.ELEMENTARY_CHARGE_Aps / (ana.Constants.HBAR_eV**2))
        * (gs * gv) / (2.0 * np.pi)
    )

    x = eps_grid - eps_edge
    expected = pref * x * (-ana.dfermi_de(eps_grid, mu, kBT))

    assert np.allclose(model_vals, expected, rtol=1e-4, atol=1e-10)

def test_kane2d_density_large_alpha_has_linear_in_x_scaling_after_rescaling():
    gv, gs = 2, 2
    meff = 0.5
    eps_edge = 0.1
    mu = 0.3
    kBT = 0.05
    alpha = 1.0e6

    model = ana.Kane2D(
        gv=gv, gs=gs, meff=meff, alpha=alpha,
        band_type=ana.Band.CONDUCTION,
        eps_edge=eps_edge,
    )
    # away from the band edge so that alpha*x >> 1 really holds
    eps_grid = np.linspace(eps_edge + 1e-2, eps_edge + 0.5, 101)
    model_vals = np.array([model.density(eps, mu, kBT) for eps in eps_grid])

    pref = (gs * gv) / (2.0 * np.pi) * meff
    x = eps_grid - eps_edge
    expected = pref * (2.0 * x) * ana.fermi_eps(eps_grid, mu, kBT)

    # divide by alpha before comparing asymptotic shape
    # Since density ~ pref * (1 + 2 alpha x) * f, dividing by alpha gives
    # pref * ((1/alpha) + 2x) * f -> pref * 2x * f in the large-alpha limit
    assert np.allclose(model_vals / alpha, expected, rtol=1e-4, atol=1e-10)

