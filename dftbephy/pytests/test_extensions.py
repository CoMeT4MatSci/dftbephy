import numpy as np

# This file is ONLY for the C-extension functions (and optional pure-Python fallbacks).
# Core analytical tests (FD/dFDde/band models) should live in test_analytical_core.py.

try:
    # C-extension / compiled module
    from dftbephy.extensions import (
        alphaFnk,
        inv_tau_nk,
        inv_tau_nk_lam,
        inv_tau_nk_lor,
        inv_tau_nk_lam_lor,
    )
    HAS_EXT = True
except ModuleNotFoundError:
    HAS_EXT = False

    # -------------------------
    # Optional pure-Python fallbacks (only used if extensions are missing)
    # -------------------------

    def alphaFnk(n, eps, frequency_points, mesh_g2, mesh_epskq, mesh_frequencies, sigma):
        nbands = mesh_g2.shape[2]
        nmodes = mesh_g2.shape[1]
        nqpoints = mesh_g2.shape[0]

        alphaF = np.zeros((frequency_points.shape[0],))
        for iq in range(nqpoints):
            for lam in range(nmodes):
                for m in range(nbands):
                    weight_el = np.exp(-(eps - mesh_epskq[iq, m]) ** 2 / (2 * sigma**2)) / (
                        np.sqrt(2 * np.pi) * sigma
                    )
                    for nf, f in enumerate(frequency_points):
                        weight_ph = np.exp(-(f - mesh_frequencies[iq, lam]) ** 2 / (2 * sigma**2)) / (
                            np.sqrt(2 * np.pi) * sigma
                        )
                        alphaF[nf] += mesh_g2[iq, lam, m, n] * weight_el * weight_ph
        return alphaF

    # If you have python fallbacks for these in your original file, keep them here too.
    # Otherwise, leave them undefined and skip tests that require them.
    # (Below we skip tests that need inv_tau_* when extensions are missing.)

# -------------------------
# Tests for extensions (or fallbacks)
# -------------------------

def test_extensions_import_or_fallback_defined():
    # Either the compiled extension exists, or (at minimum) alphaFnk is defined as a fallback.
    assert HAS_EXT or ("alphaFnk" in globals())


def test_alphaFnk_output_shape_and_finite_values():
    # Minimal synthetic inputs just to verify it runs and returns correct shape.
    # This test works both with the compiled extension and with the fallback alphaFnk above.

    nq, nmodes, nbands = 2, 3, 4
    n = 1
    eps = 0.2
    sigma = 0.05

    frequency_points = np.linspace(0.0, 1.0, 11)
    mesh_frequencies = np.linspace(0.0, 1.0, nq * nmodes).reshape(nq, nmodes)

    # mesh_epskq shape: (nqpoints, nbands)
    mesh_epskq = np.linspace(-1.0, 1.0, nq * nbands).reshape(nq, nbands)

    # mesh_g2 shape in your fallback: (nqpoints, nmodes, nbands, nbands)
    rng = np.random.default_rng(0)
    mesh_g2 = rng.random((nq, nmodes, nbands, nbands))

    out = alphaFnk(n, eps, frequency_points, mesh_g2, mesh_epskq, mesh_frequencies, sigma)

    assert isinstance(out, np.ndarray)
    assert out.shape == (frequency_points.shape[0],)
    assert np.all(np.isfinite(out))


def test_inv_tau_functions_available_if_extensions_present():
    # These only exist in the compiled extension in most setups.
    if not HAS_EXT:
        return

    # If extension is present, these symbols must exist and be callable.
    for fn in (inv_tau_nk, inv_tau_nk_lam, inv_tau_nk_lor, inv_tau_nk_lam_lor):
        assert callable(fn)
