# Band Non-parabolicity

The conduction and valence bands in semiconductors are commonly approximated as parabolic near the band edges. In practice, however, the band structure of many materials departs from this simple parabolic description, particularly at high carrier concentrations {cite}`Lundstrom2000,Walsh2019`. Non-parabolicity in the electronic dispersion is often described using the Kane model {cite}`Maassen2020`,

```{math}
:label: Eq:Kane
\varepsilon + \alpha\varepsilon^2 = \frac{\hbar^2 \vec{k}^2}{2m^*},
```

where $m^*$ is the effective mass and $\alpha$ quantifies the degree of non-parabolicity which differs for each band. Larger $\alpha$ corresponds to stronger deviations from a parabolic dispersion and increased band flattening. In the Kane model, this band flattening enhances the density of states by a factor of $(1+2\alpha\varepsilon)$ compared to the parabolic approximation. For the charge carrier densities, we obtain

```{math}
:label: Eq:N_C-Kane
n_c = g_s g_v \frac{m^*}{2\pi \hbar^2} \int d\varepsilon (1+2\alpha\varepsilon)[f^0(\varepsilon; \mu,T) -f^0(\varepsilon;\varepsilon_F,0)],
```

where $g_s$ and $g_v$ are spin and valley degeneracies, repectively. Depending on the high-symmetry point where the band edge resides, valley degeneracies must be taken into account.
The non-parabolicity also changes the energy dependence of the conductivity integrands,

```{math}
:label: Eq:sigma-Kane
\sigma_{2D} = g_s g_v \frac{e^2}{2 \pi \hbar^2} \int d\varepsilon \left[-\frac{\partial f^0(\varepsilon(\vec{k}); \mu, T)}{\partial \varepsilon}\right] \frac{\varepsilon + \alpha \varepsilon^2}{1 + 2\alpha\varepsilon } \tau_{n}(\vec{k}),
```
where $\tau$ is the $k$-dependent relaxation time. Here, it is restricted to a single band. If multiple bands are present, their contributions must be summed. The 2D conductivity is obtained as the average of the in-plane components of the conductivity tensor, $(\sigma_{xx} + \sigma_{yy})/2$. Details of the derivation can be found in the supporting information file of our paper {cite}`Unsal2025`.

In the parabolic band limit, both the charge carrier density and the electrical conductivity scale linearly with temperature in 2D when the relaxation time $\tau$ is taken to be constant. Under these conditions, their ratio, i.e. the mobility, is independent of both temperature and carrier concentration and can be obtained analytically from the standard 2D density of states and the corresponding Fermi–Dirac integrals for a parabolic dispersion, yielding $\mu_{\text{2D,parabolic}} = e^2 \tau / m$ {cite}`Hicks-Dresselhaus-2D`.

In contrast, within the Kane-band model and constant $\tau$, any temperature dependence of the conductivity arises only from the non-parabolic band structure, which renders the conductivity integrand non-linear in temperature. As a result, the mobility becomes intrinsically temperature dependent even in the absence of temperature-dependent scattering and this makes the Kane-band model particularly suitable for analyzing the temperature dependence of carrier mobility under the constant relaxation time approximation.

Furthermore, for large-scale systems, calculating the mobility from first-principles methods is computationally demanding. Evaluating its temperature- or carrier-concentration dependence requires repeating these calculations over different temperatures and doping levels which further increases the computational cost. In this context, analytical models like the Kane-band model, provide a valuable starting point for gaining insight into the temperature and carrier concentration dependency of the mobility before undertaking more expensive first-principles calculations.


*References*
```{bibliography}
:style: unsrtalpha
:filter: docname in docnames
```
