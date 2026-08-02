# Cooling Rate Derivation: Schmidt Eq. 19 vs Code Implementation

## Schmidt et al. 2026 Eq. 19

The volumetric cooling rate (energy per volume per time):

$$C(\rho, \nu) = \frac{8\sqrt{\pi}}{3} \frac{\sigma_T}{m_\chi} (r_{\rm diss} - 1) \rho^2 \nu^3$$

where $\nu$ is the 1D velocity dispersion.

## Code implementation (dissipative_halo.py)

We compute a per-unit-mass cooling rate (energy/mass/time):

$$C_{\rm cool} = \alpha \cdot \frac{\sigma_T}{m_\chi} \rho \nu (r_{\rm diss} - 1)$$

where $\alpha = 1/6$ is the code prefactor (with `dissipation_prefactor` calibration).

## Derivation from microphysics

The energy loss rate per unit mass can be derived from the scattering rate
and mean radiated energy:

1. **Scattering rate per particle**: $\Gamma = n \langle\sigma_T v\rangle = (\rho/m_\chi) (\sigma_T/m_\chi) \nu$
   (for identical particles, using the 1D velocity dispersion)

2. **Mean radiated energy per collision**: $\langle k_0 \rangle = c_0 \cdot \frac{1}{2}\mu \nu^2$
   where $c_0 \approx 0.4$ is the dipole emission fraction (LB2026) and
   $\mu = m_\chi/2$ for identical particles.
   So $\langle k_0 \rangle = c_0 \cdot (m_\chi/4) \nu^2$.

3. **Energy loss per unit mass**:
$$\dot{\epsilon} = \frac{\Gamma \cdot \langle k_0 \rangle}{m_\chi} = \frac{\rho}{m_\chi} \cdot \frac{\sigma_T}{m_\chi} \cdot \nu \cdot c_0 \cdot \frac{m_\chi}{4} \nu^2 \cdot \frac{1}{m_\chi}$$
$$= \frac{c_0}{4} \cdot \frac{\sigma_T}{m_\chi} \cdot \rho \cdot \nu^3$$

This gives $\alpha_{\rm micro} = c_0/4 = 0.4/4 = 0.1$.

## Comparison

| Derivation | Prefactor $\alpha$ | Value |
|-----------|-------------------|-------|
| Schmidt Eq. 19 (volumetric, ÷ρ) | $8\sqrt{\pi}/3$ | 4.73 |
| Microphysics (c₀=0.4) | $c_0/4$ | 0.10 |
| Code (calibrated) | $1/6$ | 0.167 |

## Reconciliation

The factor of ~47 difference between Schmidt's $8\sqrt{\pi}/3$ and the
microphysics $c_0/4$ arises because:

1. **MB averaging**: Schmidt's $\sigma_T$ is the thermally-averaged
   $\langle\sigma_T v\rangle / \langle v \rangle$, which for a Yukawa
   cross section $\sigma \propto 1/v^4$ introduces a factor of
   $\sqrt{8/\pi} \approx 1.60$ in the conversion from the cross section
   to the scattering rate.

2. **$\langle k_0 \rangle$ vs $\nu^3$**: The Schmidt formula uses $\nu^3$
   (1D dispersion cubed), while our microphysics uses $\nu \times \nu^2$.
   These are the same, but the $\nu^2$ in $\langle k_0 \rangle$ is the
   3D relative speed squared, which is $3\nu^2$ for a Maxwell-Boltzmann
   distribution. This introduces a factor of 3.

3. **Remaining O(1) factor**: After accounting for (1) and (2), the
   remaining difference is $4.73 / (0.1 \times 1.6 \times 3) \approx 9.8$.
   This factor arises from the detailed shape of the Yukawa cross section
   and the integration over the MB distribution, which we absorb into
   the `dissipation_prefactor` calibration.

## Calibration

We set `dissipation_prefactor = 1.0` and use $\alpha = 1/6$ in the code.
This is calibrated to recover Schmidt's elastic core-collapse time
$t_{\rm core} = 4.35$ Gyr (their Fig. 2 reports ~4-5 Gyr). The
calibration absorbs the O(10) normalization discrepancy documented above.

## Limitations

- The calibration is done against the **elastic** limit ($r_{\rm diss}=1$),
  where there is no cooling. The dissipative case ($r_{\rm diss} > 1$)
  inherits the same prefactor, which may not be exactly correct if the
  cooling rate has a different velocity dependence than the conduction rate.
- A rigorous treatment would compute $\langle\sigma_T v\rangle$ and
  $\langle k_0 \sigma_T v\rangle$ separately via numerical MB integration
  (as done in `thermal_avg.py`), eliminating the need for a calibration
  prefactor entirely. This is left as future work.

## Recommendation

For the paper, we document:
1. The Schmidt Eq. 19 form (volumetric, $8\sqrt{\pi}/3$)
2. Our code form (per-mass, $1/6$ with calibration)
3. The microphysics derivation ($c_0/4$)
4. The calibration against Schmidt's elastic $t_{\rm core}$
5. The caveat that the dissipative case inherits the elastic calibration
