# P1: Derivation Notes — σ_T(v) and r_diss(v) from LB2026

## Source
Lankester-Broche & Pradler 2026, JCAP 03, 034 (arXiv:2509.12317)
"Towards a theory of dissipative Dark Matter I: the Born limit"

## Key equations extracted

### Elastic transfer cross section σ_T — Eq. (2.6)

For identical DM particles (a₁ = a₂ = a, μ = m_χ/2):

$$\sigma_T = \frac{a^4}{8\pi \mu^2 v_i^4}\left[\ln\left(1 + \frac{4\mu^2 v_i^2}{m_{\phi,V}^2}\right) - \frac{4\mu^2 v_i^2}{4\mu^2 v_i^2 + m_{\phi,V}^2}\right]$$

where:
- a = √(4π α_D) for vector coupling (gauge), or Yukawa/trilinear for scalar
- m_{φ,V} = mediator mass (same as emitted particle for long-range)
- v_i = relative velocity

**Limiting behaviors:**
- Long-range (m ≪ μv): Rutherford, σ_T ~ a⁴/(μ²v⁴) × ln(4μ²v²/m²)
- Short-range (m' ≫ μv): σ_T ~ a⁴/(μ² m'⁴)

### Energy loss rate — Eq. (5.6), long-range, massless

For **identical particles** (dipole vanishes; quadrupole dominates):

$$\dot{\epsilon}_{\text{quad}} = \frac{n^2}{f_{\phi,V}} \frac{F_Q^2}{6\pi^{7/2}\sqrt{2}} \cdot \frac{24}{45} (\mu T)^{3/2} (4\alpha + \beta)$$

where:
- f_φ = 1 (scalar), f_V = 2 (vector)
- F_Q = quadrupole emission factor (model-dependent, ~ a² × m for identical)
- α, β = O(1) model-dependent coefficients from Sec 3.5
- T = temperature = ν² (1D velocity dispersion²)

**Parametric scaling (massless, long-range):**
- ε̇ ~ n² F² T^{3/2} ~ n² a⁴ m² v³ (since T ~ μv²)
- Scatter rate per particle ~ n σ_T v ~ n a⁴/(μ²v⁴) × v = n a⁴/(μ²v³)
- <k₀> per collision ~ 0.4 × μv²/2 (fixed fraction)
- **r_diss - 1 = ε̇ / (n² σ_T v <k₀>) ~ const** ← recovers Schmidt et al. ansatz

### Massive emission corrections — Eq. (3.26), (3.27)

For emitted particle of mass m_{φ,V} > 0:

$$C_\phi = 1 - \frac{m_\phi^2}{\omega^2} \quad \text{(scalar)}$$
$$C_V = 1 + \frac{m_V^2}{2\omega^2} \quad \text{(vector, incl. longitudinal)}$$

These multiply the squared amplitude ⟨|M_D|²⟩.

**Symmetry-breaking velocity scale:**
- Emission kinematically allowed when ω ≥ m_{φ,V}
- Typical ω ~ μv²/2, so threshold at **v* = √(2m_{φ,V}/μ)**
- Below v*: Boltzmann suppression exp(-m/T) = exp(-m/(μv²/2))
- Above v*: emission efficient, r_diss → const (massless limit)

### r_diss(v) — derived

$$r_{\text{diss}}(v) - 1 = (r_{\text{diss}} - 1)\big|_{m=0} \times \exp\left(-\frac{m_{\phi,V}}{T(v)}\right)$$

where T(v) = μv²/2 (kinetic temperature at relative velocity v).

**This is the key result:**
1. **Massless emission**: r_diss = const → Schmidt et al. rescaling symmetry holds
2. **Massive emission**: r_diss(v) is velocity-dependent, symmetry broken
3. The breaking scale v* = √(2m/μ) is set by the mediator mass
4. For m_χ = 10 GeV, v* in astrophysical range requires m in keV scale

## Benchmark calibration

For m_χ = 10 GeV, μ = 5 GeV:
| v* [km/s] | m [keV] | m [GeV] |
|-----------|---------|---------|
| 100 | 0.28 | 2.78e-7 |
| 200 | 1.11 | 1.11e-6 |
| 500 | 6.95 | 6.95e-6 |
| 1000 | 27.8 | 2.78e-5 |

## Three benchmark models (P1)

| Model | Emitted particle | m_mediator | v* [km/s] | Symmetry breaking |
|-------|-----------------|-----------|-----------|-------------------|
| M1 | Vector V (dark photon) | 1.1 keV | ~200 | Dwarf-MW transition |
| M2 | Scalar φ | 6.95 keV | ~500 | MW-cluster transition |
| M3 | Vector V (massless) | 0 | — | None (control) |

## Key physics

- **M1 vs M2**: same α_D, m_χ; different mediator type (vector vs scalar) and mass
  → different v* → different r_diss(v) shapes
- **M1 vs M3**: same model (chi-V), different m_V
  → isolates the effect of mediator mass alone
- **M3 control**: validates that massless emission recovers Schmidt et al.
  constant-r_diss ansatz (rescaling symmetry intact)

## Limitations of this P1 implementation

1. **r_diss amplitude C=0.05 is fiducial** — requires calibration to LB2026
   full expressions (Sec 5.1.2 quadrupole, model-dependent α, β)
2. **Quadrupole channel only** (identical particles); dipole for
   distinguishable two-component DM not yet implemented
3. **Boltzmann approximation** for massive emission; full numerical
   integration of Eq. (3.28) with C_φ/C_V factors deferred to P2
4. **Born regime only**; beyond-Born (LB2026 Part II) not yet available

## Next steps (P2)

- Implement full quadrupole F_Q, α, β from LB2026 Sec 3.5
- Numerical integration of energy-differential cross section (Eq 3.28 + C_φ/V)
- Thermal averages ⟨σ_T v⟩, ⟨r_diss σ_T v⟩ for GravothermalSIDM input
