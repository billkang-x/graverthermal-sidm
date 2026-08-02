"""Quick test: compute projected mass ratio for NFW with various R_MAX."""
import numpy as np
trapz = np.trapezoid

# NFW profile parameters
r_s = 0.085  # kpc
rho_0 = 10.0  # Msun/pc^3

def rho_nfw(r_pc):
    x = r_pc / (r_s * 1000.0)
    return rho_0 / (x * (1+x)**2)

def projected_enclosed_mass_code(r_2d_pc, r_max_pc=200000.0, n=10001):
    """Use the formula from rescale.py projected_enclosed_mass.

    M_2D(r_2d) = 2*pi * integral_0^{r_2d} r*rho dr
               + r_2d * integral_{r_2d}^{inf} 2r*rho/sqrt(r^2-r_2d^2) dr
    """
    # Part 1
    r_in = np.linspace(0.5, r_2d_pc, n // 2)
    rho_in = rho_nfw(r_in)
    M_in = 2 * np.pi * trapz(r_in * rho_in, r_in)
    # Part 2
    r_out = np.linspace(r_2d_pc * (1+1e-6), r_max_pc, n)
    rho_out = rho_nfw(r_out)
    integrand = 2 * r_out * rho_out / np.sqrt(r_out**2 - r_2d_pc**2)
    M_out = r_2d_pc * trapz(integrand, r_out)
    return M_in + M_out

print("=" * 60)
print("NFW projected enclosed mass (rescale.py formula)")
print("=" * 60)
for r_max_kpc in [0.34, 1.0, 5.0, 50.0, 200.0]:
    r_max_pc = r_max_kpc * 1000
    M_20 = projected_enclosed_mass_code(20.0, r_max_pc)
    M_90 = projected_enclosed_mass_code(90.0, r_max_pc)
    print(f"  R_MAX={r_max_kpc:6.2f} kpc: M_20={M_20:.4e}, M_90={M_90:.4e}, "
          f"ratio={M_20/M_90:.6f}")
print(f"\nExpected fluid_ratio_init = 0.319457")
