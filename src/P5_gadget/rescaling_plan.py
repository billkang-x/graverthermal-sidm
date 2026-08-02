#!/usr/bin/env python3
"""
方案1: 重标度对称性 N-body 验证

核心思想:
  我们的物理空间参数 (r_s=0.085 kpc, ρ_0=10 Msun/pc³) 太致密，
  直接 N-body 模拟时间步太小 (dt~1e-8)，无法在合理时间内完成。

  利用 Schmidt et al. 2026 附录 G 的重标度对称性:
    r → λ r,  m → μ m,  v → sqrt(μ/λ) v,  t → sqrt(λ³/μ) t
    ρ → (μ/λ³) ρ,  σ_T/m → (λ²/μ) σ_T/m

  将物理空间映射到模拟空间 (r_s=3.6 kpc, ρ_0=7.09e-3 Msun/pc³)，
  在模拟空间运行 N-body（动力学时间更长，时间步更大），
  然后将结果映射回来。

重标度参数:
  λ = r_s_phys / r_s_sim = 0.085 / 3.6 = 0.023611
  μ = (ρ_0_phys / ρ_0_sim) * λ³ = (10 / 0.00709) * 1.316e-5 = 0.018565

  速度: v_phys = 0.8867 * v_sim
  时间: t_phys = 0.02663 * t_sim  (即 t_sim = 37.55 * t_phys)
  截面: σ_phys = 0.03003 * σ_sim  (即 σ_sim = 33.31 * σ_phys)

模拟空间的动力学时间: 0.0498 Gyr (是物理空间的 37.5 倍)
"""
import numpy as np

# === 重标度参数 ===
R_S_PHYS = 0.085    # kpc
RHO_0_PHYS = 10.0   # Msun/pc^3
R_S_SIM = 3.6       # kpc
RHO_0_SIM = 7.09e-3 # Msun/pc^3

LAMBDA = R_S_PHYS / R_S_SIM
MU = (RHO_0_PHYS / RHO_0_SIM) * LAMBDA**3
V_SCALE = np.sqrt(MU / LAMBDA)        # v_phys = V_SCALE * v_sim
T_SCALE = np.sqrt(LAMBDA**3 / MU)     # t_phys = T_SCALE * t_sim
SIGMA_SCALE = LAMBDA**2 / MU          # σ_phys = SIGMA_SCALE * σ_sim

print("=" * 70)
print("重标度对称性 N-body 验证方案")
print("=" * 70)
print(f"\n物理空间: r_s={R_S_PHYS} kpc, ρ_0={RHO_0_PHYS} Msun/pc³")
print(f"模拟空间: r_s={R_S_SIM} kpc, ρ_0={RHO_0_SIM} Msun/pc³")
print(f"\nλ = {LAMBDA:.6f}")
print(f"μ = {MU:.6f}")
print(f"v_scale = {V_SCALE:.4f} (v_phys = v_scale * v_sim)")
print(f"t_scale = {T_SCALE:.6f} (t_phys = t_scale * t_sim)")
print(f"σ_scale = {SIGMA_SCALE:.6f} (σ_phys = σ_scale * σ_sim)")

# === 3 个验证点的参数转换 ===
POINTS = [
    ("P1_elastic_control", 0.1,   1.0,  0.68),
    ("P2_m3_low_sigma",    0.005, 1.05, 0.07),
    ("P3_m3_high_sigma",   0.22,  1.05, 0.10),
]

print("\n" + "=" * 70)
print("3 个验证点的参数转换 (物理空间 → 模拟空间)")
print("=" * 70)
print(f"\n{'点名':<25} {'σ_phys':>8} {'σ_sim':>10} {'r_diss':>8} {'t_phys(Gyr)':>12} {'t_sim(Gyr)':>12} {'t_sim(code)':>12}")
print("-" * 90)

sim_points = []
for name, sigma_phys, r_diss, t_phys in POINTS:
    sigma_sim = sigma_phys / SIGMA_SCALE  # σ_sim = σ_phys / σ_scale
    t_sim = t_phys / T_SCALE              # t_sim = t_phys / t_scale
    # In code units: t_unit = 0.978 Gyr
    t_sim_code = t_sim / 0.978
    print(f"{name:<25} {sigma_phys:>8.3f} {sigma_sim:>10.4f} {r_diss:>8.2f} {t_phys:>12.4f} {t_sim:>12.4f} {t_sim_code:>12.4f}")
    sim_points.append((name, sigma_sim, r_diss, t_sim, t_sim_code))

print("\n" + "=" * 70)
print("关键优势分析")
print("=" * 70)

# 物理空间动力学时间
G_cgs = 6.674e-8
rho_phys_cgs = RHO_0_PHYS * 2e33 / (3.086e18)**3
t_dyn_phys = 1.0 / np.sqrt(4 * np.pi * rho_phys_cgs * G_cgs)
rho_sim_cgs = RHO_0_SIM * 2e33 / (3.086e18)**3
t_dyn_sim = 1.0 / np.sqrt(4 * np.pi * rho_sim_cgs * G_cgs)

print(f"\n物理空间动力学时间: {t_dyn_phys/3.156e16:.6f} Gyr")
print(f"模拟空间动力学时间: {t_dyn_sim/3.156e16:.6f} Gyr")
print(f"比值: {t_dyn_sim/t_dyn_phys:.1f}x (模拟空间更慢)")

print(f"\n物理空间 N_body 步数估算 (dt~1e-8, t_target~0.07-0.68):")
for name, _, _, t_phys, _ in POINTS:
    n_steps = t_phys / (1e-8 * 0.978)  # code units
    print(f"  {name}: ~{n_steps:.1e} 步 (不可行)")

print(f"\n模拟空间 N_body 步数估算:")
# 模拟空间的特征速度
v_circ_sim = np.sqrt(4.302e-3 * 4 * np.pi * RHO_0_SIM * 1e9 * R_S_SIM**3 *
                      (np.log1p(4) - 4/5) / (4 * R_S_SIM))
# 特征加速度 a ~ v_circ^2 / r_s
a_char = v_circ_sim**2 / R_S_SIM
# 时间步 dt ~ sqrt(epsilon / a) ~ 0.1 * sqrt(r_s / a) ~ r_s / v_circ
dt_est = 0.01 * R_S_SIM / v_circ_sim  # 1% of dynamical time
for name, sigma_sim, _, t_sim, t_sim_code in sim_points:
    n_steps = t_sim_code / dt_est
    print(f"  {name}: σ_sim={sigma_sim:.2f}, t_sim={t_sim:.3f} Gyr, "
          f"dt~{dt_est:.4f}, ~{n_steps:.0f} 步")

print(f"\n模拟空间 v_circ at r=4*r_s: {v_circ_sim:.1f} km/s")
print(f"物理空间 v_circ at r=4*r_s: {v_circ_sim * V_SCALE:.1f} km/s")
