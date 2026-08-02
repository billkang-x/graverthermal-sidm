"""
Re-run fluid model at t_phys corresponding to our N-body t_code=1.0.

Our N-body runs go to t_code=1.0 (Gadget4 code units).
With UnitTime=0.978 Gyr and T_SCALE=0.02663 (rescaling):
  t_phys = t_code * 0.978 * T_SCALE = 0.02604 Gyr

This is much shorter than the original target times (0.07-1.0 Gyr).
We need fluid model predictions at this earlier time for proper comparison.

Also computes intermediate checkpoints at:
  t_code=0.5 -> t_phys = 0.01302 Gyr
  t_code=1.0 -> t_phys = 0.02604 Gyr

Output: data/P5_nbody_verify/fluid_predictions_nbody_time.csv
"""
from __future__ import annotations
import os, sys, time
import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ_ROOT = os.path.normpath(os.path.join(_HERE, '..', '..'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'cross_sections'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'fluid_runner'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'external', 'gravothermalsidm'))
sys.path.insert(0, os.path.join(_PROJ_ROOT, 'src', 'P3_rescaling'))

import matplotlib
matplotlib.use('Agg')

from astropy import units as ut
from astropy import constants as ct

from dsidm_models import benchmark_models, sigma_T_born, r_diss
from thermal_avg import effective_sigma_m_and_rdiss
from dissipative_halo import DissipativeHalo
from SourcePy.evolve import Halo
from SourcePy.record import HaloRecord
from rescale import projected_enclosed_mass

# Rescaling parameters
LAMBDA = 0.085 / 3.6
MU = (10.0 / 7.09e-3) * LAMBDA**3
T_SCALE = np.sqrt(LAMBDA**3 / MU)  # t_phys = T_SCALE * t_sim

# N-body time checkpoints (t_code -> t_phys in Gyr)
# t_sim = t_code * 0.978 Gyr (UnitTime)
# t_phys = t_sim * T_SCALE
T_CODE_CHECKPOINTS = [0.0, 0.5, 1.0]
T_PHYS_CHECKPOINTS = [tc * 0.978 * T_SCALE for tc in T_CODE_CHECKPOINTS]

# Test points - same as compute_fluid_predictions.py
TEST_POINTS = [
    {"name": "P1_elastic_control", "model_key": "elastic",
     "sigma_m_100": 0.1, "r_diss_const": 1.0,
     "r_s_kpc": 0.085, "rho_0_msun_pc3": 10.0},
    {"name": "P2_m3_low_sigma", "model_key": "const_rdiss_1p05",
     "sigma_m_100": 0.005, "r_diss_const": 1.05,
     "r_s_kpc": 0.085, "rho_0_msun_pc3": 10.0},
    {"name": "P3_m3_high_sigma", "model_key": "const_rdiss_1p05",
     "sigma_m_100": 0.220, "r_diss_const": 1.05,
     "r_s_kpc": 0.085, "rho_0_msun_pc3": 10.0},
]

N_SHELLS = 100
T_EPSILON = 1e-2
R_EPSILON = 1e-12
RHO_FACTOR_END = 1000.0
W_UNITS = 100.0
MAX_STEPS = 30000


def calibrate_alpha_D(p_model, target_sigma, v_ref=W_UNITS):
    sig_now = sigma_T_born(np.array([v_ref]), p_model)[0]
    if sig_now <= 0:
        return p_model
    ratio = target_sigma / sig_now
    import dataclasses
    return dataclasses.replace(p_model, alpha_D=p_model.alpha_D * np.sqrt(ratio))


def run_one_point_with_checkpoints(point, t_phys_checkpoints):
    """Run fluid model, recording mass ratio at specified physical times."""
    name = point["name"]
    model_key = point["model_key"]
    r_s_kpc = point["r_s_kpc"]
    rho_0 = point["rho_0_msun_pc3"]
    sigma_m = point["sigma_m_100"]

    print(f"\n=== {name} (model={model_key}, sigma/m={sigma_m}) ===")

    out_dir = os.path.join(_PROJ_ROOT, "data", "P5_nbody_verify", "_tmp_nbody_time", name)
    os.makedirs(out_dir, exist_ok=True)
    for f in os.listdir(out_dir):
        if f.endswith('.pkl'):
            os.remove(os.path.join(out_dir, f))

    rec = HaloRecord(out_dir)

    try:
        if model_key == 'const_rdiss_1p05':
            sig_fn = lambda T: np.full(np.atleast_1d(T).shape, sigma_m, dtype=float)
            rd_fn = lambda T: np.full(np.atleast_1d(T).shape, 1.05, dtype=float)
            evo = DissipativeHalo(rec,
                                   sigma_m_eff_callable=sig_fn,
                                   rdiss_eff_callable=rd_fn,
                                   flag_dissipation=True,
                                   profile='NFW', r_s=r_s_kpc, rho_s=rho_0,
                                   sigma_m_with_units=sigma_m, w_units=W_UNITS,
                                   n_shells=N_SHELLS, r_max=50.0, r_min=0.02,
                                   flag_hydrostatic_initial=True,
                                   flag_timestep_use_relaxation=True,
                                   flag_timestep_use_energy=True)
        elif model_key == 'elastic':
            sig_fn = lambda T: np.full(np.atleast_1d(T).shape, sigma_m, dtype=float)
            rd_fn = lambda T: np.full(np.atleast_1d(T).shape, 1.0, dtype=float)
            evo = DissipativeHalo(rec,
                                   sigma_m_eff_callable=sig_fn,
                                   rdiss_eff_callable=rd_fn,
                                   flag_dissipation=False,
                                   profile='NFW', r_s=r_s_kpc, rho_s=rho_0,
                                   sigma_m_with_units=sigma_m, w_units=W_UNITS,
                                   n_shells=N_SHELLS, r_max=50.0, r_min=0.02,
                                   flag_hydrostatic_initial=True,
                                   flag_timestep_use_relaxation=True,
                                   flag_timestep_use_energy=True)
        else:
            raise ValueError(f"Unknown model_key: {model_key}")

        evo.t_epsilon = T_EPSILON
        evo.r_epsilon = R_EPSILON

        if evo.t == 0:
            if evo.flag_hydrostatic_initial:
                evo.hydrostatic_adjustment()
            evo.save_halo()

        # Compute scale factors
        scale_r_kpc = evo.scale_r.to('kpc').value
        scale_rho = evo.scale_rho.to('Msun/pc**3').value
        scale_t_gyr = evo.scale_t.to('Gyr').value

        # Convert physical time checkpoints to dimensionless
        t_dimless_checkpoints = [t / scale_t_gyr for t in t_phys_checkpoints]

        print(f"  scale_t = {scale_t_gyr:.6f} Gyr")
        print(f"  Checkpoints (Gyr -> dimless):")
        for t_gyr, t_dim in zip(t_phys_checkpoints, t_dimless_checkpoints):
            print(f"    t_phys={t_gyr:.6f} Gyr -> t_dimless={t_dim:.2f}")

        # Compute initial ratio
        list_files, _ = rec.glob_pickle_files()
        data_init = rec.get_halo_state_pickled(file_halo=list_files[-1])
        r_kpc_init = data_init['r'] * scale_r_kpc
        rho_arr_init = data_init['rho'] * scale_rho
        M_in_init = projected_enclosed_mass(r_kpc_init, rho_arr_init, 0.02,
                                             r_unit='kpc', rho_unit='Msun_pc3')
        M_out_init = projected_enclosed_mass(r_kpc_init, rho_arr_init, 0.09,
                                              r_unit='kpc', rho_unit='Msun_pc3')
        ratio_init = M_in_init / M_out_init if M_out_init > 0 else np.nan
        print(f"  Initial ratio: {ratio_init:.6f}")

        # Evolve, recording at each checkpoint
        results = []
        checkpoint_idx = 1  # skip 0 (initial)
        t_target_dimless = t_dimless_checkpoints[-1]

        t_start_wall = time.time()
        for step in range(MAX_STEPS):
            if time.time() - t_start_wall > 300:  # 5 min budget
                print(f"  Wall budget hit at step {step}, t_dimless={evo.t:.4f}")
                break
            if evo.t >= t_target_dimless:
                break
            try:
                evo.conduct_heat()
                evo.hydrostatic_adjustment()
            except Exception as e:
                print(f"  Evolve failed at step {step}: {e}")
                break
            if np.any(np.isnan(evo.rho)) or np.any(np.isnan(evo.r)):
                print(f"  NaN at step {step}")
                break

            # Check if we've passed a checkpoint
            while (checkpoint_idx < len(t_dimless_checkpoints) and
                   evo.t >= t_dimless_checkpoints[checkpoint_idx]):
                evo.save_halo()
                list_files, _ = rec.glob_pickle_files()
                data = rec.get_halo_state_pickled(file_halo=list_files[-1])
                r_kpc = data['r'] * scale_r_kpc
                rho_arr = data['rho'] * scale_rho
                M_in = projected_enclosed_mass(r_kpc, rho_arr, 0.02,
                                                r_unit='kpc', rho_unit='Msun_pc3')
                M_out = projected_enclosed_mass(r_kpc, rho_arr, 0.09,
                                                 r_unit='kpc', rho_unit='Msun_pc3')
                ratio = M_in / M_out if M_out > 0 else np.nan
                t_gyr = evo.t * scale_t_gyr
                t_code = t_gyr / (0.978 * T_SCALE)

                print(f"  Checkpoint {checkpoint_idx}: t_dimless={evo.t:.2f}, "
                      f"t_phys={t_gyr:.6f} Gyr, ratio={ratio:.6f}")

                results.append({
                    "name": name,
                    "model_key": model_key,
                    "sigma_m_100": sigma_m,
                    "r_diss": point["r_diss_const"],
                    "r_s_kpc": r_s_kpc,
                    "rho_0_msun_pc3": rho_0,
                    "t_code": t_code,
                    "t_phys_gyr": t_gyr,
                    "fluid_ratio": ratio,
                    "fluid_ratio_init": ratio_init,
                    "fluid_m_inner_msun": M_in,
                    "fluid_m_outer_msun": M_out,
                })
                checkpoint_idx += 1

        # Final state
        evo.save_halo()
        list_files, _ = rec.glob_pickle_files()
        data_final = rec.get_halo_state_pickled(file_halo=list_files[-1])
        r_kpc_final = data_final['r'] * scale_r_kpc
        rho_arr_final = data_final['rho'] * scale_rho
        M_in_final = projected_enclosed_mass(r_kpc_final, rho_arr_final, 0.02,
                                              r_unit='kpc', rho_unit='Msun_pc3')
        M_out_final = projected_enclosed_mass(r_kpc_final, rho_arr_final, 0.09,
                                               r_unit='kpc', rho_unit='Msun_pc3')
        ratio_final = M_in_final / M_out_final if M_out_final > 0 else np.nan
        t_actual_gyr = evo.t * scale_t_gyr
        t_code_actual = t_actual_gyr / (0.978 * T_SCALE)

        print(f"  Final: t_phys={t_actual_gyr:.6f} Gyr (t_code={t_code_actual:.4f}), "
              f"ratio={ratio_final:.6f}")

        # Add final state if not already captured
        if not results or abs(results[-1]['t_phys_gyr'] - t_actual_gyr) > 0.001:
            results.append({
                "name": name,
                "model_key": model_key,
                "sigma_m_100": sigma_m,
                "r_diss": point["r_diss_const"],
                "r_s_kpc": r_s_kpc,
                "rho_0_msun_pc3": rho_0,
                "t_code": t_code_actual,
                "t_phys_gyr": t_actual_gyr,
                "fluid_ratio": ratio_final,
                "fluid_ratio_init": ratio_init,
                "fluid_m_inner_msun": M_in_final,
                "fluid_m_outer_msun": M_out_final,
            })

        # Add initial state
        results.insert(0, {
            "name": name,
            "model_key": model_key,
            "sigma_m_100": sigma_m,
            "r_diss": point["r_diss_const"],
            "r_s_kpc": r_s_kpc,
            "rho_0_msun_pc3": rho_0,
            "t_code": 0.0,
            "t_phys_gyr": 0.0,
            "fluid_ratio": ratio_init,
            "fluid_ratio_init": ratio_init,
            "fluid_m_inner_msun": M_in_init,
            "fluid_m_outer_msun": M_out_init,
        })

        return results

    except Exception as e:
        import traceback
        print(f"  ERROR: {e}")
        traceback.print_exc()
        return []


def main():
    print("=" * 70)
    print("Fluid model at N-body time checkpoints")
    print("=" * 70)
    print(f"T_SCALE = {T_SCALE:.6f}")
    print(f"N-body time checkpoints:")
    for tc, tp in zip(T_CODE_CHECKPOINTS, T_PHYS_CHECKPOINTS):
        print(f"  t_code={tc:.1f} -> t_phys={tp:.6f} Gyr = {tp*1000:.2f} Myr")

    all_results = []
    for point in TEST_POINTS:
        results = run_one_point_with_checkpoints(point, T_PHYS_CHECKPOINTS)
        all_results.extend(results)

    df = pd.DataFrame(all_results)
    out_csv = os.path.join(_PROJ_ROOT, "data", "P5_nbody_verify",
                           "fluid_predictions_nbody_time.csv")
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"\nSaved {len(df)} predictions to {out_csv}")
    print()
    print(df[["name", "t_code", "t_phys_gyr", "fluid_ratio", "fluid_ratio_init"]].to_string(index=False))


if __name__ == "__main__":
    main()
