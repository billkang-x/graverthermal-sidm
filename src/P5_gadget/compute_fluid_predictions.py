"""
Compute fluid model predictions for the 5 N-body verification points.

Uses the same DissipativeHalo API as run_grid_scan.py.

Output: data/P5_nbody_verify/fluid_predictions.csv
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
from rescale import (
    projected_enclosed_mass,
    M_20PC, M_90PC, MASS_RATIO_OBS, MASS_RATIO_ERR,
)

# Test points - matching setup_nbody_verification.py
TEST_POINTS = [
    {"name": "P1_elastic_control", "model_key": "elastic",
     "sigma_m_100": 0.1, "r_diss_const": 1.0,
     "r_s_kpc": 0.085, "rho_0_msun_pc3": 10.0, "t_evo_gyr": 1.0},
    {"name": "P2_m3_low_sigma", "model_key": "const_rdiss_1p05",
     "sigma_m_100": 0.005, "r_diss_const": 1.05,
     "r_s_kpc": 0.085, "rho_0_msun_pc3": 10.0, "t_evo_gyr": 0.07},
    {"name": "P3_m3_high_sigma", "model_key": "const_rdiss_1p05",
     "sigma_m_100": 0.220, "r_diss_const": 1.05,
     "r_s_kpc": 0.085, "rho_0_msun_pc3": 10.0, "t_evo_gyr": 0.10},
    {"name": "P4_m1_low_sigma", "model_key": "M1_dark_photon_massive",
     "sigma_m_100": 0.005, "r_diss_const": None,
     "r_s_kpc": 0.085, "rho_0_msun_pc3": 10.0, "t_evo_gyr": 0.05},
    {"name": "P5_m1_high_sigma", "model_key": "M1_dark_photon_massive",
     "sigma_m_100": 0.165, "r_diss_const": None,
     "r_s_kpc": 0.085, "rho_0_msun_pc3": 10.0, "t_evo_gyr": 0.05},
]

N_SHELLS = 100
T_EPSILON = 1e-2
R_EPSILON = 1e-12
RHO_FACTOR_END = 1000.0
W_UNITS = 100.0
MAX_STEPS = 30000


def calibrate_alpha_D(p_model, target_sigma, v_ref=W_UNITS):
    """Rescale alpha_D so sigma_T_born(v_ref) = target_sigma."""
    sig_now = sigma_T_born(np.array([v_ref]), p_model)[0]
    if sig_now <= 0:
        return p_model
    ratio = target_sigma / sig_now
    import dataclasses
    return dataclasses.replace(p_model, alpha_D=p_model.alpha_D * np.sqrt(ratio))


def run_one_point(point):
    """Run fluid model for one test point, return final mass ratio + initial ratio."""
    name = point["name"]
    model_key = point["model_key"]
    r_s_kpc = point["r_s_kpc"]
    rho_0 = point["rho_0_msun_pc3"]
    sigma_m = point["sigma_m_100"]
    t_target_gyr = point["t_evo_gyr"]
    print(f"\n=== {name} (model={model_key}, sigma/m={sigma_m}) ===")

    # Create temp directory for HaloRecord
    out_dir = os.path.join(_PROJ_ROOT, "data", "P5_nbody_verify", "_tmp", name)
    os.makedirs(out_dir, exist_ok=True)
    # Clean any previous pickle files
    for f in os.listdir(out_dir):
        if f.endswith('.pkl'):
            os.remove(os.path.join(out_dir, f))

    rec = HaloRecord(out_dir)

    try:
        if model_key in ('M1_dark_photon_massive', 'M2_scalar_phi_massive'):
            bm = benchmark_models()
            p_model = bm[model_key]
            # Calibrate alpha_D to give sigma_m at v=100 km/s
            p_model = calibrate_alpha_D(p_model, sigma_m, v_ref=W_UNITS)
            # quad() requires scalar return; wrap to extract single value
            def _sig_scalar(v):
                return float(np.atleast_1d(sigma_T_born(np.atleast_1d(v), p_model))[0])
            def _rd_scalar(v):
                return float(np.atleast_1d(r_diss(np.atleast_1d(v), p_model))[0])
            sigma_m_eff, rdiss_eff, _, _ = effective_sigma_m_and_rdiss(
                _sig_scalar,
                _rd_scalar,
                np.logspace(1.5, 6.5, 40),
            )
            evo = DissipativeHalo(rec,
                                   sigma_m_eff_callable=sigma_m_eff,
                                   rdiss_eff_callable=rdiss_eff,
                                   flag_dissipation=True,
                                   profile='NFW', r_s=r_s_kpc, rho_s=rho_0,
                                   sigma_m_with_units=sigma_m, w_units=W_UNITS,
                                   n_shells=N_SHELLS, r_max=50.0, r_min=0.02,
                                   flag_hydrostatic_initial=True,
                                   flag_timestep_use_relaxation=True,
                                   flag_timestep_use_energy=True)
        elif model_key == 'const_rdiss_1p05':
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
                                   flag_dissipation=False,  # no dissipation
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

        # Initial state
        if evo.t == 0:
            if evo.flag_hydrostatic_initial:
                evo.hydrostatic_adjustment()
            evo.save_halo()

        # Compute initial ratio
        list_files, _ = rec.glob_pickle_files()
        data_init = rec.get_halo_state_pickled(file_halo=list_files[-1])
        scale_r_kpc = evo.scale_r.to('kpc').value
        scale_rho = evo.scale_rho.to('Msun/pc**3').value
        r_kpc_init = data_init['r'] * scale_r_kpc
        rho_arr_init = data_init['rho'] * scale_rho
        M_in_init = projected_enclosed_mass(r_kpc_init, rho_arr_init, 0.02,
                                             r_unit='kpc', rho_unit='Msun_pc3')
        M_out_init = projected_enclosed_mass(r_kpc_init, rho_arr_init, 0.09,
                                              r_unit='kpc', rho_unit='Msun_pc3')
        ratio_init = M_in_init / M_out_init if M_out_init > 0 else np.nan
        print(f"  Initial ratio: {ratio_init:.6f}  (M_in={M_in_init:.3e}, M_out={M_out_init:.3e})")

        # Evolve to t_target (in Gyr). Convert to dimensionless time.
        scale_t_gyr = evo.scale_t.to('Gyr').value
        t_target_dimless = t_target_gyr / scale_t_gyr
        print(f"  scale_t = {scale_t_gyr:.4f} Gyr, target t_dimless = {t_target_dimless:.4f}")

        # Evolve
        t_start_wall = time.time()
        n_steps_actual = 0
        for step in range(MAX_STEPS):
            if time.time() - t_start_wall > 120:  # 2 min wall budget
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
            n_steps_actual = step + 1

        # Final ratio
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

        print(f"  Evolved {n_steps_actual} steps to t={evo.t:.4f} (dimless) = {t_actual_gyr:.4f} Gyr")
        print(f"  Final ratio: {ratio_final:.6f}  (M_in={M_in_final:.3e}, M_out={M_out_final:.3e})")

        return {
            "name": name,
            "model_key": model_key,
            "sigma_m_100": sigma_m,
            "r_diss": point["r_diss_const"],
            "r_s_kpc": r_s_kpc,
            "rho_0_msun_pc3": rho_0,
            "t_target_gyr": t_target_gyr,
            "t_actual_gyr": t_actual_gyr,
            "n_steps": n_steps_actual,
            "fluid_ratio_init": ratio_init,
            "fluid_ratio_final": ratio_final,
            "fluid_m_inner_msun": M_in_final,
            "fluid_m_outer_msun": M_out_final,
        }

    except Exception as e:
        import traceback
        print(f"  ERROR: {e}")
        traceback.print_exc()
        return {
            "name": name,
            "model_key": model_key,
            "sigma_m_100": sigma_m,
            "r_diss": point["r_diss_const"],
            "r_s_kpc": r_s_kpc,
            "rho_0_msun_pc3": rho_0,
            "t_target_gyr": t_target_gyr,
            "error": str(e),
        }


def main():
    print("=" * 70)
    print("Fluid model predictions for 5 N-body verification points")
    print("=" * 70)

    results = []
    for point in TEST_POINTS:
        r = run_one_point(point)
        results.append(r)

    df = pd.DataFrame(results)
    out_csv = os.path.join(_PROJ_ROOT, "data", "P5_nbody_verify", "fluid_predictions.csv")
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"\nSaved {len(df)} predictions to {out_csv}")
    print()
    print(df[["name", "model_key", "sigma_m_100", "fluid_ratio_init", "fluid_ratio_final"]].to_string(index=False))


if __name__ == "__main__":
    main()
