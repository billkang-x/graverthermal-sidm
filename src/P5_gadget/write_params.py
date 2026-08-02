"""Write a complete, proper Gadget4 params.txt for each test point.

The previous params.txt was missing many required fields, causing Gadget4 to
fail at startup (MaxMemSize=0 -> mymalloc_fullinfo fails to allocate).
"""
import paramiko
import os

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

# Per-point parameters
POINTS = [
    ("P1_elastic_control", 0.1,   1.0,  "Elastic control: no dissipation"),
    ("P2_m3_low_sigma",    0.005, 1.05, "M3 low sigma/m: weak dissipative cooling"),
    ("P3_m3_high_sigma",   0.220, 1.05, "M3 high sigma/m: strong cooling"),
]


def gen_params_content(name, sigma_m, r_diss, description):
    """Generate a complete Gadget4 params.txt with all required fields."""
    return f"""%paramfile for N-body verification
% Point: {name}
% Description: {description}
% SIDM sigma/m = {sigma_m} cm^2/g (compile-time via Config.sh)
% SIDM r_diss = {r_diss} (compile-time via Config.sh)

InitCondFile      ic
OutputDir         output
SnapshotFileBase  snap
RestartDir        restart
OutputListFilename output_list.txt
OutputListOn       0

TimeMax           1.0
TimeBegin         0.0
TimeBetSnapshot   0.1
TimeBetStatistics 0.01
TimeOfFirstSnapshot 0.0
CpuTimeBetRestartFile 36000.0

% Timestep parameters
ErrTolIntAccuracy  0.025
ErrTolTheta         0.5
ErrTolThetaMax      0.5
ErrTolForce         0.005
ErrTolForceAcc      0.005
MaxRMSDisplacementFac 0.2
MaxSizeTimestep     0.05
MinSizeTimestep     0.0
TypeOfOpeningCriterion 1
CourantFac         0.15
DesNumNgb          32
MaxNumNgbDeviation 2
ArtBulkViscConst   1.0
TopNodeFactor      3.0
ActivePartFracForNewDomainDecomp 0.01
CpuTreeDomainUpdate 0.1
ComovingIntegrationOn 0
MaxMemSize         1000.0

% I/O format
ICFormat           3
SnapFormat         3
NumFilesPerSnapshot 1
MaxFilesWithConcurrentIO 1

% Softening (in code units = kpc; we want 0.0005 kpc = 0.5 pc)
SofteningClassOfPartType0  1
SofteningClassOfPartType1  1
SofteningClassOfPartType2  1
SofteningComovingClass0    0.0005
SofteningMaxPhysClass0     0.0005
SofteningComovingClass1    0.0005
SofteningMaxPhysClass1     0.0005
SofteningComovingClass2    0.0005
SofteningMaxPhysClass2     0.0005

% Units (kpc, Msun, km/s)
UnitLength_in_cm          3.085678e21
UnitMass_in_g             1.989e43
UnitVelocity_in_cm_per_s  1.0e5
UnitLuminosity_in_erg_s   1.0
UnitEnergy_in_ergs        1.0
GravityConstantInternal  0.017778279

% Cosmology (no cosmology - isolated halo)
MinGasTemp                100.0
HubbleParam               0.7
Omega0                    1.0
OmegaBaryon               0.0
OmegaLambda               0.0
Hubble                     0.0

% Box (non-periodic)
BoxSize                   0.5
PeriodicBoundariesOn      0

% Cooling (turned off via Config.sh but parameters required)
InitGasTemp              1000.0
MinEgySpec               0.05
TreecoolFile             treecool.txt
"""


def run(client, cmd, timeout=60):
    print(f"$ {cmd[:200]}", flush=True)
    stdin, stdout, stderr = client.exec_command(cmd, timeout=timeout)
    out = stdout.read().decode(errors='replace').rstrip()
    err = stderr.read().decode(errors='replace').rstrip()
    if out:
        print(out[:3000], flush=True)
    if err:
        print(f"  [stderr] {err[:2000]}", flush=True)
    return out, err


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected.\n", flush=True)

    sftp = client.open_sftp()
    base_remote = "/public3/home/scg7816/dsidm_project/nbody_verify"

    for name, sigma_m, r_diss, desc in POINTS:
        params_content = gen_params_content(name, sigma_m, r_diss, desc)
        local_path = f"D:/graverthermal-sidm/src/P5_gadget/ic_scripts/{name}_params_v2.txt"
        with open(local_path, 'w') as f:
            f.write(params_content)
        sftp.put(local_path, f"{base_remote}/{name}/params.txt")
        print(f"  {name}: params.txt uploaded ({len(params_content)} chars)")

    # Also create the output_list.txt (empty - we use TimeBetSnapshot instead)
    for name, _, _, _ in POINTS:
        try:
            with sftp.open(f"{base_remote}/{name}/output_list.txt", 'w') as f:
                f.write("")  # empty file
        except Exception as e:
            print(f"  {name}: warning - couldn't write output_list.txt: {e}")

    # Also create treecool.txt placeholder (required by CoolingFile but unused
    # since SIDM-only run has no gas)
    for name, _, _, _ in POINTS:
        try:
            with sftp.open(f"{base_remote}/{name}/treecool.txt", 'w') as f:
                f.write("# empty treecool - gas cooling not used\n")
        except Exception as e:
            print(f"  {name}: warning - couldn't write treecool.txt: {e}")

    sftp.close()

    # Verify
    print("\n=== Verify ===", flush=True)
    for name, _, _, _ in POINTS:
        run(client, f"head -10 ~/dsidm_project/nbody_verify/{name}/params.txt")

    client.close()
    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
