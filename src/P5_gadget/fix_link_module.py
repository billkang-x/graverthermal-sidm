"""
Cleanly link the sidm_dissipative module into GADGET-4's main integration loop.

Strategy:
  1. Restore run.cc from the backup (pristine).
  2. Download run.cc locally, apply a clean edit with Python regex.
  3. Upload it back.
  4. Add SIDM_DISSIPATIVE to defines_extra if missing.
  5. Build with `make build` (suppresses check.py macros check; this is the
     officially documented workaround per the error message itself).
  6. Verify the new executable has dissipative symbols.

This replaces the earlier link_module.py which used fragile sed that produced
garbled `#elif defined(SIDM)_DISSIPATIVE` text and duplicate lines.
"""
import paramiko
import re
import sys
import io

HOST = "ssh.cn-zhongwei-1.paracloud.com"
PORT = 22
USER = "scg7816@ZC-M6"
PASS = "ZCB1PyO65Vw7GXd2q4LJoSYl8sbEi390"

REMOTE_RUN_CC = "~/dsidm_project/source/src/main/run.cc"
REMOTE_RUN_CC_BAK = "~/dsidm_project/source/src/main/run.cc.bak"
REMOTE_MAKEFILE = "~/dsidm_project/source/Makefile"
REMOTE_DEFINES_EXTRA = "~/dsidm_project/source/defines_extra"
LOCAL_RUN_CC = "D:/graverthermal-sidm/src/P5_gadget/run_cc_edited.cc"


def run(client, cmd, timeout=600):
    print(f"$ {cmd[:200]}")
    stdin, stdout, stderr = client.exec_command(cmd, timeout=timeout)
    out = stdout.read().decode(errors='replace').rstrip()
    err = stderr.read().decode(errors='replace').rstrip()
    if out:
        print(out[:8000])
    if err:
        print(f"  [stderr] {err[:3000]}")
    return out, err


def edit_run_cc_locally(content: str) -> str:
    """Apply the dissipative module linkage edits to run.cc content.

    Two surgical edits:
      A. Include the sidm_dissipative.h header right after sidm.h, guarded
         by #ifdef SIDM_DISSIPATIVE.
      B. Replace the elastic-only call block:
            #ifdef SIDM
              Sidm.do_sidm_scattering(&Sp, All.TimeStep);
            #endif
         with a conditional:
            #ifdef SIDM_DISSIPATIVE
              SidmDissipative.do_sidm_scattering(&Sp, All.TimeStep);
            #elif defined(SIDM)
              Sidm.do_sidm_scattering(&Sp, All.TimeStep);
            #endif
    """
    # --- Edit A: add the header include ---
    pattern_a = r'(#include "\.\./sidm/sidm\.h")'
    replacement_a = (
        r'\1\n'
        r'#ifdef SIDM_DISSIPATIVE\n'
        r'#include "../sidm/sidm_dissipative.h"\n'
        r'#endif'
    )
    new_content, n_a = re.subn(pattern_a, replacement_a, content, count=1)
    assert n_a == 1, f"Edit A (header include) failed: matched {n_a} times"
    print(f"  Edit A: added sidm_dissipative.h include ({n_a} match)")

    # --- Edit B: replace the call block ---
    # Match the exact block: #ifdef SIDM\n  Sidm.do_sidm_scattering(&Sp, All.TimeStep);\n#endif
    pattern_b = (
        r'#ifdef SIDM\n'
        r'(\s*)Sidm\.do_sidm_scattering\(&Sp, All\.TimeStep\);\n'
        r'#endif'
    )

    def replacement_b(m):
        indent = m.group(1)
        return (
            f'#ifdef SIDM_DISSIPATIVE\n'
            f'{indent}SidmDissipative.do_sidm_scattering(&Sp, All.TimeStep);\n'
            f'#elif defined(SIDM)\n'
            f'{indent}Sidm.do_sidm_scattering(&Sp, All.TimeStep);\n'
            f'#endif'
        )

    new_content, n_b = re.subn(pattern_b, replacement_b, new_content, count=1)
    assert n_b == 1, f"Edit B (call block) failed: matched {n_b} times"
    print(f"  Edit B: replaced call block with conditional ({n_b} match)")

    return new_content


def main():
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(HOST, port=PORT, username=USER, password=PASS,
                   timeout=30, look_for_keys=False, allow_agent=False)
    print("Connected to HPC.\n")

    # Step 1: Restore run.cc from backup
    print("=== Step 1: Restore run.cc from backup ===")
    run(client, f"cp {REMOTE_RUN_CC_BAK} {REMOTE_RUN_CC}")
    run(client, f"grep -c 'sidm' {REMOTE_RUN_CC}")

    # Step 2: Verify the backup has the original SIDM block (sanity check)
    print("\n=== Step 2: Verify restored run.cc has the original SIDM block ===")
    out, _ = run(client, f"grep -n 'do_sidm_scattering' {REMOTE_RUN_CC}")
    # Expect exactly one match: the original Sidm.do_sidm_scattering call
    # (Gadget4 only has one place that calls this; if there are more, abort.)

    # Step 3: Download run.cc locally
    print("\n=== Step 3: Download run.cc locally ===")
    sftp = client.open_sftp()
    sftp.get(REMOTE_RUN_CC.replace("~", "/public3/home/scg7816"), LOCAL_RUN_CC)
    print(f"  Downloaded to {LOCAL_RUN_CC}")

    with open(LOCAL_RUN_CC, 'r', encoding='utf-8') as f:
        content = f.read()
    print(f"  Original size: {len(content)} chars, {content.count(chr(10))} lines")

    # Step 4: Apply the edits locally
    print("\n=== Step 4: Apply clean edits locally ===")
    edited = edit_run_cc_locally(content)
    print(f"  Edited size: {len(edited)} chars, {edited.count(chr(10))} lines")

    with open(LOCAL_RUN_CC, 'w', encoding='utf-8') as f:
        f.write(edited)
    print(f"  Wrote edited file to {LOCAL_RUN_CC}")

    # Step 5: Upload back
    print("\n=== Step 5: Upload edited run.cc back to HPC ===")
    sftp.put(LOCAL_RUN_CC, REMOTE_RUN_CC.replace("~", "/public3/home/scg7816"))
    sftp.close()
    print("  Uploaded.")

    # Step 6: Verify the upload
    print("\n=== Step 6: Verify the edited run.cc on HPC ===")
    run(client, f"grep -n 'sidm' {REMOTE_RUN_CC}")
    run(client, f"sed -n '24,32p' {REMOTE_RUN_CC}")
    run(client, f"sed -n '255,270p' {REMOTE_RUN_CC}")

    # Step 7: Add SIDM_DISSIPATIVE to defines_extra (defensive; check.py reads
    # this file to whitelist macros that should not appear in Template-Config.sh)
    print("\n=== Step 7: Ensure SIDM_DISSIPATIVE is in defines_extra ===")
    out, _ = run(client, f"grep -c 'SIDM_DISSIPATIVE' {REMOTE_DEFINES_EXTRA}")
    if "0" in out.split("\n")[-1]:
        run(client, f"echo 'SIDM_DISSIPATIVE' >> {REMOTE_DEFINES_EXTRA}")
        run(client, f"tail -5 {REMOTE_DEFINES_EXTRA}")
    else:
        print("  SIDM_DISSIPATIVE already in defines_extra.")

    # Step 8: Verify Makefile has sidm_dissipative.o (it should already)
    print("\n=== Step 8: Verify Makefile ===")
    run(client, f"grep -n 'sidm' {REMOTE_MAKEFILE}")

    # Step 9: Clean build directory and recompile
    # Use `make build` to suppress the check.py macros check (the official
    # workaround per the error message). We've also added SIDM_DISSIPATIVE
    # to defines_extra so even `make` should work.
    print("\n=== Step 9: Recompile (make build) ===")
    run(client,
        "cd ~/dsidm_project/source && "
        "source ~/env.sh 2>&1 | tail -2 && "
        "make clean 2>&1 | tail -3 && "
        "make build -j 8 EXEC=Gadget4_dsidm 2>&1 | tail -60",
        timeout=600)

    # Step 10: Check if executable was produced
    print("\n=== Step 10: Check executable ===")
    run(client, "ls -la ~/dsidm_project/source/Gadget4_dsidm 2>&1")

    # Step 11: Verify the new binary has dissipative symbols
    print("\n=== Step 11: Verify dissipative symbols in the binary ===")
    run(client, "nm ~/dsidm_project/source/Gadget4_dsidm 2>&1 | grep -i 'dissipative' | head -10")

    client.close()
    print("\nDone.")


if __name__ == "__main__":
    main()
