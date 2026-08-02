#!/usr/bin/env python3
"""Generate a Gadget4 binary IC file (format 1) from the existing HDF5 IC.

This is more reliable than HDF5 format since the working reference job uses
binary ICs. The binary format has a 256-byte header followed by particle
data blocks.

Gadget2 binary format (NTYPES=3, GADGET2_HEADER defined):
  Header (256 bytes):
    npart[6] (uint32) — particles per type (we use 3 elements since NTYPES=3)
    mass[6] (double)
    time (double)
    redshift (double)
    flag_sfr (int)
    flag_feedback (int)
    npartTotal[6] (uint32)
    flag_cooling (int)
    num_files (int)
    BoxSize (double)
    Omega0 (double)
    OmegaLambda (double)
    HubbleParam (double)
    flag_stellarage (int)
    flag_metals (int)
    npartTotalHighWord[6] (uint32)
    flag_entropy_instead_u (int)
    ... (padding to 256 bytes)
  Then blocks (only for non-zero types):
    Coordinates (N, 3) float
    Velocities (N, 3) float
    ParticleIDs (N,) uint32/int64
    Masses (N,) float (only if mass=0 in header, otherwise use header mass)

Since our header sets MassTable[1]=624.67, we don't need to write Masses block.
Actually, Gadget2 format writes Masses for ALL types with npart>0 always
when using format 1 (legacy). Let me write all blocks to be safe.
"""
import h5py
import numpy as np
import struct
import os

INPUT_HDF5 = "D:/graverthermal-sidm/data/P5_nbody_verify/ics/ic_fixed_ntypes3.hdf5"
OUTPUT_BIN = "D:/graverthermal-sidm/data/P5_nbody_verify/ics/ic.dat"


def make_gadget2_header(npart, mass, time, redshift, boxsize,
                         omega0, omega_l, hubble, num_files=1):
    """Build a 256-byte Gadget2 IC header as bytes (NTYPES=6 array size for compat)."""
    # Pad arrays to 6 elements for header compatibility
    npart6 = np.zeros(6, dtype=np.uint32)
    npart6[:3] = npart[:3]
    mass6 = np.zeros(6, dtype=np.float64)
    mass6[:3] = mass[:3]

    flag_sfr = 0
    flag_feedback = 0
    npartTotal6 = npart6.copy()
    flag_cooling = 0
    flag_stellarage = 0
    flag_metals = 0
    npartTotalHW6 = np.zeros(6, dtype=np.uint32)
    flag_entropy = 0

    # Format: header struct as defined in Gadget2
    # See: src/data/allvars.h struct header
    fields = [
        ('6I', npart6),         # npart[6]
        ('6d', mass6),          # mass[6]
        ('d', time),
        ('d', redshift),
        ('i', flag_sfr),
        ('i', flag_feedback),
        ('6I', npartTotal6),
        ('i', flag_cooling),
        ('i', num_files),
        ('d', boxsize),
        ('d', omega0),
        ('d', omega_l),
        ('d', hubble),
        ('i', flag_stellarage),
        ('i', flag_metals),
        ('6I', npartTotalHW6),
        ('i', flag_entropy),
    ]

    # Build the buffer
    buf = b''
    for fmt, data in fields:
        buf += struct.pack(fmt, *np.asarray(data).ravel())

    # Pad to 256 bytes
    if len(buf) > 256:
        raise RuntimeError(f"Header too big: {len(buf)} bytes")
    buf += b'\x00' * (256 - len(buf))

    return buf


def main():
    print(f"Reading: {INPUT_HDF5}")
    with h5py.File(INPUT_HDF5, 'r') as f:
        coords = f['PartType1/Coordinates'][:].astype(np.float32)
        vels = f['PartType1/Velocities'][:].astype(np.float32)
        pids = f['PartType1/ParticleIDs'][:].astype(np.uint32)
        masses = f['PartType1/Masses'][:].astype(np.float32)
        header = dict(f['Header'].attrs)

    n = coords.shape[0]
    print(f"  N = {n} particles")
    print(f"  Coord range: [{coords.min():.4f}, {coords.max():.4f}]")
    print(f"  Vel range: [{vels.min():.4f}, {vels.max():.4f}]")
    print(f"  Mass: {masses[0]}")

    # Header arrays for NTYPES=3
    npart = np.array([0, n, 0], dtype=np.uint32)
    mass_table = np.array([0.0, float(masses[0]), 0.0], dtype=np.float64)

    print(f"\nWriting: {OUTPUT_BIN}")
    with open(OUTPUT_BIN, 'wb') as f:
        # Build the header bytes
        hdr_buf = make_gadget2_header(
            npart, mass_table,
            time=0.0, redshift=0.0,
            boxsize=float(header.get('BoxSize', 0.68)),
            omega0=float(header.get('Omega0', 1.0)),
            omega_l=float(header.get('OmegaLambda', 0.0)),
            hubble=float(header.get('HubbleParam', 0.7)),
            num_files=1,
        )
        # Write header block: blksize, data, blksize
        f.write(struct.pack('I', 256))
        f.write(hdr_buf)
        f.write(struct.pack('I', 256))
        print(f"  Header block: 256 + 8 = 264 bytes")

        # Helper to write a block with size markers
        def write_block(name, data_bytes):
            n = len(data_bytes)
            f.write(struct.pack('I', n))
            f.write(data_bytes)
            f.write(struct.pack('I', n))
            print(f"  {name} block: {n} + 8 = {n+8} bytes")

        # Block 2: Coordinates (N, 3) float32
        write_block("Coordinates", coords.tobytes())

        # Block 3: Velocities (N, 3) float32
        write_block("Velocities", vels.tobytes())

        # Block 4: ParticleIDs (N,) uint32
        write_block("ParticleIDs", pids.tobytes())

        # Note: Masses block omitted since header sets mass[1] != 0
        # (Gadget2 uses header mass when it's nonzero)

    print(f"\nTotal file size: {os.path.getsize(OUTPUT_BIN)} bytes")
    print("Done.")


if __name__ == "__main__":
    main()
