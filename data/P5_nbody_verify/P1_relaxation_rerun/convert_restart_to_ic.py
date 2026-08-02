#!/usr/bin/env python3
"""Convert a Gadget4 binary snapshot into a binary IC usable as InitCondFile.

Phase A of the P1 relaxation rerun writes a snapshot at t=0.5 code units. To
use that particle distribution as the t=0 initial condition for Phase B, we
rewrite it as a Gadget2-binary IC (ICFormat=1) with:

  - Time header field reset to 0.0
  - Redshift reset to 0.0
  - Particle IDs preserved (Gadget4 may have reordered them)
  - All other header fields copied through

The reader mirrors src/P5_gadget/read_binary_snap.py so the block layout is
guaranteed to match what the Gadget4 writer produced.

Usage:
    python convert_restart_to_ic.py <input_snapshot> <output_ic>

Default paths:
    input  = output_relax/snapshot_000
    output = ic_equilibrium.dat
"""
import struct
import sys
import os
import numpy as np

HEADER_SIZE = 256


def read_block(f):
    sz = struct.unpack('I', f.read(4))[0]
    data = f.read(sz)
    sz2 = struct.unpack('I', f.read(4))[0]
    if sz != sz2:
        raise ValueError(f"Block size mismatch: {sz} vs {sz2}")
    return data


def parse_header(hdr_buf):
    fmt = '<6I6d2d2i6I2i4d2i6Ii'
    sz = struct.calcsize(fmt)
    vals = struct.unpack(fmt, hdr_buf[:sz])
    h = {}
    h['npart'] = list(vals[0:6])
    h['mass'] = list(vals[6:12])
    h['time'] = vals[12]
    h['redshift'] = vals[13]
    h['flag_sfr'] = vals[14]
    h['flag_feedback'] = vals[15]
    h['npartTotal'] = list(vals[16:22])
    h['flag_cooling'] = vals[22]
    h['num_files'] = vals[23]
    h['BoxSize'] = vals[24]
    h['Omega0'] = vals[25]
    h['OmegaLambda'] = vals[26]
    h['HubbleParam'] = vals[27]
    h['flag_stellarage'] = vals[28]
    h['flag_metals'] = vals[29]
    h['npartTotalHW'] = list(vals[30:36])
    h['flag_entropy'] = vals[36]
    return h, hdr_buf[sz:]


def make_gadget2_header(h):
    """Rebuild the 256-byte header, forcing time=0 and redshift=0."""
    npart = np.array(h['npart'], dtype=np.uint32)
    mass = np.array(h['mass'], dtype=np.float64)
    fields = [
        ('6I', npart),
        ('6d', mass),
        ('d', 0.0),              # time = 0  (new IC)
        ('d', 0.0),              # redshift = 0
        ('i', h['flag_sfr']),
        ('i', h['flag_feedback']),
        ('6I', npart.copy()),    # npartTotal
        ('i', h['flag_cooling']),
        ('i', 1),                # num_files
        ('d', h['BoxSize']),
        ('d', h['Omega0']),
        ('d', h['OmegaLambda']),
        ('d', h['HubbleParam']),
        ('i', h['flag_stellarage']),
        ('i', h['flag_metals']),
        ('6I', np.array(h['npartTotalHW'], dtype=np.uint32)),
        ('i', h['flag_entropy']),
    ]
    buf = b''
    for fmt, data in fields:
        buf += struct.pack(fmt, *np.asarray(data).ravel())
    buf += b'\x00' * (HEADER_SIZE - len(buf))
    return buf


def write_block(f, name, data_bytes):
    f.write(struct.pack('I', len(data_bytes)))
    f.write(data_bytes)
    f.write(struct.pack('I', len(data_bytes)))


def main():
    in_path  = sys.argv[1] if len(sys.argv) > 1 else 'output_relax/snapshot_000'
    out_path = sys.argv[2] if len(sys.argv) > 2 else 'ic_equilibrium.dat'

    if not os.path.exists(in_path):
        print(f"ERROR: input snapshot not found: {in_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Converting {in_path} -> {out_path}")

    with open(in_path, 'rb') as fin:
        hdr_buf = read_block(fin)
        h, _ = parse_header(hdr_buf)

        n_total = int(sum(h['npart']))
        print(f"  Header: time={h['time']:.4f} -> 0.0")
        print(f"  npart = {h['npart']}, total = {n_total}")
        print(f"  mass  = {h['mass']}")

        # Read remaining blocks raw so we can copy them verbatim.
        # Gadget4 writes, in order: COORDS, VEL, IDS, (MASS if needed),
        # (U if gas), POT (if EVALPOTENTIAL).
        remaining_blocks = []
        while True:
            pos = fin.tell()
            four = fin.read(4)
            if len(four) < 4:
                break
            fin.seek(pos)
            try:
                blk = read_block(fin)
            except Exception:
                break
            remaining_blocks.append(blk)
            print(f"  Read block: {len(blk)} bytes")

    # Write out as IC: new header (time=0) + same data blocks.
    with open(out_path, 'wb') as fout:
        hdr = make_gadget2_header(h)
        fout.write(struct.pack('I', HEADER_SIZE))
        fout.write(hdr)
        fout.write(struct.pack('I', HEADER_SIZE))
        for blk in remaining_blocks:
            write_block(fout, 'data', blk)

    print(f"\nWrote {out_path}: {os.path.getsize(out_path)/1e6:.2f} MB")
    print("Header time reset to 0.0 — ready to use as InitCondFile for Phase B.")


if __name__ == '__main__':
    main()
