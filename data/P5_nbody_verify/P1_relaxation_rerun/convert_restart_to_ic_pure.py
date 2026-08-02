#!/usr/bin/env python3
"""Convert a Gadget4 binary snapshot into a binary IC (no numpy dependency).

Pure-Python version using only struct + os, so it runs on HPC's bare python3.
"""
import struct
import sys
import os

HEADER_SIZE = 256


def read_block(f):
    sz = struct.unpack('I', f.read(4))[0]
    data = f.read(sz)
    sz2 = struct.unpack('I', f.read(4))[0]
    if sz != sz2:
        raise ValueError("Block size mismatch: %d vs %d" % (sz, sz2))
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
    return h


def make_gadget2_header(h):
    """Rebuild 256-byte header, forcing time=0 and redshift=0.

    Snapshot headers may contain garbage in fields beyond the standard IC
    layout (Gadget4 stores metadata in padding). We force all flags to 0
    and use npart for npartTotal to produce a clean IC header.
    """
    import struct as st
    buf = b''
    # npart[6] uint32 - only keep types 0-2, zero out 3-5 (snapshot garbage)
    npart_clean = list(h['npart'][:3]) + [0, 0, 0]
    buf += st.pack('<6I', *npart_clean)
    # mass[6] double - only keep types 0-2, zero out 3-5
    mass_clean = list(h['mass'][:3]) + [0.0, 0.0, 0.0]
    buf += st.pack('<6d', *mass_clean)
    # time=0, redshift=0
    buf += st.pack('<2d', 0.0, 0.0)
    # flag_sfr=0, flag_feedback=0
    buf += st.pack('<2i', 0, 0)
    # npartTotal[6] = npart_clean (clean copy)
    buf += st.pack('<6I', *npart_clean)
    # flag_cooling=0, num_files=1
    buf += st.pack('<2i', 0, 1)
    # BoxSize, Omega0=0, OmegaLambda=0, HubbleParam=1.0
    buf += st.pack('<4d', h['BoxSize'], 0.0, 0.0, 1.0)
    # flag_stellarage=0, flag_metals=0
    buf += st.pack('<2i', 0, 0)
    # npartTotalHW[6] = all zeros
    buf += st.pack('<6I', 0, 0, 0, 0, 0, 0)
    # flag_entropy=0 (CRITICAL: must be 0 for IC)
    buf += st.pack('<i', 0)
    # pad to 256 with zeros (clears any snapshot metadata in padding)
    buf += b'\x00' * (HEADER_SIZE - len(buf))
    return buf


def write_block(f, data_bytes):
    f.write(struct.pack('I', len(data_bytes)))
    f.write(data_bytes)
    f.write(struct.pack('I', len(data_bytes)))


def main():
    in_path  = sys.argv[1] if len(sys.argv) > 1 else 'output_relax/snapshot_000'
    out_path = sys.argv[2] if len(sys.argv) > 2 else 'ic_equilibrium.dat'

    if not os.path.exists(in_path):
        print("ERROR: input snapshot not found: %s" % in_path, file=sys.stderr)
        sys.exit(1)

    print("Converting %s -> %s" % (in_path, out_path))

    with open(in_path, 'rb') as fin:
        hdr_buf = read_block(fin)
        h = parse_header(hdr_buf)

        n_total = int(sum(h['npart']))
        print("  Header: time=%.4f -> 0.0" % h['time'])
        print("  npart = %s, total = %d" % (h['npart'], n_total))
        print("  mass  = %s" % h['mass'])

        # Read remaining blocks raw
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
            print("  Read block: %d bytes" % len(blk))

    # Write out as IC: new header (time=0) + same data blocks
    with open(out_path, 'wb') as fout:
        hdr = make_gadget2_header(h)
        fout.write(struct.pack('I', HEADER_SIZE))
        fout.write(hdr)
        fout.write(struct.pack('I', HEADER_SIZE))
        for blk in remaining_blocks:
            write_block(fout, blk)

    fsize = os.path.getsize(out_path)
    print("\nWrote %s: %.2f MB" % (out_path, fsize / 1e6))
    print("Header time reset to 0.0 - ready to use as InitCondFile for Phase B.")


if __name__ == '__main__':
    main()
