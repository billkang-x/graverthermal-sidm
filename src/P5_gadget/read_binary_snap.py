#!/usr/bin/env python3
"""Reader for Gadget4 binary snapshot files (SnapFormat=1, GADGET2_HEADER).

Binary format:
  For each block:
    [4-byte int N] [N bytes of data] [4-byte int N]

  Header block (256 bytes):
    npart[6] int32      : particle count per type
    mass[6]  double     : mass per type (0 means per-particle)
    time     double
    redshift double
    flag_sfr int32
    flag_feedback int32
    npartTotal[6] uint32  : total num particles (low word)
    flag_cooling int32
    num_files int32
    BoxSize  double
    Omega0   double
    OmegaLambda double
    HubbleParam double
    flag_stellarage int32
    flag_metals int32
    npartTotalHW[6] uint32  : high word for npartTotal
    flag_entropy int32
    (filler to 256 bytes)

  Data blocks (in order):
    Coordinates  npart_total x 3   float32
    Velocities   npart_total x 3   float32
    ParticleIDs   npart_total       uint32
    Masses       (only types with MassTable==0)   float32
    InternalEnergy (gas only)      float32
    (optional) Potential, Acceleration, etc.
"""
import struct
import numpy as np

HEADER_SIZE = 256
NTYPES = 6  # Gadget2 header uses 6 types


def read_block(f):
    """Read a Fortran-style block: [int N][N bytes][int N]."""
    sz = struct.unpack('I', f.read(4))[0]
    data = f.read(sz)
    sz2 = struct.unpack('I', f.read(4))[0]
    if sz != sz2:
        raise ValueError(f"Block size mismatch: {sz} vs {sz2}")
    return data


def parse_header(hdr_buf):
    """Parse 256-byte Gadget2 header."""
    if len(hdr_buf) < HEADER_SIZE:
        raise ValueError(f"Header too small: {len(hdr_buf)}")
    # Layout (matching gen_ic_binary.py):
    #   npart[6] uint32, mass[6] double, time double, redshift double,
    #   flag_sfr int32, flag_feedback int32,
    #   npartTotal[6] uint32, flag_cooling int32, num_files int32,
    #   BoxSize double, Omega0 double, OmegaLambda double, HubbleParam double,
    #   flag_stellarage int32, flag_metals int32,
    #   npartTotalHW[6] uint32, flag_entropy int32, padding to 256 bytes
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


def read_binary_snapshot(path, types_to_read=(1,)):
    """Read a Gadget4 binary snapshot.

    Returns dict with: header, coords[type] -> (N,3) array,
                       vels[type] -> (N,3) array,
                       ids[type] -> (N,) array,
                       masses[type] -> (N,) array
    """
    with open(path, 'rb') as f:
        # Header block
        hdr_buf = read_block(f)
        h = parse_header(hdr_buf)
        print(f"  Header: time={h['time']:.4f}, redshift={h['redshift']:.4f}")
        print(f"    npart={list(h['npart'])}")
        print(f"    mass={list(h['mass'])}")
        print(f"    npartTotal={list(h['npartTotal'])}")
        print(f"    BoxSize={h['BoxSize']}")

        npart = list(h['npart'])
        # Sanity check: some Gadget4 snapshots have garbage in npart[3:6]
        # due to header format differences. Only keep types 0-2 (gas, DM, gas2).
        # Also zero out small garbage values in type 3-5.
        for i in range(3, 6):
            npart[i] = 0
            h['npart'][i] = 0
        npart_total = sum(npart)
        if npart_total == 0:
            raise ValueError("Empty snapshot (npart=0)")

        # Coordinates
        coords_buf = read_block(f)
        coords_all = np.frombuffer(coords_buf, dtype='<f4')
        # Verify: coords should be npart_total * 3
        expected_floats = npart_total * 3
        if len(coords_all) != expected_floats:
            # Try to infer actual npart from data size
            actual_npart = len(coords_all) // 3
            if actual_npart * 3 == len(coords_all):
                print(f"  Warning: header npart_total={npart_total} but data has {actual_npart} particles")
                print(f"    Adjusting npart to match data. This may affect type slicing.")
                # Redistribute: keep types 0-2 as-is, zero everything else
                npart_total = actual_npart
        coords_all = coords_all.reshape(npart_total, 3)
        # Velocities
        vels_buf = read_block(f)
        vels_all = np.frombuffer(vels_buf, dtype='<f4').reshape(npart_total, 3)
        # IDs
        ids_buf = read_block(f)
        ids_all = np.frombuffer(ids_buf, dtype='<u4')

        # Per-type slicing
        result = {'header': h}
        offset = 0
        for t in range(NTYPES):
            n = npart[t]
            if n == 0:
                continue
            c_slice = slice(offset, offset + n)
            result.setdefault('coords', {})[t] = np.array(coords_all[c_slice], dtype=np.float64)
            result.setdefault('vels', {})[t] = np.array(vels_all[c_slice], dtype=np.float64)
            if len(ids_all) == npart_total:
                result.setdefault('ids', {})[t] = np.array(ids_all[c_slice])
            # Per-type masses: either from header or per-particle
            if h['mass'][t] > 0:
                result.setdefault('masses', {})[t] = np.full(n, h['mass'][t], dtype=np.float64)
            offset += n

        # If any type has MassTable==0, there is a Masses block next
        types_with_per_particle_mass = [t for t in range(NTYPES)
                                         if npart[t] > 0 and h['mass'][t] == 0]
        if types_with_per_particle_mass:
            try:
                mass_buf = read_block(f)
                mass_all = np.frombuffer(mass_buf, dtype='<f4')
                # Masses are written per type, in order, but only for types with MassTable==0
                offset = 0
                for t in range(NTYPES):
                    n = npart[t]
                    if n == 0:
                        continue
                    if h['mass'][t] == 0:
                        result.setdefault('masses', {})[t] = np.array(
                            mass_all[offset:offset+n], dtype=np.float64)
                        offset += n
            except Exception as e:
                print(f"  Warning: could not read Masses block: {e}")

        return result


def read_snapshot(path, types_to_read=(1,)):
    """Auto-detect HDF5 vs binary by file header bytes."""
    with open(path, 'rb') as f:
        magic = f.read(8)
    if magic[:8] == b'\x89HDF\r\n\x1a\n':
        # HDF5
        import h5py
        with h5py.File(path, 'r') as ff:
            h = dict(ff['Header'].attrs)
            result = {'header': h}
            for t in types_to_read:
                key = f'PartType{t}'
                if key in ff:
                    grp = ff[key]
                    result.setdefault('coords', {})[t] = grp['Coordinates'][:]
                    result.setdefault('vels', {})[t] = grp['Velocities'][:]
                    if 'Masses' in grp:
                        result.setdefault('masses', {})[t] = grp['Masses'][:]
                    elif 'MassTable' in h and h['MassTable'][t] > 0:
                        result.setdefault('masses', {})[t] = np.full(len(grp['Coordinates']),
                                                                      h['MassTable'][t])
                    if 'ParticleIDs' in grp:
                        result.setdefault('ids', {})[t] = grp['ParticleIDs'][:]
            return result
    else:
        # Binary Gadget2 format - first 4 bytes are block size = 256
        return read_binary_snapshot(path, types_to_read)


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print("Usage: read_binary_snap.py <snapshot_file>")
        sys.exit(1)
    snap = read_snapshot(sys.argv[1])
    h = snap['header']
    print(f"\nSummary:")
    for t in sorted(snap.get('coords', {}).keys()):
        c = snap['coords'][t]
        m = snap.get('masses', {}).get(t, None)
        print(f"  PartType{t}: N={len(c)}, mass sum={m.sum() if m is not None else 'N/A'}")
        r = np.sqrt((c**2).sum(axis=1))
        print(f"    r: min={r.min():.4f}, max={r.max():.4f}, mean={r.mean():.4f}")
