#!/bin/bash
# 编译环境 - 与ode28项目一致
if [[ -f /public1/soft/modules/module.sh ]]; then
  source /public1/soft/modules/module.sh
fi

module purge || true
module load gcc/12.2 mpi/openmpi/4.1.5-gcc12.2 gsl/2.0 hdf5/1.8.13-gcc-zyq fftw/3.3.8-mpi

# 关键：强制CC和CXX使用MPI编译器，避免回退到系统cc
export CC=mpicc
export CXX=mpicxx
export FC=mpif90

cd ~/gas_sidm_bars/src/gadget4

echo "=== Build Environment ==="
echo "CC: $CC -> $(which $CC)"
echo "CXX: $CXX -> $(which $CXX)"
echo "GCC version: $(g++ --version | head -1)"
echo ""

echo "=== Cleaning build ==="
make clean > /dev/null 2>&1

echo "=== Starting make -j8 ==="
make -j8 2>&1

echo ""
echo "=== Build Result ==="
if [ -f Gadget4 ]; then
  ls -lh Gadget4
  echo "SUCCESS: Gadget4 compiled!"
else
  echo "FAILED: Gadget4 not found"
fi
