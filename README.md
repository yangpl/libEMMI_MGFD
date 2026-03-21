# libEMMI_MGFD

`libEMMI_MGFD` is a 3D marine controlled-source electromagnetic (CSEM) modelling and inversion code based on a multigrid frequency-domain solver. The main executable supports:

- forward modelling (`mode=0`)
- full-waveform inversion (`mode=1`)
- gradient evaluation (`mode=2`)

The implementation is written in C and uses MPI so that different shots can be distributed across processes.

## Reference

Pengliang Yang and An Ping, "libEMMI MGFD: A program of marine controlled-source electromagnetic modelling and inversion using frequency-domain multigrid solver", *Computer Physics Communications* 305 (2024) 109327.  
DOI: [10.1016/j.cpc.2024.109327](https://doi.org/10.1016/j.cpc.2024.109327)

## Repository Layout

- `src/`: core solver, inversion, and I/O code
- `include/`: public headers for model, acquisition, and optimization structures
- `bin/`: compiled executable (`main`)
- `run_modelling/`: forward-modelling example and plotting scripts
- `run_inversion/`: inversion example and plotting scripts
- `run_layered/`: layered-model helper scripts
- `create_model3d/`: utilities for preparing simple 3D models
- `doc/`: paper-related notes and documentation sources
- `GUI/`: separate GUI prototype

## Requirements

- MPI C compiler and runtime, for example OpenMPI or MPICH
- standard C math library
- optional: Python 3 for plotting and comparison scripts in the example folders

The build uses `mpicc` by default.

## Build

Compile the main program from [`src/Makefile`](/home/pyang/Documents/libEMMI_MGFD/src/Makefile):

```bash
cd src
make
```

This creates [`bin/main`](/home/pyang/Documents/libEMMI_MGFD/bin/main).

## Running The Solver

The executable reads parameters from command-line `key=value` pairs:

```bash
mpirun -n <nproc> ../bin/main key=value key=value ...
```

Each MPI rank handles one shot by default. If `nproc` is larger than the number of available sources, the code aborts unless you provide an explicit `shots=` list whose length matches `nproc`.

## Required Inputs

At minimum, the solver expects:

- a resistivity model for `rho11`, `rho22`, and `rho33` via `frho11=`, `frho22=`, `frho33=`
- source geometry file `fsrc=`
- receiver geometry file `frec=`
- source-receiver connection table `fsrcrec=`
- mesh and survey extents (`nx`, `ny`, `nz`, `dx`, `dy`, `dz`, `ox`, `oy`, `oz`)
- at least one frequency in `freqs=`

### Resistivity Model Format

`frho11`, `frho22`, and `frho33` are raw binary files of `float32` values with exactly `nx*ny*nz` samples each. The code reads them in C order as one value per cell.

For isotropic runs, point all three parameters to the same file.

### Acquisition File Formats

The example files in [`run_modelling/`](/home/pyang/Documents/libEMMI_MGFD/run_modelling) show the expected ASCII layout.

`sources.txt`

```text
x y z azimuth dip isrc
```

`receivers.txt`

```text
x y z azimuth dip irec
```

`src_rec_table.txt`

```text
isrc irec
```

Coordinates are in meters. `azimuth` and `dip` describe source or receiver orientation. The first line is treated as a header and skipped.

## Common Parameters

The code accepts many runtime parameters. The most commonly used ones are:

- `mode`: `0` modelling, `1` inversion, `2` gradient only
- `freqs`: comma-separated frequencies in Hz
- `chsrc`: active source channels, for example `Ex`, `Ey`, `Ez`
- `chrec`: active receiver channels, for example `Ex,Ey,Hx,Hy`
- `nx,ny,nz`: number of cells in the input model
- `dx,dy,dz`: input cell sizes
- `ox,oy,oz`: model origin
- `n1,n2,n3`: multigrid mesh size used internally
- `icheck=1`: also write full field volumes `Ex.bin`, `Ey.bin`, `Ez.bin`, `Hx.bin`, `Hy.bin`, `Hz.bin`
- `shots`: optional explicit list of source indices assigned to MPI ranks

Additional modelling and inversion controls are implemented in [`src/emf_init_free.c`](/home/pyang/Documents/libEMMI_MGFD/src/emf_init_free.c), [`src/acq_init_free.c`](/home/pyang/Documents/libEMMI_MGFD/src/acq_init_free.c), [`src/do_modelling.c`](/home/pyang/Documents/libEMMI_MGFD/src/do_modelling.c), and [`src/do_fwi.c`](/home/pyang/Documents/libEMMI_MGFD/src/do_fwi.c).

## Forward-Modelling Example

[`run_modelling/run.sh`](/home/pyang/Documents/libEMMI_MGFD/run_modelling/run.sh) is the simplest starting point:

```bash
cd run_modelling
bash run.sh
```

That script:

- writes an `inputpar.txt` file with a small test setup
- runs `mpirun -n 1 ../bin/main ...`
- produces `emf_0001.txt`
- with `icheck=1`, also writes `Ex.bin`, `Ey.bin`, `Ez.bin`, `Hx.bin`, `Hy.bin`, `Hz.bin`

The output CSEM data are written as ASCII files named:

```text
emf_0001.txt
emf_0002.txt
...
```

one file per shot / MPI rank. Each file contains:

```text
iTx iRx chrec frequency/Hz Real{E/H} Imag{E/H}
```

## Inversion Example

[`run_inversion/run.sh`](/home/pyang/Documents/libEMMI_MGFD/run_inversion/run.sh) demonstrates a multi-shot inversion:

```bash
cd run_inversion
bash run.sh
```

Important inversion-specific inputs:

- `fbathy=`: bathymetry grid in binary format
- `niter=`: maximum number of iterations
- `preco=`: preconditioning mode
- `npar=` and `idxpar=`: number and identity of inverted parameters
- `bound=1`, `minpar=`, `maxpar=`: bounds in physical resistivity
- `gamma1=` and `gamma2=`: Tikhonov and TV regularization weights
- `dsmute=` and `drmute=`: source/receiver muting radii for the gradient

Typical inversion outputs include:

- `iterate.txt`: iteration log
- `rmse_misfit.txt`: misfit history
- `rmse_XXXX.txt`: per-shot residual summary
- `Rh_iteration_XX`, `Rv_iteration_XX`: binary model snapshots by iteration

The inversion code expects observed data in shot-based ASCII files named `emf_XXXX.txt`, using the same layout as the forward-modelling output.

## Reproducing The Bundled Inversion Test

Inside [`run_inversion/`](/home/pyang/Documents/libEMMI_MGFD/run_inversion), the included files support a synthetic workflow:

1. Generate observed data from the true model by running with `mode=0` and `frho11=frho frho22=frho frho33=frho`.
2. Start inversion from the initial model by running with `mode=1` and `frho11=frho_init frho22=frho_init frho33=frho_init`.

The example `run.sh` already contains a configured inversion command; adjust it if you want to regenerate synthetic observations first.

## Notes

- The code is MPI-oriented. In practice, set `nproc` to the number of shots you want to process in parallel.
- If `icheck=1`, full-field binary outputs can be large.
- The top-level worktree currently contains generated files and build artifacts; example directories are intended to be run in place.

## Acknowledgements

Pengliang Yang acknowledges Dieter Werthmüller for sustained feedback during the development of the multigrid CSEM code and for the related `emg3d` ecosystem used for comparison.

## Contact

Pengliang Yang   
Email: `ypl.2100@gmail.com`
