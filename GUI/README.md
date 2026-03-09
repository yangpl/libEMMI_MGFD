# libEMMI_MGFD GTK Launcher

This is a GTK3 interface to run `libEMMI_MGFD` by editing `inputpar.txt` and launching the executable directly with MPI.

## Build

```bash
cd GUI
make
```

## Run

```bash
cd GUI
./gui
```

## What It Does

- Select project directory, run directory, and the `bin/main` executable.
- Edit `inputpar.txt` directly in the main text panel.
- Load and save `inputpar.txt` from the GUI.
- Set the MPI rank count in the GUI.
- Run the solver with:
  ```bash
  mpirun -n <ranks> <binary> $(cat inputpar.txt)
  ```
  in the selected run directory.
- Stream stdout and stderr to the log panel.
- Stop a running job with `SIGTERM`.
- Run an optional plot command manually in the selected run directory.

## Default Paths

By default the launcher points to:

- project directory: `libEMMI_MGFD`
- run directory: `libEMMI_MGFD/run_modelling`
- binary: `libEMMI_MGFD/bin/main`
