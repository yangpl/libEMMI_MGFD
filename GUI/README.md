# libEMMI_MGFD GTK Launcher

This is a GTK3 interface to run `libEMMI_MGFD` with an editable `run.sh`, following the same general workflow as the `SMIwiz` launcher.

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

- Select project directory, run directory, and `run.sh`.
- Load and edit `run.sh` directly in the main text panel.
- Edit common modelling and inversion parameters in a form.
- Sync controls:
  - `Script -> Form` parses `key=value` lines from the script into form fields.
  - `Form -> Script` writes form values back into the script text.
- Save `run.sh` from the GUI.
- Run `bash run.sh` directly from the selected run directory while streaming stdout/stderr to the log panel.
- Stop a running job with `SIGTERM`.
- Run an optional plot command manually in the selected run directory.
