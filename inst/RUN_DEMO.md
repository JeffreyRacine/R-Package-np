# Demo Run Guide

This directory supports three execution modes:

- `serial`: plain `np` scripts (`*_serial.R`)
- `attach`: `mpiexec` + `npRmpi.init(mode="attach", ...)` scripts (`*_npRmpi_attach.R`)
- `profile`: `mpiexec` + manual broadcast scripts (`*_npRmpi_profile.R`) that require `Rprofile`

## One-time setup for profile mode

The package ships the profile template at:

- `inst/Rprofile`

Use one of these options before running profile mode:

1. Copy to the job working directory as `.Rprofile`.
2. Or set `R_PROFILE_USER` to that file path.

Example:

```bash
R_PROFILE_USER="$(Rscript -e 'cat(system.file("Rprofile", package="npRmpi"))')"
export R_PROFILE_USER
```

## Run everything

From this directory:

```bash
./runall
```

This runs:

1. `serial`
2. `attach` for `NP=2,3,4`
3. `profile` for `NP=2,3,4`

Results are written under:

- `serial/`
- `n_2_attach/`, `n_3_attach/`, `n_4_attach/`
- `n_2_profile/`, `n_3_profile/`, `n_4_profile/`

## Run a single mode manually

Examples:

```bash
mkdir -p serial && (cd serial && make -f ../makefile MODE=serial)
mkdir -p n_2_attach && (cd n_2_attach && make -f ../makefile MODE=attach NP=2)
mkdir -p n_2_profile && (cd n_2_profile && make -f ../makefile MODE=profile NP=2)
```

## Tiny smoke runs

You can reduce sample size in many demos via:

```bash
NP_DEMO_N=80
```

Example:

```bash
mkdir -p n_2_profile
cd n_2_profile
NP_DEMO_N=80 make -f ../makefile MODE=profile NP=2 DEMOS=npcdensml
```
