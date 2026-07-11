# MPI setup and runtime modes

npRmpi 0.70-5 is release-validated with MPICH. Open MPI runtime validation is
explicitly waived for this release; an Open-MPI-linked build should not be
treated as having the same supported runtime contract.

## 1. Interactive session/spawn mode

Use session mode when working in an interactive R process and allowing npRmpi
to spawn its worker pool:

```r
library(npRmpi)
npRmpi.init(nslaves = 1)

# ... run np* calls ...

npRmpi.quit(force = TRUE)
```

On macOS, npRmpi enables worker reuse by default. Consequently,
`npRmpi.quit()` without `force = TRUE` may retain the pool for another call to
`npRmpi.init()`. Use `force = TRUE` when the R session is finished or when a
clean lifecycle test is required.

`mode = "auto"` selects attach mode when the process is already part of an MPI
world with more than one rank; otherwise it selects spawn mode.

## 2. Cluster/batch attach mode

Prelaunch the MPI world:

```bash
mpiexec -n 4 Rscript foo.R
```

Inside `foo.R`:

```r
library(npRmpi)
is_master <- npRmpi.init(mode = "attach", autodispatch = TRUE,
                         np.messages = FALSE)

if (isTRUE(is_master)) {
  # ... run np* calls ...
  npRmpi.quit(mode = "attach")
}
```

Non-master ranks enter the package worker loop. The estimator body and attach
close belong on the master rank. Use a fresh `mpiexec` launch for each
independent attach job and verify that no job-scoped MPI/R processes remain.

## 3. Profile/manual-broadcast mode

Profile mode is available for workflows that explicitly broadcast commands.
Set `R_PROFILE_USER` to the installed npRmpi `Rprofile`, launch the desired MPI
world, and use `np.mpi.initialize()` plus `mpi.bcast.*` calls in the script.
Do not also call `npRmpi.init()` inside that profile-managed workflow.

## 4. MPICH network-interface diagnostics

Most local MPICH runs need no interface override. When MPICH/OFI diagnostics
show that an explicit interface is required, apply it to the command being
tested rather than treating it as a package default:

```bash
FI_TCP_IFACE=en0 Rscript your-script.R
```

Interface names are host-specific. A loopback-only diagnostic can use:

```bash
FI_TCP_IFACE=lo0 FI_PROVIDER=tcp FI_SOCKETS_IFACE=lo0 \
  Rscript your-script.R
```

Loopback settings are not multi-host release evidence.

## 5. Cleanup and positive claims

- Run each backend and independent attach/profile world in a fresh R process.
- Use `force = TRUE` when a spawned session must leave no reusable pool.
- After failures or timeouts, identify and clean only route-scoped processes.
- Do not call a route validated merely because package loading or compilation
  succeeded; use an installed tarball and substantive runtime sentinels.
- A teardown exit 137 is ignorable only after substantive pass output and proof
  that no route-scoped workers remain.
