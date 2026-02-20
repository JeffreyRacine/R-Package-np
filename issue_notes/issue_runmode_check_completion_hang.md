# Run-mode `R CMD check --as-cran` completion instability

Observed repeatedly during `man/run` mode checks in `npRmpi`:

- `npRmpi-Ex.Rout` reaches footer (`### * <FOOTER> ... quit('no')`) indicating examples completed.
- `00check.log` can remain truncated at `* checking examples ...` and the check driver may terminate unexpectedly.

This appears related to a known interaction between `Rmpi` and `npRmpi` in check-time process management, not a specific single example failure once current doc patches are applied.

TODO:

1. Reproduce with a minimal package/examples harness outside `npRmpi`.
2. Compare behavior with and without slave reuse (`npRmpi.reuse.slaves`).
3. Compare `R CMD check` with `--no-examples` + standalone examples run to isolate teardown stage.
4. Capture and inspect master/slave process exit paths during check finalization.
