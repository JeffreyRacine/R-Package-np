# Conditional Plot Parity TODO

## Status

- `npcdens` and `npcdist` generalized/adaptive nearest-neighbor plot/bootstrap
  routes were repaired to restore `np`/`npRmpi` north-star parity in the
  audited smoke matrix.
- The remaining visible mismatch is concentrated in fixed-bandwidth
  local-polynomial conditional-density bootstrap (`npcdens`, `bwtype = "fixed"`,
  `regtype = "lp"`, especially bootstrap `inid`).

## Quick Diagnosis

- This residual does not appear to be a semantic routing divergence.
- It is concentrated on evaluation rows where the local-polynomial moment system
  is nearly singular; when the system is well-conditioned, `np` and `npRmpi`
  bootstrap rows agree to negligible tolerance.
- Existing forensic artifacts are under `/tmp`, notably:
  - `/tmp/target_np_lp_inid_20260323.rds`
  - `/tmp/target_nprmpi_lp_inid_20260323.rds`
  - `/tmp/reghelper_np_row181_20260323.rds`
  - `/tmp/reghelper_nprmpi_row181_20260323.rds`

## Working Interpretation

- The current evidence points to numerical amplification in the fixed LP
  bootstrap solve rather than a remaining plot/helper contract mismatch.
- If this is revisited, treat it as a stabilization problem for near-singular
  local-polynomial bootstrap systems, not as a generalized `npRmpi` plot parity
  failure.

## Revisit Trigger

- Reopen only if we want to add a numerically stabilized fixed-LP bootstrap
  solve/ridge policy that is explicitly matched to `np`, or if a new reproducible
  mismatch appears away from near-singular rows.
