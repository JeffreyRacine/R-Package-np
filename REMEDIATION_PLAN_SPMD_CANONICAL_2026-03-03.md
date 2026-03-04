# npRmpi Canonical Remediation Plan: Session SPMD Unification (2026-03-03)

## Canonical Status
1. This is the active remediation plan for session/attach/profile unification.
2. `nslaves=0` is removed permanently from `npRmpi`.
3. `npRmpi` runtime contract is MPI-only (`nslaves>=1`).
4. Remediation implementation phases are complete; follow-on work is maintenance/performance tuning under the same gate policy.

## Execution Status (Updated 2026-03-03/04)
Completed checkpoint tranches:
1. Phase 1 (`nslaves=0` removal): complete (`2fc9317`).
2. Phase 2 (SPMD control-plane scaffolding): complete (`905ffd9`, `6b6be81`, `1fece7e`, `3bf4914`).
3. Phase 3 (`npregbw` LL/LP CV migration + `O(n^2)` helper excision): complete (`24dbd4b` and subsequent route/contract gates).
4. Phase 4 (core-family migration progress): complete for `npudens*`, `npudist*`, `npcdens*`, `npcdist*`, `npscoef*`, `npindex*`, `npplreg*`.
5. Phase 5 (hot-path legacy asymmetry decommission): complete for core estimator/CV hot paths via typed locked opcode handlers and guard-validated route parity.
6. Phase 6 (non-core migration): tranche A complete for `npqreg*`, `npconmode*`, `npksum*`, `npregiv*`, `npregivderiv*`, `npcmstest`, `npqcmstest`, `npdeneqtest`, `npdeptest`, `npsdeptest`, `npsymtest`, and `npunitest` via typed locked opcodes.
7. Phase 6 (non-core migration): tranche B complete for `npsigtest*` typed locked opcode routing plus session/attach/profile fast-fail contract coverage.
8. Phase 6 closeout audit: complete for active autodispatch call heads; no remaining `manual_distributed_call` usage in estimator/test paths.
9. Phase 6 closeout reliability hardening: complete for session lifecycle reseed/sequence reset, fixed-length worker type packets, LC CVLS stable fast-path routing, and `npsigtest` invalid-index fast-fail pre-dispatch validation.
10. Phase 6 final gate closeout: complete (all-flags session route subprocess contract, session/attach/profile/manual validators, serial `np` anchor, phase-23 paired-seed performance/equivalence evidence, and tarball-first `R CMD build/check --as-cran` with expected legacy NOTE/WARNING only).
11. Phase 6 post-closeout revalidation: complete (phase-24 rerun of session/attach/profile/manual route validators with clean orphan-process scans).
12. Phase 6 terminal confirmation: complete (phase-25 rerun of preflight scans, full route matrix, all-flags subprocess contract tests, serial `np-master` anchor smoke, and tarball-first `R CMD build/check --as-cran` with expected legacy NOTE/WARNING only and no check errors).
13. Post-closeout packaging hygiene: complete (`.Rbuildignore` now excludes top-level remediation markdown files from source tarball; `R CMD check --as-cran` no longer reports non-standard top-level files warning; remaining warning is environment-side `checkbashisms` availability only).
14. Local as-cran tooling closure: complete for this machine using a shellcheck-backed `checkbashisms` compatibility shim on `PATH` (scoped to real shell scripts), yielding `R CMD check --as-cran` status `1 NOTE` only (CRAN incoming version-jump note).
15. Phase-26 assurance rerun: complete (fresh preflight scans, session/attach/profile/manual route validators, all-flags subprocess contract pass, serial `np-master` anchor pass, plain `--as-cran` status `1 WARNING, 1 NOTE`, shimmed `--as-cran` status `1 NOTE`, and clean orphan scans).
16. Phase-27 full assurance rerun: complete in a fresh artifact root with repeated route matrix pass (`session`/`attach`/`profile`/`manual`), subprocess contract pass (`PASS 38`), serial `np-master` anchor pass, plain vs shimmed `--as-cran` status confirmation (`1 WARNING, 1 NOTE` vs `1 NOTE`), and clean orphan scans.
17. Phase-28 invariant and non-core contract closeout: complete (static scans for `nslaves=0`, `manual_distributed_call`, and `O(n^2)` helper references; patched `test-spmd-step-contract.R` collective-step subprocess tests to run through broadcasted SPMD entry; non-core contract tests pass with zero nonzero-fail markers; full route matrix pass; plain vs shimmed `--as-cran` status confirmation; clean orphan scans).
18. Phase-29 terminal reconfirmation: complete (fresh invariant scans, patched non-core contract tests with zero nonzero-fail markers, route matrix pass, subprocess contract `PASS 38`, serial `np-master` anchor pass, shimmed `--as-cran` status `1 NOTE`, plain `--as-cran` status in expected warning+note band, and clean orphan scans).
19. Phase-30 deterministic local check harness: complete (`issue_notes/run_as_cran_local_clean.sh` added; shellcheck-clean; produces reproducible local tarball-first `--as-cran` status `1 NOTE` with `future file timestamps ... OK` and `top-level files ... OK` on this machine).
20. Phase-31 final gate pack: complete (hardened `issue_notes/run_spmd_final_gate_pack.sh` with heartbeat-safe step execution, fresh full gate matrix pass, `FINAL_GATE_PACK_OK` token, and clean orphan scans in a new artifact root).
21. Phase-32 terminal rerun: complete (fresh preflight, full route matrix and contract pass markers, serial anchor pass, shimmed note-only as-cran status, and `FINAL_GATE_PACK_OK` token in a new artifact root).
22. Phase-33 archival rerun: complete (fresh archival gate-pack rerun with full route/contract/anchor/as-cran checks and `FINAL_GATE_PACK_OK` token).
23. Final evidence index generated: `/tmp/spmd_canonical_20260304_0001/FINAL_EVIDENCE_INDEX_20260304.md` (aggregated artifact existence + terminal closure tokens).
24. Permanent handoff record added: `/Users/jracine/Development/np-npRmpi/SPMD_REMEDIATION_FINAL_HANDOFF_2026-03-04.md`.
25. Local closure milestone tag created: `spmd-remediation-complete-2026-03-04` -> `abde3bf`.
26. Closure certificate added: `/Users/jracine/Development/np-npRmpi/SPMD_REMEDIATION_CLOSURE_CERTIFICATE_2026-03-04.md`.
27. Machine-readable closure manifest added: `/Users/jracine/Development/np-npRmpi/SPMD_REMEDIATION_CLOSURE_MANIFEST_2026-03-04.json`.
28. Post-certificate local tag created: `spmd-remediation-post-cert-closeout-2026-03-04` -> `6e86ce5`.
29. Manifest verifier script added and executed: `/Users/jracine/Development/np-npRmpi/issue_notes/verify_spmd_closure_manifest.sh` -> `/tmp/spmd_canonical_20260304_0001/phase34_manifest_verify_20260304_013133` with `CLOSURE_MANIFEST_VERIFY_OK`.
30. Portable closure bundle exporter added and executed: `/Users/jracine/Development/np-npRmpi/issue_notes/export_spmd_closure_bundle.sh` -> `/tmp/spmd_canonical_20260304_0001/phase35_closure_bundle_20260304_013242` with `CLOSURE_BUNDLE_EXPORT_OK`.
31. Full closure audit runner added and executed: `/Users/jracine/Development/np-npRmpi/issue_notes/run_spmd_closure_full_audit.sh` -> `/tmp/spmd_canonical_20260304_0001/phase36_full_audit_20260304_013355` with `SPMD_CLOSURE_FULL_AUDIT_OK`.
32. Compressed closure archive packager added and executed: `/Users/jracine/Development/np-npRmpi/issue_notes/package_spmd_closure_archive.sh` -> `/tmp/spmd_canonical_20260304_0001/phase37_closure_archive_20260304_013513` with `CLOSURE_ARCHIVE_PACKAGE_OK`.
33. Final no-code reproducibility rerun complete: full audit rerun at `/tmp/spmd_canonical_20260304_0001/phase38_full_audit_20260304_013613` with `SPMD_CLOSURE_FULL_AUDIT_OK` and downstream token chain confirmation.
34. Repo-local closeout workflow note added: `/Users/jracine/Development/np-npRmpi/SPMD_REMEDIATION_CLOSEOUT_WORKFLOW_2026-03-04.md`.
35. Final-freeze local tag created: `spmd-remediation-final-freeze-2026-03-04` -> `870762c`.
36. Closure manifest/certificate refreshed to current head and tag map.
37. Manifest verification rerun complete: `/tmp/spmd_canonical_20260304_0001/phase40_manifest_verify_20260304_013945` with `CLOSURE_MANIFEST_VERIFY_OK`.
38. Full-audit rerun complete after refresh: `/tmp/spmd_canonical_20260304_0001/phase41_full_audit_20260304_013952` with `SPMD_CLOSURE_FULL_AUDIT_OK`.
39. Phase-47 stagewise terminal gate-pack rerun complete in a fresh artifact root with `FINAL_GATE_PACK_OK`: `/tmp/spmd_canonical_20260304_0001/phase47_stagewise_gatepack_20260304_014609`.
40. Packaging-hygiene follow-up complete: `.Rbuildignore` updated to exclude top-level `SPMD_REMEDIATION_*` closure files from source tarball; local deterministic as-cran harness returned `Status: 1 NOTE` (CRAN incoming/version note only).
41. Manifest verification rerun complete after phase-47 updates: `/tmp/spmd_canonical_20260304_0001/phase48_manifest_verify_20260304_015705` with `CLOSURE_MANIFEST_VERIFY_OK`.
42. Full-audit rerun complete after phase-47 updates: `/tmp/spmd_canonical_20260304_0001/phase49_full_audit_20260304_015705` with `SPMD_CLOSURE_FULL_AUDIT_OK`.
43. Closure-manifest verifier hardened to validate all manifest-declared phases (including phase47 with configurable as-cran log relative path), plus robust `FAILNZ_*` whitespace-tolerant zero checks.
44. Dynamic verifier/fullaudit reruns complete after verifier hardening: `/tmp/spmd_canonical_20260304_0001/phase53_manifest_verify_dynamic_20260304_015851` (`CLOSURE_MANIFEST_VERIFY_OK`) and `/tmp/spmd_canonical_20260304_0001/phase54_full_audit_dynamic_20260304_015854` (`SPMD_CLOSURE_FULL_AUDIT_OK`).
45. Final local closure anchor tag created: `spmd-remediation-final-closure-2026-03-04` -> `006127a`.
46. Post-tag closure reruns complete: `/tmp/spmd_canonical_20260304_0001/phase57_manifest_verify_tagged_20260304_020048` (`CLOSURE_MANIFEST_VERIFY_OK`) and `/tmp/spmd_canonical_20260304_0001/phase58_full_audit_tagged_20260304_020052` (`SPMD_CLOSURE_FULL_AUDIT_OK`).
47. Closure manifest updated to include `final_closure` tag metadata (`spmd-remediation-final-closure-2026-03-04` -> `006127a`).
48. Final-anchor closure reruns complete with explicit final-tag verification: `/tmp/spmd_canonical_20260304_0001/phase59_manifest_verify_finalanchor_20260304_020248` (`CLOSURE_MANIFEST_VERIFY_OK`) and `/tmp/spmd_canonical_20260304_0001/phase60_full_audit_finalanchor_20260304_020249` (`SPMD_CLOSURE_FULL_AUDIT_OK`).
49. Post-closeout maintenance tranche-1 completed: unreachable `.npRmpi_spmd_next_seq` removed; route-gated validation recorded in `/tmp/spmd_tranche1_deadscaffold_20260304_103024`.
50. Post-closeout maintenance tranche-2 completed: `mpi.comm.spawn` scalar positive-integer `nslaves` hardening plus new contract test (`test-rcomm-arg-contract.R`), with route-gated validation in `/tmp/spmd_tranche2_modsweep_20260304_103653`.
51. Post-closeout maintenance tranche-3 completed: modernization micro-scan harness added (`issue_notes/run_modernization_micro_scan.sh`, shellcheck-clean) with inventory artifact `/tmp/spmd_tranche3_modscan_20260304_104115`.
52. Post-closeout maintenance tranche-4 completed: modernization micro-scan harness refined to reduce false positives (token-boundary `eval/parse`, active-call-only `library/require`) with artifact `/tmp/spmd_tranche4_modscan_refine_20260304_104220`.
53. Post-closeout maintenance tranche-5 completed: stale commented `parse(...)` fragments removed from `R/Rcoll.R`, route-gated validation recorded in `/tmp/spmd_tranche5_gate_20260304_104321`, and refreshed micro-scan artifact `/tmp/spmd_tranche5_modscan_cleanup_20260304_104312` (`scan_parse=0`).
54. Post-closeout maintenance tranche-6 completed: IV-family type-stability hardening replaced classifier `sapply(..., is.numeric)` with typed `vapply(..., logical(1))` in `R/npregiv.R` and `R/npregivderiv.R`, with route-gated validation and refreshed modernization scan (`scan_sapply=67`) in dedicated tranche artifacts.
55. Post-closeout maintenance tranche-7 completed: `W.lp` classifier hardening in `R/util.R` (`xdat.col.numeric`: `sapply` -> typed `vapply(logical(1))`) with route-gated validation (`session`/`attach`/`profile`/`manual` pass tokens), installed-package `npreg`/`npsigtest` tests passing, and refreshed scan summary (`scan_sapply=66`).
56. Post-closeout maintenance tranche-8 completed: modernization scan harness I/O hardening in `issue_notes/run_modernization_micro_scan.sh` (explicit `ROOT_DIR`/`OUT_DIR` arg support and fail-fast root guard), validated by `shellcheck` and explicit out-dir scan smoke.
57. Post-closeout maintenance tranche-9 completed: plot-engine factor-classifier type-stability hardening in `R/np.plot.engine.{pl,sc,cond,con}bandwidth.R` (`all.isFactor`: `sapply`/`unlist(sapply)` -> typed `vapply(logical(1))`), with plot contract tests, session plot smoke, attach/profile/manual route validators, and refreshed scan summary (`scan_sapply=62`).
58. Post-closeout maintenance tranche-10 completed: installed-package test-hygiene hardening for `test-npreg-arg-contract.R` using namespace-safe lookup (`getFromNamespace(\"npreg.rbandwidth\", \"npRmpi\")`), with installed test pass and full route validator pass.

Latest checkpoint commits:
1. `8f5e959` (density/distribution opcode classification)
2. `47bbc49` (locked typed handlers for migrated core opcodes)
3. `388d643` (contract test realignment)
4. `9d163f6` (lock typed opcode for `npindexbw`)
5. `f4c83f0` (lock typed opcode for `npindex`)
6. `c34ff68` (lock typed opcodes for remaining core estimators)
7. `138401a` (fix forwarded-dot argument resolution in autodispatch)
8. `519955b` (lock non-core autodispatch families to typed SPMD opcodes)
9. `55226a0` (`npsigtest` locked opcode migration + contract and route coverage)
10. `13ddf0a` (phase-22 closeout hardening: session seq reset, fixed worker type packets, LC CVLS stable path, and `npsigtest` fast-fail pre-dispatch validation)
11. `3f307f7` (test contract hardening: initialize/cleanup MPI pool in `test-npsigtest-spmd-contract.R` for tarball-check reliability)
12. `a6ce628` (canonical plan checkpoint update for phase-22/23 closeout)
13. `0e161a4` (canonical status update marking remediation implementation complete and phase-23 baseline artifact root)
14. `965b44e` (append final canonical checkpoint hash for completed remediation record)
15. `cdde2a2` (phase-24 canonical closeout artifact update)
16. `ca987f7` (phase-25 terminal confirmation artifacts)
17. `c6a968c` (record local note-only as-cran closeout with compatibility shim evidence)
18. `7d571b7` (record phase-26 assurance rerun artifacts)
19. `4b120b4` (record phase-27 full assurance rerun artifacts)
20. `6b3b522` (broadcast-safe step-contract subprocess tests + phase-28 canonical closeout record)
21. `c2dce64` (record phase-29 terminal reconfirmation artifacts)
22. `795f34a` (add deterministic local as-cran harness script and record phase-30 closeout)
23. `635ddbb` (post-closeout maintenance tranche-1 dead-scaffold cull + remediation-doc refinement)
24. `bc6cbe9` (post-closeout maintenance tranche-2 `mpi.comm.spawn` contract hardening + tests)
25. post-closeout maintenance tranche-3 (modernization micro-scan harness + doc ledger update)
26. post-closeout maintenance tranche-4 (micro-scan false-positive reduction + doc ledger update)
27. post-closeout maintenance tranche-5 (Rcoll parse-comment cleanup + route/scans ledger update)
28. post-closeout maintenance tranche-6 (npregiv/npregivderiv classifier `vapply` hardening + route/scans ledger update)
29. post-closeout maintenance tranche-7 (`W.lp` classifier `vapply` hardening + route/scans ledger update)
30. post-closeout maintenance tranche-8 (modernization scan harness arg/IO hardening + validation ledger update)
31. post-closeout maintenance tranche-9 (plot-engine `all.isFactor` classifier `vapply` hardening + route/scans ledger update)
32. post-closeout maintenance tranche-10 (installed test namespace-lookup hardening for `npreg` arg-contract + route/scans ledger update)

Latest artifact roots:
1. `/tmp/spmd_canonical_20260304_0001/phase10_density_opcode_timeout_20260303_201934`
2. `/tmp/spmd_canonical_20260304_0001/phase11_locked_opcode_handlers_20260303_202742`
3. `/tmp/spmd_canonical_20260304_0001/phase12_individual_tests_20260303_204052`
4. `/tmp/spmd_canonical_20260304_0001/phase19_preflight_20260303_210733`
5. `/tmp/spmd_canonical_20260304_0001/phase20_noncore_opcode_locks_20260303_211050`
6. `/tmp/spmd_canonical_20260304_0001/phase21_npsigtest_trancheB_20260303_222056`
7. `/tmp/spmd_canonical_20260304_0001/phase22_closeout_20260303_223604`
8. `/tmp/spmd_canonical_20260304_0001/phase23_finalize_20260303_2330`
9. `/tmp/spmd_canonical_20260304_0001/phase24_finalclose_20260303_233544`
10. `/tmp/spmd_canonical_20260304_0001/phase25_terminal_20260303_233815`
11. `/tmp/spmd_canonical_20260304_0001/phase26_assurance_20260303_235554`
12. `/tmp/spmd_canonical_20260304_0001/phase27_fullsuite_20260304_000422`
13. `/tmp/spmd_canonical_20260304_0001/phase28_invariants_20260304_001450`
14. `/tmp/spmd_canonical_20260304_0001/phase29_terminal_20260304_003554`
15. `/tmp/spmd_canonical_20260304_0001/phase30_clean_harness_20260304_004703`
16. `/tmp/spmd_canonical_20260304_0001/phase31_manual_gatepack_20260304_005852`
17. `/tmp/spmd_canonical_20260304_0001/phase32_final_gatepack_20260304_010940`
18. `/tmp/spmd_canonical_20260304_0001/phase33_archival_gatepack_20260304_011716`
19. `/tmp/spmd_canonical_20260304_0001/FINAL_EVIDENCE_INDEX_20260304.md`
20. `/Users/jracine/Development/np-npRmpi/SPMD_REMEDIATION_FINAL_HANDOFF_2026-03-04.md`
21. `/Users/jracine/Development/np-npRmpi/SPMD_REMEDIATION_CLOSURE_CERTIFICATE_2026-03-04.md`
22. `/Users/jracine/Development/np-npRmpi/SPMD_REMEDIATION_CLOSURE_MANIFEST_2026-03-04.json`
23. `/Users/jracine/Development/np-npRmpi/issue_notes/verify_spmd_closure_manifest.sh`
24. `/tmp/spmd_canonical_20260304_0001/phase34_manifest_verify_20260304_013133`
25. `/Users/jracine/Development/np-npRmpi/issue_notes/export_spmd_closure_bundle.sh`
26. `/tmp/spmd_canonical_20260304_0001/phase35_closure_bundle_20260304_013242`
27. `/Users/jracine/Development/np-npRmpi/issue_notes/run_spmd_closure_full_audit.sh`
28. `/tmp/spmd_canonical_20260304_0001/phase36_full_audit_20260304_013355`
29. `/Users/jracine/Development/np-npRmpi/issue_notes/package_spmd_closure_archive.sh`
30. `/tmp/spmd_canonical_20260304_0001/phase37_closure_archive_20260304_013513`
31. `/tmp/spmd_canonical_20260304_0001/phase38_full_audit_20260304_013613`
32. `/Users/jracine/Development/np-npRmpi/SPMD_REMEDIATION_CLOSEOUT_WORKFLOW_2026-03-04.md`
33. `/tmp/spmd_canonical_20260304_0001/phase40_manifest_verify_20260304_013945`
34. `/tmp/spmd_canonical_20260304_0001/phase41_full_audit_20260304_013952`
35. `/tmp/spmd_canonical_20260304_0001/phase47_stagewise_gatepack_20260304_014609`
36. `/tmp/spmd_canonical_20260304_0001/phase48_manifest_verify_20260304_015705`
37. `/tmp/spmd_canonical_20260304_0001/phase49_full_audit_20260304_015705`
38. `/tmp/spmd_canonical_20260304_0001/phase53_manifest_verify_dynamic_20260304_015851`
39. `/tmp/spmd_canonical_20260304_0001/phase54_full_audit_dynamic_20260304_015854`
40. `/tmp/spmd_canonical_20260304_0001/phase57_manifest_verify_tagged_20260304_020048`
41. `/tmp/spmd_canonical_20260304_0001/phase58_full_audit_tagged_20260304_020052`
42. `/tmp/spmd_canonical_20260304_0001/phase59_manifest_verify_finalanchor_20260304_020248`
43. `/tmp/spmd_canonical_20260304_0001/phase60_full_audit_finalanchor_20260304_020249`
44. `/tmp/spmd_tranche1_deadscaffold_20260304_103024`
45. `/tmp/spmd_tranche2_modsweep_20260304_103653`
46. `/tmp/spmd_tranche3_modscan_20260304_104115`
47. `/tmp/spmd_tranche4_modscan_refine_20260304_104220`
48. `/tmp/spmd_tranche5_modscan_cleanup_20260304_104312`
49. `/tmp/spmd_tranche5_gate_20260304_104321`
50. `/tmp/spmd_tranche6_iv_vapply_20260304_104740`
51. `/tmp/spmd_tranche6_modscan_postvapply_20260304_104841`
52. `/tmp/spmd_tranche7_wlp_vapply_20260304_105424`
53. `/tmp/spmd_tranche8_scanio_20260304_105725`
54. `/tmp/spmd_tranche9_plot_vapply_20260304_110304`
55. `/tmp/spmd_tranche10_npreg_arg_contract_20260304_111022`

## Objective
Keep user-facing workflow unchanged (`npreg(...)`, `npregbw(...)`, etc.) while making internal execution rank-symmetric SPMD for MPI-sensitive paths in all modes:
1. `session/spawn`
2. `attach`
3. `profile/manual-broadcast`

Core family coverage for this plan:
1. `npreg*`
2. `npudens*`
3. `npcdens*`
4. `npudist*`
5. `npcdist*`
6. `npscoef*`
7. `npindex*`
8. `npplreg*`

Important scope clarification:
1. The core-family list is a migration priority list, not an allowlist.
2. All exported `npRmpi` functions (including non-core families such as `npsigtest`, `npsymtest`, `npcopula`) are expected to run.
3. Non-core families are currently supported under the same runtime contract and are migrated to typed locked opcodes in Phase 6.

## Non-Negotiable Invariants
1. No silent serial fallback when MPI route is selected.
2. No runtime `np::` bridge from `npRmpi`.
3. No masked failures in demos/tests.
4. No default drift unless explicitly approved.
5. One risk axis per tranche.
6. No `O(n^2)` helper-path reintroduction in estimation/CV internals.
7. Helper policy: prefer direct distributed kernels in MPI-sensitive CV/estimation hot paths.

## Baseline and Artifacts
1. Baseline commit for active branch: `45f0e2cbe81dcb4fb40e4d8b01ce8c54e23eeaf6`.
2. All tranche artifacts go under:
   - `/tmp/spmd_canonical_<timestamp>/phase0`
   - `/tmp/spmd_canonical_<timestamp>/phase1`
   - `/tmp/spmd_canonical_<timestamp>/phase2`
   - `/tmp/spmd_canonical_<timestamp>/phase3`
   - `/tmp/spmd_canonical_<timestamp>/phase4`
   - `/tmp/spmd_canonical_<timestamp>/phase5`
   - `/tmp/spmd_canonical_<timestamp>/phase6`

## Required Pre-Flight (Before Every Tranche)
1. `git status --short` for:
   - `/Users/jracine/Development`
   - `/Users/jracine/Development/np-master`
   - `/Users/jracine/Development/np-npRmpi`
2. Demo masking scan:
   - `rg -n --glob 'demo/*.R' '(^|[^#])\\btry\\s*\\(|(^|[^#])\\btryCatch\\s*\\('`
3. Orphan worker scan:
   - `ps aux | rg 'slavedaemon\\.R|Rslaves\\.sh'`
4. Route sanity (`npRmpi`, tiny smokes):
   - session/spawn (`nslaves=1`)
   - attach (`mpiexec -n 2`)

## Route/Smoke Commands (Existing Tools Only)
1. Manual-broadcast validator:
   - `Rscript --no-save /Users/jracine/Development/np-npRmpi/issue_notes/validate_route_manual_broadcast.R`
2. Attach validator:
   - `FI_TCP_IFACE=en0 mpiexec -n 2 Rscript --no-save /Users/jracine/Development/np-npRmpi/issue_notes/validate_route_attach.R`
3. Profile validator:
   - `R_PROFILE_USER=/Users/jracine/Development/np-npRmpi/inst/Rprofile FI_TCP_IFACE=en0 mpiexec -n 2 Rscript --no-save /Users/jracine/Development/np-npRmpi/issue_notes/validate_route_profile.R`
4. Session tiny smoke (`nslaves=1`) token check:
   - `Rscript --no-save -e 'suppressPackageStartupMessages(library(npRmpi)); npRmpi.init(nslaves=1, quiet=TRUE); on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE); set.seed(1); x<-runif(80); y<-rnorm(80); bw<-npregbw(y~x, regtype=\"lc\", bwmethod=\"cv.ls\", nmulti=1); fit<-npreg(bws=bw); stopifnot(inherits(fit,\"npregression\")); cat(\"SESSION_N1_OK\\n\")'`

## Target Internal Contract (SPMD Step Engine)
For migrated MPI-sensitive opcodes:
1. Master creates envelope: `{seq_id, opcode, args_ref, timeout_class}`.
2. Master broadcasts envelope once per step.
3. All ranks enter same opcode handler.
4. All ranks execute matching collective cadence for that opcode.
5. Completion uses explicit ACK/ERR semantics.
6. Sequence mismatch or timeout is hard-fail with route/opcode diagnostics.

## Phase Plan

## Phase 0: Baseline Freeze
Scope:
1. Snapshot current clean baseline and execution environment.
2. Pin reference logs for route, parity, and performance.

Validation:
1. Route validators all pass.
2. Targeted MPI-sensitive tests pass:
   - `tests/testthat/test-session-routing-subprocess-contract.R`
   - `tests/testthat/test-session-arg-contract.R`
   - `tests/testthat/test-jksum-gating-smoke.R`
   - `tests/testthat/test-ll-lp-degree1-parity.R`

Acceptance:
1. `/tmp` phase0 bundle includes commands, raw logs, and summary table.

## Phase 1: Remove `nslaves=0` Permanently
Scope:
1. Enforce `nslaves>=1` in runtime lifecycle.
2. Remove remaining active `nslaves==0` branches in runtime paths.
3. Update docs/tests/messages to reflect MPI-only contract.

Validation:
1. Negative contract test:
   - `npRmpi.init(nslaves=0)` hard-fails with explicit remediation (`use np for serial`).
2. Positive route check:
   - session smoke with `nslaves=1` passes.
3. Attach/profile/manual validators pass.
4. No orphan worker processes.

Acceptance:
1. No active runtime path supports `nslaves=0`.
2. User-facing message contract is intentional and tested.

## Phase 2: Control-Plane Scaffolding (No Method Logic Change)
Scope:
1. Add envelope structure and `seq_id` management.
2. Add typed opcode registry + dispatcher shell.
3. Add ACK/ERR completion scaffolding.
4. Add timeout/divergence diagnostics.

Validation:
1. One tiny test opcode works in session/attach/profile.
2. Deliberate sequence mismatch fails fast with explicit diagnostics.
3. Existing validators and targeted tests still pass.

Acceptance:
1. Infrastructure exists without altering estimator math.
2. No route regressions.

## Phase 3: Migrate `npregbw` LL/LP CVLS/CVAIC First
Scope:
1. Move highest-risk CV route to typed opcode path.
2. Ensure strict collective-order symmetry in session path.
3. Remove legacy fallback for migrated calls only.
4. Excise known `O(n^2)` helper-path usage from migrated LL/LP CVLS routes.

Validation:
1. Historical session hang repro no longer hangs.
2. Attach/profile behavior preserved.
3. Numerical parity:
   - objective value parity (`fval`) where applicable
   - `num.feval` parity/expected behavior
4. Decision-tier performance gate:
   - interleaved paired `times=25` screening
   - escalate to `times>=100` when not obvious
5. Static helper-excision gate on migrated paths:
   - `rg -n 'np_reg_cv_ls_stable_ll_glp|stable_ll_glp|O\\(n\\^2\\)' /Users/jracine/Development/np-npRmpi/src /Users/jracine/Development/np-npRmpi/R`
   - only acceptable matches are comments/archival notes, not active CV route logic.

Acceptance:
1. Session route becomes stable for LL/LP CVLS/CVAIC path.
2. No numerical regression beyond declared tolerance.
3. No active `O(n^2)` helper dependency remains in the migrated LL/LP route.

## Phase 4: Migrate Remaining Core Families
Scope:
1. Migrate sensitive paths in:
   - `npudens*`
   - `npudist*`
   - `npcdens*`
   - `npcdist*`
   - `npscoef*`
   - `npindex*`
   - `npplreg*`
2. Keep migration opcode-by-opcode.

Validation:
1. Route matrix pass for each migrated family.
2. Family-level parity/performance evidence retained per tranche.

Acceptance:
1. No migrated hot path relies on legacy asymmetric dispatch.

## Phase 5: Decommission Legacy Asymmetry in Hot Paths
Scope:
1. Remove free-form worker eval for estimator/CV hot paths.
2. Keep only controlled opcode execution for migrated families.
3. Finalize docs and route guidance.

Validation:
1. Full route matrix pass:
   - `np` serial anchor (`np-master`)
   - `npRmpi` session/spawn (`nslaves=1` minimum + `nslaves>1` where relevant)
   - `npRmpi` attach (`mpiexec -n 2`)
   - `npRmpi` profile/manual-broadcast (`mpiexec -n 2`)
2. Demo triplet fail-fast parity checks pass.

Acceptance:
1. Session is operationally SPMD-equivalent to attach/profile for migrated hot paths.

## Phase 6: Extend to Remaining Non-Core Function Families
Scope:
1. Migrate additional MPI-sensitive non-core/orchestrator/helper routes not in the core family set.
2. Preserve orchestrator boundary policy:
   - orchestrators remain local unless whole-wrapper dispatch is explicitly route-validated.
3. Keep one-risk-axis micro-tranches and retain route-safe behavior in attach/profile/session.

Validation:
1. Route matrix pass for each touched non-core family.
2. Numerical parity and condition-contract stability for each touched API.
3. Performance gate applied where touched path is hot/iterative.

Acceptance:
1. Touched non-core routes are SPMD-safe under session/attach/profile without masked fallback behavior.

## Gate Pack Required At Every Phase
1. Route gate:
   - session, attach, profile/manual validators.
2. Numerical gate:
   - objective and bandwidth parity (or bounded/justified drift).
3. Performance gate:
   - paired interleaved seeds,
   - declared MEI (absolute + relative),
   - mean + median + tail + spread checks,
   - multiplicity-adjusted beneficiary superiority + unaffected equivalence.
4. Demo/test gate:
   - demo masking scan clean,
   - targeted tests pass,
   - subprocess timeout guards used where appropriate.
5. Cleanup gate:
   - no orphan `slavedaemon.R`/`Rslaves.sh`.
6. Complexity gate:
   - no active introduction of `O(n^2)` helper storage/paths in estimation/CV internals.

## Commit and Rollback Policy
1. Commit only after full phase gates pass.
2. Commit message must include:
   - changed scope,
   - validated routes/methods,
   - artifact path.
3. Immediate rollback of tranche if any hard gate fails:
   - hang/deadlock,
   - numerical drift beyond tolerance,
   - performance regression beyond MEI,
   - condition-message regression not explicitly accepted.

## Immediate Next Action
1. Keep route gate (`session`, `attach`, `profile/manual`) and orphan-cleanup checks mandatory for any post-remediation patch.
2. Use `/tmp/spmd_canonical_20260304_0001/phase23_finalize_20260303_2330` as the post-remediation baseline for any future tranche.

## Deferred Backlog (Not In Current Tranche)
1. Enable `coef=TRUE` smooth-coefficient plot error bands in `np` and `npRmpi` without method remapping.
2. Phase A (low risk): support asymptotic coefficient bands by wiring `gerr[, coef.index]` into `np.plot.engine.scbandwidth.R` (both repos), with fail-fast for unsupported combinations.
3. Phase B (separate tranche): implement coefficient-specific bootstrap targets for `scbandwidth` plots (current bootstrap helper centers on mean-path targets).
4. Keep the current `coef` + plot-errors disable guard until Phase A tests pass across `np` serial and `npRmpi` session/attach/profile routes.
