# SPMD Remediation Closeout Workflow (2026-03-04)

Use this to re-validate closure quickly on this machine.

## Preconditions
1. Repo: `/Users/jracine/Development/np-npRmpi`
2. Base artifacts root exists: `/tmp/spmd_canonical_20260304_0001`

## Commands
1. Full audit rerun:
```bash
cd /Users/jracine/Development/np-npRmpi
./issue_notes/run_spmd_closure_full_audit.sh /tmp/spmd_canonical_20260304_0001/phaseXX_full_audit_$(date +%Y%m%d_%H%M%S)
```
Expected token: `SPMD_CLOSURE_FULL_AUDIT_OK`.

2. Bundle export:
```bash
cd /Users/jracine/Development/np-npRmpi
./issue_notes/export_spmd_closure_bundle.sh /tmp/spmd_canonical_20260304_0001/phaseXX_closure_bundle_$(date +%Y%m%d_%H%M%S)
```
Expected token: `CLOSURE_BUNDLE_EXPORT_OK`.

3. Archive package:
```bash
cd /Users/jracine/Development/np-npRmpi
./issue_notes/package_spmd_closure_archive.sh /tmp/spmd_canonical_20260304_0001/phaseXX_closure_archive_$(date +%Y%m%d_%H%M%S)
```
Expected token: `CLOSURE_ARCHIVE_PACKAGE_OK`.

## Canonical token chain
1. `FINAL_GATE_PACK_OK`
2. `CLOSURE_MANIFEST_VERIFY_OK`
3. `CLOSURE_BUNDLE_EXPORT_OK`
4. `SPMD_CLOSURE_FULL_AUDIT_OK`
5. `CLOSURE_ARCHIVE_PACKAGE_OK`

## Anchors
1. Canonical plan: `/Users/jracine/Development/np-npRmpi/REMEDIATION_PLAN_SPMD_CANONICAL_2026-03-03.md`
2. Final handoff: `/Users/jracine/Development/np-npRmpi/SPMD_REMEDIATION_FINAL_HANDOFF_2026-03-04.md`
3. Certificate: `/Users/jracine/Development/np-npRmpi/SPMD_REMEDIATION_CLOSURE_CERTIFICATE_2026-03-04.md`
4. Manifest: `/Users/jracine/Development/np-npRmpi/SPMD_REMEDIATION_CLOSURE_MANIFEST_2026-03-04.json`
5. Evidence index: `/tmp/spmd_canonical_20260304_0001/FINAL_EVIDENCE_INDEX_20260304.md`
