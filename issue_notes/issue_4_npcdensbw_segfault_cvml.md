# Issue #4: npcdensbw cv.ls then cv.ml segfault (factor y)

**Type:** Bug (crash)

## Report
Calling `npcdensbw(..., bwmethod="cv.ls")` followed by `npcdensbw(..., bwmethod="cv.ml")` with factor `y` can crash the R session (segfault). Reproduction from issue:

```
bw <- npcdensbw(low ~ lwt, bwmethod="cv.ls", data=birthwt, nmulti=1)
bw <- npcdensbw(low ~ lwt, bwmethod="cv.ml", data=birthwt, nmulti=1)
```

## Local reproduction
Running the same script in this repo triggers a hard crash (signal 6) when the second call is executed.

## Likely cause
The crash occurs inside C routine `np_density_conditional_bw`. The second call may reuse global/static state that is not fully reinitialized when switching methods (cv.ls → cv.ml) with factor response data.

## Proposed fix (not yet implemented)
- Investigate C globals in `src/np.c` for `np_density_conditional_bw` and ensure per-call state is initialized and freed regardless of bwmethod.
- Consider explicit reset of any static buffers or global pointers used by the conditional density bandwidth routine.
- If no immediate fix, add a guard in R to prevent switching bwmethod without a new session or explicit reset (fail-fast with message).

## Risk assessment
High. Requires C-level debugging and careful memory/state audit. Not pursued yet.

## Status
Reproduced; no patch proposed yet.
