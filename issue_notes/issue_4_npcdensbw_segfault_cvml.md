# Issue #4 (alt): npcdensbw cv.ml segfault after cv.ls

**Status:** Resolved (see issue_4_npcdensbw_cvls_cvml_segfault.md)

This is the same underlying defect. Fix is in `src/np.c` (reset externals/pointers for conditional density between calls).
