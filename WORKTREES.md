# Worktrees And Pushing

## Current Layout

This repo is used with two worktrees:
- `/Users/jracine/Development/np-master` for the serial `np` package (`master` branch).
- `/Users/jracine/Development/np-npRmpi` for the MPI `npRmpi` package (`npRmpi` branch).

## Push `master` To GitHub

From the `np-master` worktree:

```bash
cd /Users/jracine/Development/np-master
git status -sb
git push origin master
```

## Push `npRmpi` To GitHub

From the `np-npRmpi` worktree:

```bash
cd /Users/jracine/Development/np-npRmpi
git status -sb
git push origin npRmpi
```

## Should You Keep Two Folders?

Yes, it is safer and more convenient to keep two worktrees:
- You can build/test both branches without switching or stashing.
- Each worktree stays clean and focused on its package identity.
If you prefer a single folder, you can remove one worktree and switch branches in place, but it is easier to make mistakes when switching.

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
