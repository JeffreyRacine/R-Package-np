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
