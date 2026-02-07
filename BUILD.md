# Build (np)

## Build Tarball

```bash
R CMD build .
```

This produces `np_0.60-20.tar.gz` in the same directory.

## Install

```bash
R CMD INSTALL np_0.60-20.tar.gz
```

## Quick Load Check

```bash
R -q -e 'library(np); sessionInfo()'
```
