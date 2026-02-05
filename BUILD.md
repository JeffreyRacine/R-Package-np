# Build (np)

## Build Tarball

```bash
cd /Users/jracine/Development/np-master
R CMD build np
```

This produces `np_0.60-21.tar.gz` in the same directory.

## Install

```bash
cd /Users/jracine/Development/np-master
R CMD INSTALL np_0.60-21.tar.gz
```

## Quick Load Check

```bash
R -q -e 'library(np); sessionInfo()'
```
