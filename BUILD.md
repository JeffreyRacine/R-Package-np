# Build (np)

## Build Tarball

```bash
cd /Users/jracine/Development
R CMD build np-master
```

This produces `np_0.70-0.tar.gz` in `/Users/jracine/Development`.

## Install

```bash
cd /Users/jracine/Development
R CMD INSTALL np_0.70-0.tar.gz
```

## Quick Load Check

```bash
R -q -e 'library(np); sessionInfo()'
```

## Check

```bash
cd /Users/jracine/Development
R CMD check --as-cran np_0.70-0.tar.gz
```
