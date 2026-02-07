# npRmpi

This is the R package `npRmpi` (Parallel Nonparametric Kernel Methods for Mixed Datatypes) written and maintained by Jeffrey S. Racine (racinej@mcmaster.ca) and co-authored by Tristen Hayfield (tristen.hayfield@gmail.com).

## Installation

Presuming that a working implementation of MPI exists on the target machine, you can install the stable version on [CRAN](https://cran.r-project.org/package=npRmpi):

```r
install.packages('npRmpi', dependencies = TRUE)
```

Or download the [zip ball](https://github.com/JeffreyRacine/R-Package-np/zipball/npRmpi) or
[tar ball](https://github.com/JeffreyRacine/R-Package-np/tarball/npRmpi),
decompress and run `R CMD INSTALL` on it, or use the **devtools** package
to install the development version:

```r
library(devtools); install_github('JeffreyRacine/R-Package-np', ref = 'npRmpi')
```

Note also that if you wish a fast install without the building of
vignettes (or if you do not have TeX installed on your system), add
the option build_vignettes=FALSE to the install_github() call.

## MPI (MPICH via MacPorts) Quick Setup

```bash
export RMPI_TYPE=MPICH
export RMPI_INCLUDE=/opt/local/include/mpich-mp
export RMPI_LIB_PATH=/opt/local/lib/mpich-mp
export RMPI_LIBS="-L/opt/local/lib/mpich-mp -lmpi"
export CC=mpicc
export CXX=mpicxx
```

Then build/install from the repo:

```bash
R CMD build .
R CMD INSTALL npRmpi_0.60-20.tar.gz
```

See `BUILD.md` and `WORKTREES.md` in this repo for local build details.

For more information on this project please visit the maintainer's website (https://experts.mcmaster.ca/people/racinej).
