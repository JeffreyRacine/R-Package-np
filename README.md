# npRmpi

This is the R package `npRmpi' (NonpRmpiarametric Kernel Methods for Mixed Datatypes) written and maintained by Jeffrey S. Racine (racinej@mcmaster.ca) and co-authored by Tristen Hayfield (tristen.hayfield@gmail.com)

## Installation

Presuming that a working implementation of MPI exists on the target machine, you can install the stable version on [CRAN](http://cran.r-project.org/package=npRmpi):

```r
install.packages('npRmpi', dependencies = TRUE)
```

Or download the [zip ball](https://github.com/JeffreyRacine/R-Package-npRmpi/zipball/master) or [tar ball](https://github.com/JeffreyRacine/R-Package-npRmpi/tarball/master), decompress and run `R CMD INSTALL` on it, or install then use the **devtools** package to install the development version:

```r
library(devtools); install_github('R-Package-npRmpi', 'JeffreyRacine')
```

Note also that if you wish a fast install without the building of
vignettes (or if you do not have TeX installed on your system), add
the option build_vignettes=FALSE to the install_github() call.

Note that if you wish to install the MPI-enabled development version
of the package (i.e. the package `npRmpi'), you can add the option
ref='npRmpi' to the install_github call above presuming that your
system has the required MPI subsystem installed (see my homepage for
further details).

For more information on this project please visit the maintainer's website (http://www.economics.mcmaster.ca/faculty/racinej).

