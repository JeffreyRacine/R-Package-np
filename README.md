# npRmpi

This is the R package `npRmpi' (NonpRmpiarametric Kernel Methods for Mixed Datatypes) written and maintained by Jeffrey S. Racine (racinej@mcmaster.ca) and co-authored by Tristen Hayfield <hayfield@mpia.de>

## Installation

Presuming that a working implementation of MPI exists on the target machine, you can install the stable version on [CRAN](http://cran.r-project.org/package=npRmpi):

```r
install.packages('npRmpi', dependencies = TRUE)
```

Or download the [zip ball](https://github.com/JeffreyRacine/R-Package-npRmpi/zipball/master) or [tar ball](https://github.com/JeffreyRacine/R-Package-npRmpi/tarball/master), decompress and run `R CMD INSTALL` on it, or install then use the **devtools** package to install the development version:

```r
library(devtools); install_github('R-Package-npRmpi', 'JeffreyRacine')
```

Note Windows users have to first install [Rtools](http://www.murdoch-sutherland.com/Rtools/).

For more information on this project please visit the maintainer's website (http://www.economics.mcmaster.ca/faculty/racinej).

