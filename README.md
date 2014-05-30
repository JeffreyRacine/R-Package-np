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

Note Windows users have to first install
[Rtools](http://cran.r-project.org/bin/windows/Rtools), while OS X
users have to first install
[Xcode](https://itunes.apple.com/us/app/xcode/id497799835) and the
command line tools (in OS X 10.9 or higher, once you have Xcode
installed, open a terminal and run xcode-select --install).

Note that if you wish to install the MPI-enabled development version
of the package (i.e. the package `npRmpi'), you can add the option
ref='npRmpi' to the install_github call above presuming that your
system has the required MPI subsystem installed (see my homepage for
further details).

If you wish a fast install without the building of vignettes, add the
option build_vignettes=FALSE to the install_github call.

For more information on this project please visit the maintainer's website (http://www.economics.mcmaster.ca/faculty/racinej).

