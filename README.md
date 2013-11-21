# np

This is the R package `np' (Nonparametric Kernel Methods for Mixed Datatypes) written and maintained by Jeffrey S. Racine (racinej@mcmaster.ca) and co-authored by Tristen Hayfield (hayfield@mpia.de)

## Installation

You can install the stable version on [CRAN](http://cran.r-project.org/package=np):

```r
install.packages('np', dependencies = TRUE)
```

Or download the [zip ball](https://github.com/JeffreyRacine/R-Package-np/zipball/master) or [tar ball](https://github.com/JeffreyRacine/R-Package-np/tarball/master), decompress and run `R CMD INSTALL` on it, or install then use the **devtools** package to install the development version:

```r
library(devtools); install_github('R-Package-np', 'JeffreyRacine')
```

Note Windows users have to first install
[Rtools](http://cran.r-project.org/bin/windows/Rtools), while OS X
users have to first install
[Xcode](https://itunes.apple.com/us/app/xcode/id497799835) and the
command line tools (in OS X 10.9 or higher, once you have Xcode
installed, open a terminal and run xcode-select --install).

For more information on this project please visit the maintainer's website (http://www.economics.mcmaster.ca/faculty/racinej).

