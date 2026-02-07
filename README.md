# np

![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/np)

This is the R package `np` (Nonparametric Kernel Methods for Mixed Datatypes) written and maintained by Jeffrey S. Racine (racinej@mcmaster.ca) and co-authored by Tristen Hayfield (tristen.hayfield@gmail.com).

## Installation

You can install the stable version on [CRAN](https://cran.r-project.org/package=np):

```r
install.packages('np', dependencies = TRUE)
```

Or download the [zip ball](https://github.com/JeffreyRacine/R-Package-np/zipball/master)
or [tar ball](https://github.com/JeffreyRacine/R-Package-np/tarball/master),
decompress and run `R CMD INSTALL` on it.

Alternatively, you can install the development version but before
doing so Windows users have to first install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/), while macOS
users have to first install
[Xcode](https://apps.apple.com/us/app/xcode/id497799835?mt=12) and the
command line tools (once you have Xcode
installed, open a terminal and run `xcode-select --install`). Note also
that versions of e.g. Rtools are paired with versions of R so ensure
you have the latest version of R installed prior to commencing this
process.

After installing
[Rtools](https://cran.r-project.org/bin/windows/Rtools/)/[Xcode](https://apps.apple.com/us/app/xcode/id497799835?mt=12)
and **devtools** (via install.packages("devtools")), install the
development package using the following command:

```r
library(devtools); install_github('JeffreyRacine/R-Package-np')
```

Note also that if you wish a fast install without the building of
vignettes (or if you do not have TeX installed on your system), add
the option build_vignettes=FALSE to the install_github() call.

Note that if you wish to install the MPI-enabled development version
of the package (i.e. the package `npRmpi`), you can add the option
`ref='npRmpi'` to the `install_github()` call above, presuming that your
system has a working MPI subsystem installed. See `WORKTREES.md` and
`BUILD.md` in this repo for local build details.

For more information on this project please visit the maintainer's website (https://experts.mcmaster.ca/people/racinej).
