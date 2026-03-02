## R CMD check results

Local (macOS, R 4.5.2, --as-cran):

0 errors | 1 warning | 2 notes

### Notes

1. `checking CRAN incoming feasibility ... NOTE`
- Version jump from CRAN `0.60-20` to submission `0.70-1` is intentional for a major modernization release.

2. `checking for future file timestamps ... NOTE`
- Environment clock verification note (`unable to verify current time`) observed intermittently on local host.

### Warning

1. `checking top-level files ... WARNING`
- Non-standard top-level directory `archive/` retained intentionally for project provenance/history.
