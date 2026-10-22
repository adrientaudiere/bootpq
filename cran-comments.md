## First submission

This is the first submission of bootpq to CRAN.

## R CMD check results

0 errors | 0 warnings | 1 note

* "New submission", and "Version contains large components (0.1.0.9000)"
  on the development version; the latter disappears with the release.

## Examples

Each function repeats an analysis over several rarefactions, so the examples
use few permutations (`nperm = 2`) and are wrapped in `\donttest{}`.

## Test environments

* local Linux (Pop!_OS 24.04), R 4.6.1, `R CMD check --as-cran --run-donttest`
  against the CRAN version of 'MiscMetabar' (0.16.8), 2026-10-22.
* win-builder R-devel (Windows Server 2022, 2026-10-22): 1 note (new
  submission; possibly misspelled technical words in DESCRIPTION).
