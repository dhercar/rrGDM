## Test environments

-   GitHub Actions (R-CMD-check): macOS-latest (release), windows-latest (release), ubuntu-latest (devel, release, oldrel-1)
-   Local: Windows 11, R 4.2.2

## R CMD check results

0 errors \| 0 warnings \| 1 note

-   This is a new submission.

* On Windows: "Compilation used the following non-portable flag(s):
  '-Wa,-mbig-obj'". The flag is included because without it the assembler
  fails with "too many sections ... file too big". It is set only in
  src/Makevars.win.

## Notes on the package

* The package contains compiled code using 'TMB' and 'RcppEigen'
  (LinkingTo).
