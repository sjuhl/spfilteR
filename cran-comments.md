# spfilteR 2.1.0

- update vignette to include an example of the negative binomial model
- imrprovement in the console output of function `vp()`
- improvements and bug fixes in `MI.vec()` and `MI.decomp()`:
    - correctly handle missing values in each variable separately if multiple variables are supplied
    - check for variable names only inside the function and not in the global environment
    - removal of constant terms supplied to the functions
- improve the handling of missingness in `MI.resid()`
- minor adjustments to helper functions
- update tests


# Test Results:
── R CMD check results ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── spfilteR 2.1.0 ────
Duration: 54s

❯ checking for future file timestamps ... NOTE
  unable to verify current time

0 errors ✔ | 0 warnings ✔ | 1 note ✖
