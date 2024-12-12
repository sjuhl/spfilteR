# spfilteR 2.0

# NEWS.md:
- allow for unsupervised eigenvector selection in negative binomial models
    - `glmFilter()` now supports 'nb' (for negative binomial) as model type
    - adjustments in summary method and helper functions to handle negative binomial models
    - update tests for negative binomial model
- `glmFilter()` also provides McFadden's adjusted pseudo R-squared for the filtered vs. the unfiltered model
- bug fixes
    - bug fix in help pages of MI.local() and MI.vec() functions
    - resolves an error in `lmFilter()` and `glmFilter()` occurring when covariates are supplied as data.frame


# Test Results:
── R CMD check results ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── spfilteR 2.0 ────
Duration: 49s

❯ checking for future file timestamps ... NOTE
  unable to verify current time

0 errors ✔ | 0 warnings ✔ | 1 note ✖
