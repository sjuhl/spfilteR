# spfilteR 2.2.0

- new selection criteria for unsupervised eigenvector selection via stepwise regression:
  - selection based on the corrected Akaike information criterion ('AICc') now available in `glmFilter()`
  - `lmFilter()` now supports eigenvector selection based on AIC, AICc, and BIC
- lasso-based eigenvector selection now supported in `lmFilter()`
- `lmFilter()` now allows to compute conditional standard errors for regression coefficients using a partial regression framework.
- minor adjustment to the summary method
- add deviance residuals for negative binomial regression models
- add warning message that deviance residuals will become the default in `glmFilter()` in future releases
- update vignette & tests

---

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

---

# spfilteR 2.0.0

- allow for unsupervised eigenvector selection in negative binomial models
    - `glmFilter()` now supports 'nb' (for negative binomial) as model type
    - adjustments in summary method and helper functions to handle negative binomial models
    - update tests for negative binomial model
- `glmFilter()` also provides McFadden's adjusted pseudo R-squared for the filtered vs. the unfiltered model
- bug fixes
    - bug fix in help pages of `MI.local()` and `MI.vec()` functions
    - resolves an error in `lmFilter()` and `glmFilter()` occurring when covariates are supplied as data.frame
- assign variable names to output (if provided)
- improve the handling of missingness in `MI.vec()`, `MI.decomp()`, and `MI.local()`
- update tests

---

# spfilteR 1.1.5

- fix minor bug in help pages

---

# spfilteR 1.1.4

- update citation information
- fix: use isTRUE(all.equal()) instead of "==" on numeric vectors

---

# spfilteR 1.1.3

- update tests

---

# spfilteR 1.1.2

- fix broken links
- update citation information

---

# spfilteR 1.1.1

- CRAN resubmission
- improve readability of code
- update author mail address
- include citation

---

# spfilteR 1.1.0.9000

- improve readability of code
- update author mail address
- include citation

---

# spfilteR 1.1.0

- CRAN resubmission
- fix minor bug when checking 'tol' in `lmFilter()` function
- new vignette name
- add new functions
    - `MI.local()` function to calculate local Moran's I
    - `vp()` function for variation partitioning
- update reference to `MI.local()` in documentation files

---

# spfilteR 1.0.0.9000

- fix minor bug when checking 'tol' in `lmFilter()` function
- new vignette name
- add new functions
    - `MI.local()` function to calculate local Moran's I
    - `vp()` function for variation partitioning
- update reference to `MI.local()` in documentation files

---

# spfilteR 1.0.0

- include `MI.decomp()` to decompose Moran's I
- rename `MI.resid()`
- prepare for CRAN submission

---

# spfilteR 0.2.0

- `lmFilter()` and `glmFilter()` now also support unsupervised eigenvector selection based on the significance of residual autocorrelation

---

# spfilteR 0.1.0
