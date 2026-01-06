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


# Test Results:
── R CMD check results ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── spfilteR 2.2.0 ────
Duration: 55.5s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔
