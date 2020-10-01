# spfilteR

This package provides a number of useful functions that facilitate the analysis of spatially autocorrelated data based on the eigenfunction decomposition of an exogenously given connectivity matrix ***W***. The main function `getEVs` specifies a projection matrix ***M*** and symmetrizes the connectivity matrix by ***V***=1/2*(***W***+***W***') before decomposing the transformed matrix ***MVM***. If covariates are supplied, this function constructs the projection matrix by: ***M***=***I***-***X***(***X***'***X***)^-1 ***X***'. Eigenvectors obtained from this specification are not only mutually uncorrelated but also orthogonal to the covariates in ***X***. In contrast, if no covariates are supplied, the projection matrix simplifies to ***M***= ***I***-***11***'/*n*, where ***1*** is a vector of ones and *n* is the number of units. Besides the eigenvectors and eigenvalues, `getEVs` also provides the Moran coefficient associated with each eigenvector.

Subsequently, these eigenvectors can be used to perform semiparametric spatial filtering in a regression framework.
