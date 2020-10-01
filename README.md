# spfilteR

This package provides a number of useful functions that facilitate the analysis of spatially autocorrelated data based on the eigenfunction decomposition of an exogenously given connectivity matrix ***W***. The main function 'getEVs' specifies a projection matrix ***M*** and symmetrizes the connectivity matrix by ***V***=0.5*(***W***+***W***') before decomposing the transformed matrix ***MVM***.
