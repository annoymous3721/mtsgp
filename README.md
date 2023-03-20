# mtsgp
Implementation in R

# Dependencies
`install.packages(c("plgp", "mvtnorm", "dplyr","MASS","pracma","Matrix","condmixt","ggplot2"))`

# Rscripts
* Figure 2: `Rscript Figure2.R`
* Figure 3: Compare prediction MSE with SGP, MAGMA, splines
`Rscript Figure3.R`
* Figure 4: Derive MTSGP results `Rscript Figure4.R`. Save the datasets from R to compare ranking results with soft-dtw and soft-dtw divergence `Python dtw_ranking_Figure4.py`
* Table 1: Compare AUC with TLCC, splines, SGPs. `Rscript table1.R`
* Real data tables and plots: Derive MTSGP results `Rscript epCountsLoop_rbf.R`. Analyze the results with statistical tests `Rscript plot_epCountsLoop.R`. Data is too large to upload to Github, but can be requested if needed. 