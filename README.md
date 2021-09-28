# bayesmove_method_comparison
[![DOI](https://zenodo.org/badge/289030089.svg)](https://zenodo.org/badge/latestdoi/289030089)


R code associated with Cullen et al. Methods in Ecology and Evolution paper *Identifying latent behavioral states in animal movement with M4, a non-parametric Bayesian method*.

This repository stores the R scripts used to generate the simulated data, analyze the data with each of the five methods that were compared (BCPA, EMbC, HMM, segclust2d, M4), and generate figures. However, snail kite movement data are not included since this is a federally protected species under the Endangered Species Act and therefore these data are sensitive. All R scripts can be found in the 'R' folder and all simulated data and model results can be found in the 'data' folder.

To recreate the analyses presented in this paper based on the simulated tracks, R scripts should be run in the following order: 

  1. CRW Mixed Membership Simulation_weird.R
  2. CRW Mixed Membership Simulation_parametric.R
  3. HMM Simulations.R
  4. Segmentation of Behav Sim Tracks.R
  5. LDA Sim Tracks Behav.R
  6. BCPA.R
  7. HMM.R
  8. Segclust2d.R
  9. EMbC.R
  10. Method Comparison_weird.R
  11. Method Comparison_parametric.R
  12. Method Comparison_HMMsims.R

The additional .R files include a file of helper functions for the processing and wrangling of data (helper functions.R), functions used to simulate correlated random walks (CRWs; Simulation Functions.R), and functions used to iterate HMMs, segclust2d models, and EMbC models (Iterate HMMs.R, Iterate Segclust2d.R, Iterate EMbC.R). Other R scripts are included that perform a sensitivity analysis on how the data streams are discretized and compare their performance in estimating breakpoints and the number of likely states (Binning Sensitivity Analysis_*.R). Data files (.csv) are also included that were saved throughout the course of this workflow and were used for data analysis and visualization.
