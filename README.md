# bayesmove_method_comparison
[![DOI](https://zenodo.org/badge/289030089.svg)](https://zenodo.org/badge/latestdoi/289030089)


R code associated with Cullen et al. Methods in Ecology and Evolution paper *Identifying latent behavioral states in animal movement with non-parametric Bayesian methods*.

This repository stores the R scripts used to generate the simulated data, analyze the data with each of the three methods that were compared (BCPA, HMM, non-parametric Bayesian), and generate figures. However, snail kite movement data are not included since this is a federally protected species under the Endangered Species Act and therefore these data are sensitive.

To recreate the analyses presented in this paper based on the simulated tracks, R scripts should be run in the following order: 1) CRW Mixed Membership Simulation_weird.R, 2) Segmentation of Behav Sim Tracks_weird.R, 3) LDA Sim Tracks Behav_weird.R, 4) BCPA.R, 5) HMM_weird.R, and 6) Method Comparison.R. The additional .R files include a file of helper functions for the processing and wrangling of data (helper functions.R), functions used to simulate correlated random walks (CRWs; Simulation Functions.R), and functions used to iterate HMMs over 2-4 possible states for 30 different starting values (Iterate HMMs.R). Data files (.csv) are also included that were saved throughout the course of this workflow and were used for data analysis and visualization.
