# Design-Based Treatment Effect Estimation in Paired Cluster-Randomized Experiments

This repository contains the code used to complete the analyses discussed in "A General Framework for Design-Based Treatment Effect Estimation in Paired Cluster-Randomized Experiments."

## Contents:

* `analysis/`
  * `00-helper-functions.R`
    * functions for point and variance estimators in the simulations that are not regression-based estimators
  * `00-sim-helpers.R`
    * functions to run simulations
  * `01-sim-all-est.Rmd`
    * running simulations under all data settings for 20 and 200 pairs and with cluster size correlated with control outcomes or not
    * outputs: `sims_m20_n150/*`, `sims_m200_n150/*`, `sims_m20_n150_neff/*`, `sims_m200_n150_neff/*`
  * `02-cta-eda.Rmd`
    * inputs: `MS_HSdata.Rdata` (from https://github.com/manncz/aeis-aux-rct)
    * outputs: `cta-dat-clean.Rdata`
    * clean TEA data for CTAI study schools and conduct EDA
  * `03-cta-ri.Rmd`
    * inputs: `cta-dat-clean.Rdata`
    * outputs: `cta_ri/*`, `cta_plot/*`
    * runs randomization inference simulations with CTAI data and cleans simulation data for figures
  * `03-figures.Rmd`
    * inputs: `sims_m20_n150/*`,`sims_m200_n150/*`, `sims_m20_n150_neff/*`,`sims_m200_n150_neff/*`, `cta-dat-clean.Rdata`, `cta_plot/*`
    * outputs: all figures
  * `temp/`
    * `sims_m20_n150/`
    * `sims_m200_n150/`
    * `sims_m20_n150_neff/`
    * `sims_m200_n150_neff/`
    * `sims_m200_n150_neff/`
    * `cta_plot/`
    * `cta-dat-clean.Rdata/`

* `figures/`
  * all paper figures

## Notes:

* Simulations were run in the University of Michigan Great Lakes High Performance Computing environment.
* Simulations for 20 pairs were run across 36 cores with .75 GB of memory per core. They took approximately 13 hours or 17-21:43 in CPU walltime.
* Simulations for 200 pairs were were run across 36 cores with 1 GB of memory per core. They took approximately 2 days and 18 hours or 90-23:22 in CPU walltime.
