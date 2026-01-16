## Neural posterior inference for calibrating ice sheet simulators

This folder contains R code for the observing system simulation experiment (OSSE) and the real data example in the manuscript "Neural posterior inference with state-space models for calibrating ice sheet simulators".

## Folder structure
Each example is stored in one folder. 

Each folder contains scripts that are numbered in sequence to produce the results of the OSSE/real data example. Each folder also contains a `source` sub-folder, which contains auxiliary R scripts needed for running the 1D Shelfy-Stream Approximation (SSA) ice sheet model.

For the real data example, there is a `data_preprocessing` folder containing scripts that prepare Antarctic data for the case study. The scripts perform the following steps:
* Read (time-averaged) velocity data and extract a flowline from a selected point in the interior of the ice sheet towards the ocean.
* Read Antarctic (continent-wide) data and filter observations to only those that fall within the Thwaites Glacier boundaries. Datasets used include:
  * Annual surface velocity and surface elevation data from [NASA MEaSUREs ITS_LIVE](https://its-live.jpl.nasa.gov/),
  * BedMap3 data from [SCAR](https://bedmap.scar.org/),
  * BedMachine v3 data from [NSIDC](https://nsidc.org/data/nsidc-0756/versions/3),
  * Surface mass balance data from the [Regional Atmospheric Climate Model 2 (RACMO2), version 2.3](https://zenodo.org/records/3677642).
* Map 2D observations to the 1D flowline by taking the average of the nearest 4 values around each grid point along the flowline as the value at that grid point.
* For the surface velocity and surface elevation, combine annual datasets from 2010--2020 into a consolidated dataset 

## RStudio version requirements
The code in this repository was written using R version 4.4.3. 

## Package requirements 
Running the source code requires the following packages (along with their dependencies, which should be installed automatically):
1. `tensorflow` v2.14
2. `keras` v2.13.0
3. `mvtnorm` v1.2
4. `Matrix` v1.6-5
5. `tidyr` v1.3.1
6. `dplyr` v1.1.4
7. `ggplot2` v3.4.4
8. `gridExtra` v2.3
9. `gtable` v0.3.4  
10. `sf`
11. `matrixStats`
12. `FRK`       
 
