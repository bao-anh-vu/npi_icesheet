## Neural posterior inference for calibrating ice sheet simulators"

This folder contains R code for the observing system simulation experiment (OSSE) and the real data example in the manuscript "Neural posterior inference with state-space models for
calibrating ice sheet simulators".

## Folder structure
Each example is stored in one folder. 

Each folder contains scripts that are numbered in the sequence to produce the results of the OSSE/real data example. Each folder also contains a `source` sub-folder, which contains auxiliary R scripts needed for running the 1D Shelfy-Stream Approximation (SSA) ice sheet model.

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
 
