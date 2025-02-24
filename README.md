# Replication Code for "Merger Guidelines for the Labor Market"

## Overview

This repository contains the replication code for the paper **"Merger Guidelines for the Labor Market"** by **David Berger, Thomas Hasenzagl, Kyle Herkenhoff, Simon Mongey, and Eric A. Posner**. The code reproduces the results and figures presented in the paper using MATLAB.

## Requirements

To run the replication code, you need the following:

- **Software:** MATLAB (tested on MATLAB R2019b)
- **Operating System:** Tested on macOS and Linux
- **Required Toolboxes:** None

## Repository Structure

```plaintext
bhhmp_replication_code/
│── README.md                           # This file
│── data/                               # Folder for parameter values and input data
│── scripts/                            # Main scripts that replicate figures and tables
│   ├── S01_arnold.m                    # Replicates Table 2
│   ├── S02_REG.m                       # Computes required efficiency gains (REGs)
│   ├── S03_guidelines_table_wages.m    # Replicates Table 3 and D3
│   ├── S04_guidelines_table_welfare.m  # Replicates Table 4
│   ├── S05_merger_simulation_heatmap.m # Replicates Figure 1
│   ├── S06_dwp_heatmap.m               # Replicates Figure 2
│   ├── S07_aggregates.m                # Replicates Figure 4 and Table 5
│   ├── S08_symmetric_mergers_REG.m     # Computes REGs for figure B1
│   ├── S09_symmetric_mergers_figure.m  # Replicates Figure B1
│   ├── S10_REGs_percentile.m           # Replicates Figure C1
│   ├── S11_REGs_shares.m               # Replicates Figure C2
│   ├── S12_merger_fractions.m          # Replicates Figure C3
│   ├── S13_unidimensional_figure.m     # Replicates Figure D1
│   ├── S14_unidimensional_table.m      # Replicates Table D1
│   ├── S15_output_losses.m             # Replicates Figure D2
│   ├── S16_employment_losses.m         # Replicates Figure D3
│   ├── S17_error_rates.m               # Replicates Table D2
│── functions/                          # Contains MATLAB functions used across scripts
│── results/                            # Stores output files, including tables and figures
│   ├── tables/                         # Replicated tables
│   ├── figures/                        # Replicated figures
│   ├── matfiles/                       # Output .mat files from analysis
```
## Notes

Scripts `S03-S07`, `S10-S14`, and `S17` require running `S02_REG.m` first to generate the output file `productivity_gains_results.mat`. Note that this script can be computationally intensive and may require significant memory resources.

Script `S09` requires running `S08_symmetric_mergers_REG.m` first to produce the output file `productivity_gains_results_symmetric.mat`.

To replicate **Table D4**, set `param.alpha` to `1.03` in `S02_REG.m` (you can locate `param.alpha` near the top of the script where parameter initialization occurs) and then run `S03_guidelines_table_wages.m` to generate the table.

To replicate **Table D5**, change the `eta` or `theta` values in `S02_REG.m` (these parameters can be found near the top of the script where parameter initialization occurs). You must run `S02_REG.m` four times for the two different values of `eta` and `theta`. After running these, execute `S03_guidelines_table_wages.m` to generate the corresponding output tables.

For any issues or questions, please contact Thomas Hasenzagl at thomas.hasenzagl@gmail.com or any of the other authors of the paper.
