# MATLAB Code for Cross-Feeding Model Analysis

## Overview

This repository contains the MATLAB source code used to generate the modeling results presented in Fig. 2b and Fig. 2d–g of the manuscript:

**Interaction Topology Underpins the Resilience of Microbial Mutualism**

The code implements a reduced model of a two-species cross-feeding community and generates the vector field, perturbation trajectories, landscape representations, and bifurcation diagrams described in the manuscript.

## System Requirements

* MATLAB R2024b
* MATLAB Symbolic Math Toolbox
* Windows 11

The code was tested using MATLAB R2024b on Windows 11. No third-party software packages or non-standard hardware are required.

## Installation

No installation is required beyond MATLAB and the Symbolic Math Toolbox.

Download or clone the repository, place all `.m` files in the same folder, and open that folder as the current directory in MATLAB. Setup should take less than five minutes.

## Reproducing the Figures

### Figure 2b

Run:

```matlab
Fig2b
```

This script generates the two-species vector field, equilibrium points, and separatrix shown in Fig. 2b.

### Figure 2d

Run:

```matlab
Fig2d
```

The perturbation strength is controlled by `mumax_perturbed` near the beginning of `Fig2d.m`.

Run the script three times using the following values:

| Perturbation strength | `mumax_perturbed` |
| --------------------- | ----------------: |
| Weak                  |             `4.2` |
| Moderate              |             `3.0` |
| Strong                |             `1.0` |

The default setting is:

```matlab
mumax_perturbed = 4.2;
```

Therefore, running the unmodified script generates the results for the weak perturbation.

Each run generates the pulse and press population trajectories together with the corresponding normal and perturbed landscape plots.

### Figure 2e–g

Run:

```matlab
Fig2efg
```

This script generates the bifurcation diagrams for metabolite production, dilution rate, and metabolite supplementation shown in Fig. 2e–g.

## Demo and Expected Output

No external input dataset is required. The model parameters and initial conditions are predefined in the scripts and serve as a demonstration of the code.

Running the main scripts opens the corresponding results in MATLAB figure windows. The figures are not saved automatically.

Additional conditions can be explored by modifying the parameters in `getPara.m`, the perturbation strength in `Fig2d.m`, or the parameter ranges in `Fig2efg.m`.

## Runtime

Each main script completed within a few seconds on the tested system. The exact runtime may vary depending on the computer and MATLAB configuration.

## Figure Preparation

The MATLAB scripts reproduce the model-derived content of the figures. Final panel assembly and cosmetic formatting were performed in Adobe Illustrator. These graphical adjustments did not alter the calculated vector field, equilibrium positions, bifurcation results, or the population-axis positions of the landscape state markers.

## Code Availability

The code associated with this study is provided in this repository.
