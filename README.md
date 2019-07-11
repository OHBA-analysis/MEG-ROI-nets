
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## ROInets

**A Matlab package for performing leakage-robust network inference between ROIs in MEG data**

The methodology used in this pipeline is set out in 

> Colclough, G. L., Brookes, M., Smith, S. M. and Woolrich, M. W., "A symmetric multivariate leakage correction for MEG connectomes," NeuroImage 117, pp. 439-448 (2015)

This package provides tools for network analysis between regions of interest
in source-reconstructured MEG data, analysing amplitude correlations in 
band-limited power using a Gaussian Markov Random Field network model, after 
applying a symmetric multivariate orthogonalisation correction for source leakage. 

This package was originally developed by [@GilesColclough](https://github.com/GilesColclough), and is now maintained by the [OHBA Analysis](https://github.com/OHBA-analysis) group.

### What you need for it to run

 - FSL
 - FieldTrip
 - Matlab Stats + signal processing toolboxes - though you could probably change the code to make it work
 - QPAS mex files (optional)

### What you need to get started

 - Source-reconstructed resting-state MEG data (you can probably get it to work for task)
 - A set of ROIs or a spatial basis set, in the same space and resolution as the MEG data, saved as a nifti

### How to get started

This package is distributed as part of [OSL](https://ohba-analysis.github.io/osl-docs/); this is how it should be installed on your machine. 

 - The `+ROInets` folder is a Matlab package. Do not change the name of this folder, and do not add its contents on your path. All the functions inside it can be used with a "dot-syntax", by typing e.g. `ROInets.run_network_analysis`.
 - For a brief summary of what each function does, type `help ROInets.Contents`. The help text of each function should provide more information, and examples can be found in `+ROInets/+examples`.
 - The top-level function is `ROInets.run_individual_network_analysis`. View the helptext for this to view all the pipeline options. 
 - If you often use this package, and would like to avoid typing the prefix `ROInets.*` all the time, check out [`import`](https://uk.mathworks.com/help/matlab/ref/import.html).

### What paper to cite

> Colclough, G. L., Brookes, M., Smith, S. M. and Woolrich, M. W., "A symmetric multivariate leakage correction for MEG connectomes, "NeuroImage 117, pp. 439-448 (2015).

### Where to ask for help, or report bugs

For technical support or any other issues, please either [open an issue](https://github.com/OHBA-analysis/MEG-ROI-nets/issues), or contact the OHBA Analysis Group by [email](mailto:mark.woolrich@ohba.ox.ac.uk).

### An overview of the pipeline

 - Select time region for analysis
 - Band-pass filter the data
 - Find a time-course for each ROI using PCA 
 - Remove effects of leakage using an orthogonalisation process which finds the closest set of orthogonal vectors
 - Find down-sampled power envelopes
 - Perform network inference using partial correlation, L1 regularised using DP-glasso, optimised using cross-validation
 - Convert correlations to z-statistics, by scaling relative to an empirical null, generated to share the same temporal smoothness properties as the input data. 

### License

```
Copyright 2015 OHBA
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
```
